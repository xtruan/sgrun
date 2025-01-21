import os
import glob
import shutil
import json
import time
import datetime
from flask import Flask, request, abort, send_file, render_template

app = Flask(__name__)

CONDA_ENV = os.getenv('CONDA_ENV', '')

# INPUT_DIR and OUTPUT_DIR should be external storage mounted to container
# WORK_DIR is an internal volume of the container
BASE_DIR = '/data'
INPUT_DIR = BASE_DIR + '/input'
INDICATOR_FILE_NAME = '.sgrun'
WORK_DIR = BASE_DIR + '/workdir/sgtmp'
OUTPUT_DIR = BASE_DIR + '/output'
OUTPUT_FILE_BASENAME = 'alignment'

INTERACTIVE = os.getenv('PROCESS_INTERACTIVE_MODE', 'True').lower() == 'true'
ONLY_NEW_FILES = os.getenv('PROCESS_ONLY_NEW_FILES', 'True').lower() == 'true'

gQueue = []

def posix_sec_now():
    return time.time()

def posix_sec_to_datetime(posix_sec):
    return datetime.datetime.fromtimestamp(float(posix_sec)).replace(tzinfo=datetime.timezone.utc)

def datetime_now():
    return posix_sec_to_datetime(posix_sec_now())

def datetime_now_label():
    return str(datetime_now()).split('+')[0].replace(' ', '_').replace('-', '_').replace(':', '_').replace('.', '_').replace('+', '_')

# only used in INTERACTIVE mode
@app.route('/', defaults={'req_path': ''})
@app.route('/<path:req_path>')
def dir_listing(req_path):

    # Joining the base and the requested path
    abs_path = os.path.join(WORK_DIR, req_path)

    # Return 404 if path doesn't exist
    if not os.path.exists(abs_path):
        return abort(404)

    # Check if path is a file and serve
    if os.path.isfile(abs_path):
        return send_file(abs_path)

    # Show directory contents
    files = os.listdir(abs_path)
    return render_template('files.html', files=files)

# only used in INTERACTIVE mode
@app.route("/upload", methods=["POST"])
def upload_file():
    ret = {}
    # Get the file object
    files = request.files.getlist("file")

    dir = os.path.join(WORK_DIR, datetime_now_label())
    os.makedirs(dir, exist_ok=True)

    uploaded_files = []
    for file in files:
        # Save the file to disk
        uploaded_file = os.path.join(dir, file.filename)
        with open(uploaded_file, "wb") as f:
            f.write(file.read())
            uploaded_files.append(uploaded_file)

    global gQueue
    gQueue.append(uploaded_files)
    print(gQueue)

    ret["msg"] = "File uploaded successfully!"
    ret["queue"] = str(gQueue)
    return json.dumps(ret)

# used in both INTERACTIVE and headless mode
@app.route("/process")
def process():
    ret = {}
    output_name = OUTPUT_FILE_BASENAME
    
    global gQueue
    if len(gQueue) == 0:
        ret["msg"] = "Queue is empty!"
        return json.dumps(ret)
    
    files = gQueue.pop(0)
    if len(files) > 0:
        file_dir = os.path.dirname(files[0])
        files_str = " ".join(files)
        output_name = output_name + '_' + os.path.basename(files[0])
        # Run minimap2
        os.system(f"minimap2 -ax asm5 {files_str} > {file_dir}/{output_name}.sam")
        # Run samtools
        os.system(f"samtools sort {file_dir}/{output_name}.sam > {file_dir}/{output_name}.bam")
        os.system(f"samtools index {file_dir}/{output_name}.bam")
        # Run sniffles
        os.system(f"sniffles -i {file_dir}/{output_name}.bam -v {file_dir}/{output_name}.vcf")

        ret["msg"] = "Processed successfully!"
        ret["dir"] = file_dir
        return json.dumps(ret)
    
    ret["msg"] = "No files to process!"
    return json.dumps(ret)

    # 1) detect files added and run pre-processing sam/bam
    # 2) spins up in response to web call for analysis vcf

def copy_fasta_files(input_folder, fasta_folder, output_folder, move=False):
    # Find all .fasta files in the given folder
    return copy_files('*.fasta', input_folder, fasta_folder, output_folder, move)

def copy_files(glob_pattern, input_folder, search_folder, output_folder, move=False):
    # Find all matching files in the given folder
    files = glob.glob(os.path.join(search_folder, glob_pattern))
    
    # Move/Copy each matching file to the output folder while maintaining the directory structure
    destination_dir = None
    for file in files:
        # Extract the relative path from the root input folder to the file
        relative_path = os.path.relpath(file, input_folder)
        
        # Construct the destination path in the output folder
        destination_path = os.path.join(output_folder, relative_path)
        destination_path = destination_path.replace(" ", "_")
        
        # Ensure the directory structure exists in the output folder
        os.makedirs(os.path.dirname(destination_path), exist_ok=True)
        
        if move:
            # Move the file to the output folder
            print(f"Moving file {file} to {destination_path}...")
            shutil.move(file, destination_path)
        else:
            # Copy the file to the output folder
            print(f"Copying file {file} to {destination_path}...")
            shutil.copy(file, destination_path)

        destination_dir = os.path.dirname(destination_path)

    return destination_dir

def enqueue_fasta_files(fasta_folder):
    # Find all .fasta files in the given folder
    fasta_files = glob.glob(os.path.join(fasta_folder, '*.fasta'))
    
    # Print the name of each .fasta file
    for file in fasta_files:
        print("Adding " + os.path.basename(file))
    
    # Add the fasta files to the queue
    if len(fasta_files) > 0:
        global gQueue
        gQueue.append(fasta_files)

def watch_directories(input_folder, working_folder, output_folder):
    do_process = True
    
    while True:
        # Walk through the root folder
        break_all = False
        for dirpath, dirnames, _ in os.walk(input_folder):
            for dirname in dirnames:
                # Check if the directory has not been processed yet
                processed_file = os.path.join(dirpath, dirname, INDICATOR_FILE_NAME)
                #print(f"new check: {(not ONLY_NEW_FILES)} or {(not os.path.exists(processed_file))}")
                if (not ONLY_NEW_FILES) or (not os.path.exists(processed_file)):
                    # Add the directory to the set of processed directories using indicator file
                    with open(processed_file, 'w') as f:
                        f.write('start\n')
                    
                    # Clean working directory
                    shutil.rmtree(working_folder)
                    os.makedirs(working_folder, exist_ok=True)
                    # Copy fasta files to working directory
                    destination_folder = copy_fasta_files(
                        input_folder, 
                        os.path.join(dirpath, dirname), 
                        working_folder)

                    # Print the fasta files in this directory
                    if destination_folder is not None:
                        with open(processed_file, 'a') as f:
                            f.write('copied inputs\n')

                        print(f"Fasta files in directory '{destination_folder}':")
                        enqueue_fasta_files(destination_folder)
                        with open(processed_file, 'a') as f:
                            f.write('enqueued\n')
                        
                    global gQueue
                    # Process the enqueued files in the current directory
                    if do_process:
                        print("Processing enqueued files...")
                        print(str(gQueue))
                        empty_queue = (not gQueue)
                        out = json.loads(process())
                        with open(processed_file, 'a') as f:
                            f.write('processed\n')
                        if 'msg' in out:
                            print(out['msg'])
                        if 'dir' in out:
                            print(out['dir'])
                            destination_folder = copy_files(
                                OUTPUT_FILE_BASENAME + '*', 
                                working_folder, 
                                out['dir'], 
                                output_folder)
                            if destination_folder is not None:
                                with open(processed_file, 'a') as f:
                                    f.write('copied outputs\n')
                        # Break out of the loops so we explicity check for new files after processing
                        # Only do this if the queue was not empty
                        print("Processing complete!")
                        if not empty_queue:
                            print("Restart loop - process complete")
                            break_all = True
                    else:
                        # Clear the queue if we are not processing files
                        gQueue = []
                        # Set the flag to process files on the next iteration
                        do_process = True

                    with open(processed_file, 'a') as f:
                        f.write('done\n')

                else:
                    print(f"{INDICATOR_FILE_NAME} file found: '{processed_file}'")

                if break_all:
                    # Break inner for
                    break
            if break_all:
                break_all = False
                # Break outer for
                break

        # Sleep for some time before checking again
        # You can adjust this time according to your needs
        print(f"Sleeping...")
        time.sleep(20)

if __name__ == "__main__":

    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(WORK_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if INTERACTIVE:
        print("Running in interactive mode...")
        envPort = int(os.getenv('SERVER_PORT', '8050'))
        envDebug = os.getenv('DEBUG_MODE', 'True').lower() == 'true'
        app.run(debug=envDebug, host='0.0.0.0', port=envPort)

    else:
        watch_directories(INPUT_DIR, WORK_DIR, OUTPUT_DIR)