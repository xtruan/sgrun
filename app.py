import os
import glob
import time
import datetime
from flask import Flask, request, abort, send_file, render_template

app = Flask(__name__)

CONDA_ENV = os.getenv('CONDA_ENV', '')
BASE_DIR = '/data/uploads'
INTERACTIVE = os.getenv('PROCESS_INTERACTIVE_MODE', 'True').lower() == 'true'
ONLY_NEW_FILES = os.getenv('PROCESS_ONLY_NEW_FILES', 'False').lower() == 'true'

gQueue = []

def posix_sec_now():
    return time.time()

def posix_sec_to_datetime(posix_sec):
    return datetime.datetime.fromtimestamp(float(posix_sec)).replace(tzinfo=datetime.timezone.utc)

def datetime_now():
    return posix_sec_to_datetime(posix_sec_now())

def datetime_now_label():
    return str(datetime_now()).split('+')[0].replace(' ', '_').replace('-', '_').replace(':', '_').replace('.', '_').replace('+', '_')

@app.route("/upload", methods=["POST"])
def upload_file():
    # Get the file object
    files = request.files.getlist("file")

    dir = os.path.join(BASE_DIR, datetime_now_label())
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

    return "File uploaded successfully! Queue: " + str(gQueue)

@app.route("/process")
def process():
    output_name = "alignment"
    
    global gQueue
    if len(gQueue) == 0:
        return "Queue is empty!"
    
    files = gQueue.pop(0)
    if len(files) > 0:
        file_dir = os.path.dirname(files[0])
        files_str = " ".join(files)
        # Run minimap2
        os.system(f"minimap2 -ax asm5 {files_str} > {file_dir}/{output_name}.sam")
        #file_dir = os.path.join(BASE_DIR, "2023_07_21_22_49_02_746364")
        # Run samtools
        os.system(f"samtools sort {file_dir}/{output_name}.sam > {file_dir}/{output_name}.bam")
        os.system(f"samtools index {file_dir}/{output_name}.bam")
        # Run sniffles
        os.system(f"sniffles -i {file_dir}/{output_name}.bam -v {file_dir}/{output_name}.vcf")

        return "Processed successfully!"
    
    return "No files to process!"

    # 1) detect files added and run pre-processing sam/bam
    # 2) spins up in response to web call for analysis vcf

@app.route('/', defaults={'req_path': ''})
@app.route('/<path:req_path>')
def dir_listing(req_path):

    # Joining the base and the requested path
    abs_path = os.path.join(BASE_DIR, req_path)

    # Return 404 if path doesn't exist
    if not os.path.exists(abs_path):
        return abort(404)

    # Check if path is a file and serve
    if os.path.isfile(abs_path):
        return send_file(abs_path)

    # Show directory contents
    files = os.listdir(abs_path)
    return render_template('files.html', files=files)

def enqueue_fasta_files(folder):
    # Find all .fasta files in the given folder
    fasta_files = glob.glob(os.path.join(folder, '*.fasta'))
    
    # Print the name of each .fasta file
    for file in fasta_files:
        print("Adding " + os.path.basename(file))
    
    # Add the fasta files to the queue
    if len(fasta_files) > 0:
        global gQueue
        gQueue.append(fasta_files)

def watch_directories(root_folder):
    # Initialize a set to keep track of directories already processed
    processed_directories = set()

    do_process = True
    if ONLY_NEW_FILES:
        # Skip processing on the first run if we only want to process new files
        do_process = False

    
    while True:
        # Walk through the root folder
        for dirpath, dirnames, _ in os.walk(root_folder):
            for dirname in dirnames:
                # Check if the directory has not been processed yet
                if dirpath + dirname not in processed_directories:
                    # Add the directory to the set of processed directories
                    processed_directories.add(dirpath + dirname)
                    
                    # Print the fasta files in this directory
                    print(f"Fasta files in directory '{dirname}':")
                    enqueue_fasta_files(os.path.join(dirpath, dirname))
                    
        global gQueue
        # Process the enqueued files
        if do_process:
            print("Processing enqueued files...")
            print(str(gQueue))
            print(process())
        else:
            # Clear the queue if we are not processing files
            gQueue = []
            # Set the flag to process files on the next iteration
            do_process = True

        # Sleep for some time before checking again
        # You can adjust this time according to your needs
        time.sleep(5)

if __name__ == "__main__":

    os.makedirs(BASE_DIR, exist_ok=True)

    if INTERACTIVE:
        print("Running in interactive mode...")
        envPort = int(os.getenv('SERVER_PORT', '8050'))
        envDebug = os.getenv('DEBUG_MODE', 'True').lower() == 'true'
        app.run(debug=envDebug, host='0.0.0.0', port=envPort)

    else:
        watch_directories(BASE_DIR)