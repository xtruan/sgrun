import os
import glob

def print_fasta_files(folder):
    # Find all .fasta files in the given folder
    fasta_files = glob.glob(os.path.join(folder, '*.fasta'))
    
    # Print the name of each .fasta file
    for file in fasta_files:
        print(os.path.basename(file))

def watch_directories(root_folder):
    # Initialize a set to keep track of directories already processed
    processed_directories = set()
    
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
                    print_fasta_files(os.path.join(dirpath, dirname))
                    
        # Sleep for some time before checking again
        # You can adjust this time according to your needs
        time.sleep(5)

if __name__ == "__main__":
    import time
    
    # Define the root folder to watch
    root_folder = "/input"
    
    # Start watching directories
    watch_directories(root_folder)