# SGRun - Genomic Workflow Automation Tool

This application automates the process of aligning genomic sequences, sorting and indexing alignments, and detecting structural variations. The tool integrates **minimap2**, **samtools**, and **sniffles** into a streamlined workflow, all encapsulated in a Docker container for ease of use.

## Features

- **Sequence Alignment**: Uses `minimap2` for aligning genomic sequences.
- **BAM File Processing**: Sorts and indexes alignment files using `samtools`.
- **Structural Variation Detection**: Identifies structural variations with `sniffles`.
- **Dockerized Workflow**: Runs entirely in a Docker container, ensuring consistent environments across systems.

## Usage

### Docker Requirements

Ensure Docker is installed and running on your system. No additional software installation is required, as all dependencies are included in the container.

### Mount Points

When running the container, the following directories must be mounted:

- `/data/input`: Directory containing input genomic sequence files.
- `/data/output`: Directory where results will be written.
- `/data/workdir`: Working directory for temporary files.

### Run the Workflow

1. Place your input genomic sequence files in a folder in the `/data/input` directory. Note that the reference sequence file *MUST* be the first file alphabetically in it's respective folder otherwise behavior is undefined.
2. Run the Docker container with the required mounts:
   ```bash
   docker run --rm \
     -v /path/to/input:/data/input \
     -v /path/to/output:/data/output \
     -v /path/to/workdir:/data/workdir \
     ghcr.io/xtruan/sgrun:main
