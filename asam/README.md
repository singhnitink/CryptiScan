# SAM2 & VMD Apptainer Environment

This directory contains the files required to build, configure, and seamlessly run an Apptainer container supporting **SAM2** (ensemble generation), **VMD**, and **NAMD** with GPU acceleration via PyTorch.

## Overview of Files

*   `nodes.sif`: The built Apptainer container containing the entire environment (Apptainer image).
*   `vmd_sam.def`: The Apptainer definition file (recipe) used to build `nodes.sif`. It builds an Ubuntu 24.04-based container with CUDA 12.5.1, PyTorch (CUDA 12.4), SAM2, VMD, NAMD, and related plugins.
*   `notes.dat`: Contains system commands to install Apptainer on an Ubuntu-based system.
*   `nodes_launcher.sh`: A bash script designed to set up a workspace and run SAM2 ensemble generation via the container.

## Setup Instructions

### 1. Install Apptainer
If you don't have Apptainer installed, you can use the instructions from `notes.dat`:
```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

### 2. Build the Container (If needed)
If you need to rebuild `nodes.sif` or modify the container environment, use the definition file:
```bash
sudo apptainer build nodes.sif vmd_sam.def
```

## Usage

### Running SAM2 Ensemble Generation
The provided `nodes_launcher.sh` script makes it easy to start using SAM2. 

1. First, make sure the launcher is executable:
   ```bash
   chmod +x nodes_launcher.sh
   ```

2. Run the script:
   ```bash
   ./nodes_launcher.sh
   ```

### Prerequisites
Before running the script, you must manually create the working directory and place your PDB file there:
1. Create a directory. The name of this directory MUST match the `WORKDIR` variable defined inside the `nodes_launcher.sh` script (the default is `adk_sam`).
2. Place your target PDB file inside that directory, named exactly `protein.pdb` (e.g., `<WORKDIR>/protein.pdb`).

### What `nodes_launcher.sh` does:
1. Ensures the `<WORKDIR>` directory exists.
2. Automatically extracts the baseline config file (`mdcath_model.yaml`) from the container into `<WORKDIR>/`.
3. Executes SAM2 (`generate_ensemble.py`) internally using GPU acceleration (`--nv`), outputting the generated ensemble files to `<WORKDIR>/test_output.*`.

### Advanced Usage & Configuration
If you need to tweak the generation parameters (like temperature, number of frames, or device), open `nodes_launcher.sh` and edit the arguments passed to `generate_ensemble.py`:
*   `-n`: Number of frames
*   `-b`: Batch size
*   `-T`: Temperature
*   `-d`: Device (`cuda` or `cpu`)
