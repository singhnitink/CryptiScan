#!/bin/bash

CONTAINER="nodes.sif"

# Directory where config, input, and output files will be saved


WORKDIR="adk_sam"

mkdir -p $WORKDIR
echo "Please add the pdb file to to $WORKDIR"

# Copying config and input files out of the container
apptainer exec $CONTAINER cp /opt/sam2/config/mdcath_model.yaml $WORKDIR/mdcath_model.yaml


#checks
# echo "Checking sam import"
# apptainer exec $CONTAINER python3 -c "import sam; print('sam ok')"

# echo "Checking PyTorch version"
# apptainer exec $CONTAINER python3 -c "import torch; print('torch version:', torch.__version__)"

# echo "Checking CUDA availability"
# apptainer exec --nv $CONTAINER python3 -c "import torch; print('CUDA available:', torch.cuda.is_available())"


# Run SAM2 ensemble generation

echo "Running SAM2 ensemble generation..."
#change -n: for number of frames; -T: Temperature; -d: cuda or cpu
apptainer exec --nv $CONTAINER python3 /opt/sam2/scripts/generate_ensemble.py \
    -c $WORKDIR/mdcath_model.yaml \
    -i $WORKDIR/protein.pdb \
    -o $WORKDIR/test_output \
    -n 8 \
    -b 4 \
    -T 320 \
    -d cuda

echo "Done! Output files should be at $WORKDIR/test_output.*"