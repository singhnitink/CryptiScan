

We can run the mutation pipeline like this: 

python3 run_pipeline.py --pdb 1JWP --chain A --top 5 --mutate_to GLU --mode combined --scraper_python $CONDA_PREFIX/bin/python3 --modeller_python $CONDA_PREFIX/bin/python3 --check_numbering

Here, PDB ID: 1JWP; selecting top n(5) residues from cryptobank; mutating to (GLU); creating a single mutant with multiple mutations (combined); check any numbering discrepencies in the PDB
