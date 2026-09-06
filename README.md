# CryptiScan

CryptiScan is an automated computational pipeline for predicting, perturbing, and validating **cryptic binding pockets** in proteins.

```mermaid
---
config:
  layout: dagre
  theme: mc
---
flowchart TB
    A["CryptoBank Database"] -- "scrape_cryptobank.py" --> B["Top-N Cryptic Residues"]
    C["RCSB PDB"] -- "Fetch & Clean" --> D["PDB Structure File"]
    B --> E["ESM-Scan (ESM-1v Model)"]
    D --> E
    E -- "Score all 19 substitutions; pick highest score" --> G["MODELLER mutate.py"]
    G --> H["Relaxed Mutant Structure(s)"]
    H --> I["SAM2 / BioEmu Sampling"]
    I --> J["Holo-like Conformations"]
    J --> K["Pocket Analysis (MDpocket + RMSD)"]
    K --> L{"Validation Success?"}
    L -- Yes --> M["Cryptic Target Validated"]
    L -- No --> N["Iterate / Refine"]
```

## Pipeline Overview

1. **CryptoBank Scoring**: Identifies the top residues predicted to form cryptic pockets.
2. **PDB Preparation**: Fetches the experimental structure and cleans chain records.
3. **ESM-Scan (Zero-Shot Scoring)**: Evaluates all 19 alternative amino acid mutations at each cryptic residue using ESM-1v log-probabilities and selects the mutation with the highest evolutionary/fitness score.
4. **MODELLER Mutagenesis**: Constructs mutant 3D structures with sidechain energy minimization and molecular dynamics simulated annealing.
5. **SAM2 Ensemble Sampling**: Generates conformational ensembles on the mutant structures to sample cryptic pocket opening.

For instructions on launching the pipeline, see [pipeline/README.md](file:///home/nsingh/Desktop/github/CryptiScan/pipeline/README.md).
