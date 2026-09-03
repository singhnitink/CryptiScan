We will populate this file for the cryptic pocket algo.

Use the mutation code like this: 
python3 modeller_mutate.py 1dhk 1 HIS A > mutation.log 
% here 1dhk is the pdb filename 1dhk.pdb, 1 is the residue nymber to be mutated, HIS is the new reidue at residue number 1 and A is chain name
```mermaid
---
config:
  layout: elk
---
flowchart TD

    A[CryptoBank Database] -->|scrape_cryptobank.py| B[Top-N Cryptic Residues]

    C[RCSB PDB] -->|Fetch & Clean| D[PDB Structure File]

    B --> E[Choose Mutation Target]
    D --> E

    E -->|Combined Mutant| F[Mutation Specification]

    F --> G[MODELLER mutate.py]

    G --> H["Mutant Structure(s)"]

    H --> I[aSAMt/BioEmu Container]

    I --> J[Holo-like Conformations]

    J --> K[Pocket Analysis]

    K -->|MDpocket + RMSD| L[Validate Holo-like Opening]

    L --> M{Validation Success?}

    M -->|Yes| N[Target Identified]
    M -->|No| E

    classDef input stroke:#818cf8,fill:#eef2ff
    classDef process stroke:#2dd4bf,fill:#f0fdfa
    classDef output stroke:#4ade80,fill:#f0fdf4
    classDef decision stroke:#facc15,fill:#fefce8
    classDef validation stroke:#38bdf8,fill:#f0f9ff

    class A,C input
    class B,D,E,F,G,H,I,J,K process
    class L validation
    class M decision
    class N output
```
