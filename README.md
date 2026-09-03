We will populate this file for the cryptic pocket algo.

Use the mutation code like this: 
python3 modeller_mutate.py 1dhk 1 HIS A > mutation.log 
% here 1dhk is the pdb filename 1dhk.pdb, 1 is the residue nymber to be mutated, HIS is the new reidue at residue number 1 and A is chain name
```mermaid
---
config:
  layout: dagre
  theme: mc
---
flowchart TB
    A["CryptoBank Database"] -- "scrape_cryptobank.py" --> B["Top-N Cryptic Residues"]
    C["RCSB PDB"] -- Fetch & Clean --> D["PDB Structure File"]
    B --> E["check the group of the amino acid (eg. hydrophobic) then mutate it with the amino acids from other groups"]
    G["MODELLER mutate.py"] --> H["Mutant Structure(s)"]
    H --> I["aSAMt/BioEmu Container"]
    I --> J["Holo-like Conformations"]
    K["Pocket Analysis"] -- MDpocket + RMSD --> L["Validate Holo-like Opening"]
    L --> M{"Validation Success?"}
    M -- Yes --> N["Target Identified"]
    D --> A
    E --> G

     A
     B
     C
     D
     E
     G
     H
     I
     J
     K
     L
     M
     N
```
