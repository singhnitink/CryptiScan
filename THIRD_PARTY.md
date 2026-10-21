# Third-party components

The MIT licence in `LICENSE` covers the CryptiScan scripts in this repository
only. The container image built from `pipeline_local.def`, and the archived
image deposited with the release, redistribute the components below under
their own licences. Their notices are retained inside the image.

| Component | Version | Licence | Source |
|---|---|---|---|
| aSAM / SAM2 | commit `ce076ed` | Apache 2.0 | https://github.com/giacomo-janson/sam2 |
| fpocket | conda-forge build, pinned in the image | MIT | https://github.com/Discngine/fpocket |
| ESM-1v (`fair-esm`) | 2.0.0, `esm1v_t33_650M_UR90S_1` | MIT | https://github.com/facebookresearch/esm |
| ProtT5-XL-UniRef50 + CryptoBank head | `prot_t5_xl_uniref50_full_v2` | see upstream | https://huggingface.co/ThorbenF/prot_t5_xl_uniref50_full_v2 |
| MDAnalysis | 2.7.0 | GPL-2.0-or-later | https://www.mdanalysis.org |
| MDTraj | 1.10.3 | LGPL-2.1 | https://www.mdtraj.org |
| PyTorch | 2.4.0+cu124 | BSD-3-Clause | https://pytorch.org |
| NumPy, SciPy, pandas, Matplotlib | pinned in the image | BSD-3-Clause / PSF | respective projects |
| OpenMM | conda-forge build, recorded in the image | MIT; CUDA and OpenCL platforms LGPL | https://openmm.org |
| PDBFixer | conda-forge build, recorded in the image | MIT | https://github.com/openmm/pdbfixer |

The exact OpenMM and PDBFixer builds are written to
`/opt/pipeline/MUTATE_VERSION.txt` inside the image.

## Mutagenesis

Earlier versions used MODELLER for the mutagenesis step. Its licence does not
allow redistribution inside the image, so it was replaced by PDBFixer + OpenMM
(`openmm_mutate.py`). The image contains no MODELLER code, and no licence key
is required. `modeller_mutate.py` remains in the repository for reference only.

## MDAnalysis and the GPL

`clean_structure.py` imports MDAnalysis, which is GPL-2.0-or-later. Scientific
Python tools that import MDAnalysis are commonly released under permissive
licences, and that is the practice followed here. Anyone redistributing a
modified version of this repository in a context where GPL compatibility is
material should take their own legal advice.
