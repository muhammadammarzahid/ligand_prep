# LigandPrep Pipeline

**LigandPrep Pipeline** is a high-performance, multi-processed Python command-line tool designed to prepare chemical libraries for virtual screening, molecular docking, and cheminformatics workflows.

It automates the cleaning of raw molecule data, filters by drug-likeness, adjusts protonation states for biological pH, and generates energy-minimized 3D conformers while respecting input stereochemistry.

## 🚀 Features

* **Multiprocessing Support:** Automatically utilizes available CPU cores (leaving one free) for maximum efficiency.
* **Smart Deduplication:** Removes duplicates based on Canonical Isomeric SMILES, ensuring enantiomers are treated as distinct.
* **3D Chirality Awareness:** Detects if input files (SDF) have 3D coordinates and preserves the specific stereochemistry (R/S) during processing.
* **Lipinski Filtering:** Optional filtering based on the Rule of 5 (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10).
* **Standardization:** Uses `MolVS` to standardize tautomers and neutralize charges.
* **Biological Protonation:** Uses `Dimorphite-DL` to generate protonation states relevant to pH 7.0–7.4.
* **Stereoisomer Enumeration:** Generates isomers for undefined chiral centers (max 4 per molecule).
* **3D Conformation:** Generates 3D coordinates using ETKDGv3 and performs MMFF energy minimization.
* **Real-time Progress:** Features a `tqdm` progress bar and detailed molecule counting stats at every stage.

---

## 🛠️ Installation

Because this script relies on **RDKit**, it is highly recommended to use **Conda** or **Mamba** for installation.

### 1. Create a Conda Environment

Create a clean environment with Python and RDKit installed:

```bash
conda create -n ligandprep python=3.9 -c conda-forge
conda activate ligandprep
conda install -c conda-forge rdkit

```

### 2. Install Python Dependencies

Install the remaining required packages via pip:

```bash
pip install molvs dimorphite-dl tqdm

```

### 3. Setup

Save the python script as `ligand_prep.py` in your working directory.

---

## 📖 Usage

### Basic Command

The simplest usage takes an input file and an output filename. This runs the **full pipeline** (Filter → Standardize → Protonate → Stereo → 3D).

```bash
python ligand_prep.py inputs.smi output.sdf

```

* **Input formats:** `.smi`, `.smiles`, `.sdf`
* **Output format:** `.sdf` (Contains 3D coordinates and properties)

### Modes & Arguments

You can customize the pipeline using flags.

| Flag                   | Description                                                                    |
| ---------------------- | ------------------------------------------------------------------------------ |
| `--dedup-only`       | Loads, removes duplicates, and saves immediately. No processing.               |
| `--filter-only`      | Loads, removes duplicates, checks Lipinski rules, and saves. No 3D generation. |
| `--skip-filter`      | Disables the Lipinski Rule of 5 filter (keeps all molecules).                  |
| `--skip-standardize` | Skips MolVS standardization (tautomers/charge correction).                     |
| `--skip-protonate`   | Skips pH 7.4 protonation (keeps original protonation state).                   |
| `--skip-stereo`      | Skips enumeration of stereoisomers.                                            |
| `--skip-3d`          | Skips 3D embedding and minimization (saves as 2D).                             |

---

## 💡 Usage Examples

### 1. The "Docking Ready" Run

Process a raw SMILES list, clean it, and generate 3D structures for docking.

```bash
python ligand_prep.py raw_compounds.smi prepared_ligands.sdf

```

### 2. Large Library Filtering

If you have a massive library (e.g., 100k molecules) and just want to remove heavy/non-drug-like molecules before doing expensive 3D calculations:

```bash
python ligand_prep.py huge_library.sdf filtered_library.sdf --filter-only

```

### 3. Just Deduplication

Clean up a file by removing duplicate structures.

```bash
python ligand_prep.py messy_data.smi clean_unique.sdf --dedup-only

```

### 4. Custom Pipeline (No Protonation)

If your molecules are already protonated correctly (e.g. from ZINC15), but you need 3D coordinates:

```bash
python ligand_prep.py inputs.sdf output.sdf --skip-protonate

```

---

## 📊 Pipeline Logic

1. **Loader:** Reads `.sdf` or `.smi` files. Checks for 3D coordinates to assign initial stereochemistry tags.
2. **Deduplicator:** Hashes molecules to Canonical Isomeric SMILES to remove exact duplicates.
3. **Lipinski Filter:** Checks MW, LogP, HBA, HBD. (Configurable).
4. **Standardizer:** Applies `MolVS` standardization to fix tautomers and salts.
5. **Protonator:** Uses `Dimorphite-DL` to adjust hydrogens for pH 7.4.
6. **Stereo Enumerator:** Expands *undefined* chiral centers (up to 4 isomers). Centers defined in step 1 are preserved.
7. **3D Embedder:** Adds Hydrogens -> Embeds (ETKDGv3) -> MMFF94 Minimization.
8. **Writer:** Saves to SDF.

## ⚠️ Notes

* **Windows Users:** The script includes `multiprocessing.freeze_support()`, so it works natively on Windows. However, ensure you run it inside a standard terminal (CMD/PowerShell/Anaconda Prompt).
* **Performance:** The script reserves 1 CPU core for the OS and uses the rest for calculation. 3D embedding is the most time-consuming step.
* **Failures:** If a molecule fails embedding or sanitization, it is silently dropped to prevent the pipeline from crashing. The final count will reflect this.
