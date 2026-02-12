import sys
import multiprocessing
import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from molvs import Standardizer
from dimorphite_dl import protonate_smiles

# ---------------------------------------------------------
# Helper Functions (Must be top-level for multiprocessing)
# ---------------------------------------------------------

def check_lipinski(mol):
    """
    Returns True if molecule passes Lipinski's Rule of 5.
    Rules: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
    """
    if mol is None: return False
    
    try:
        # Calculate descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)

        if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
            return True
    except:
        return False
        
    return False

def process_single_molecule(task_data):
    """
    Worker function to process a single molecule based on config.
    Args:
        task_data: Tuple of (mol, name, config_dict)
    Returns:
        List of prepared RDKit Mol objects
    """
    mol, name, config = task_data
    
    # --- Step 1: Lipinski Filter ---
    if config.get('do_filter', True):
        if not check_lipinski(mol):
            return []

    # If "Filter Only" mode is active, return the original molecule immediately
    if config.get('stop_after_filter', False):
        mol.SetProp("_Name", name)
        return [mol]

    # --- Step 2: Standardization ---
    current_mol = mol
    if config.get('do_standardize', True):
        try:
            s = Standardizer()
            mol_std = s.standardize(mol)
            # charge_parent handles neutralization and salt removal
            current_mol = s.charge_parent(mol_std, skip_standardize=True)
        except:
            # If standardization fails, drop molecule
            return []

    # --- Step 3: Protonation (pH 7.0 - 7.4) ---
    protonated_mols = []
    if config.get('do_protonate', True):
        try:
            # isomericSmiles=True is CRITICAL to preserve the stereochemistry 
            # we read from the input 3D structure.
            smi = Chem.MolToSmiles(current_mol, isomericSmiles=True)
            
            # Dimorphite-DL protonation
            protonated_smiles_list = protonate_smiles(smi, ph_min=7.0, ph_max=7.4)
            protonated_mols = [Chem.MolFromSmiles(s) for s in protonated_smiles_list]
        except:
            return []
    else:
        # If skipping protonation, just use the standardized molecule
        protonated_mols = [current_mol]

    # --- Step 4: Stereoisomers & 3D Conformation ---
    prepared_mols = []
    
    for p_mol in protonated_mols:
        if p_mol is None: continue

        # Enumerate Stereoisomers?
        # If the input was 3D and we read chirality, this step will respect it 
        # (only enumerating undefined centers).
        isomers = []
        if config.get('do_stereo', True):
            opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, maxIsomers=4)
            try:
                isomers = list(EnumerateStereoisomers(p_mol, options=opts))
            except:
                continue
        else:
            isomers = [p_mol]
        
        for i, iso in enumerate(isomers):
            # Naming: Append suffix if we generated variants
            variant_suffix = f"_v{i+1}" if len(isomers) > 1 or len(protonated_mols) > 1 else ""
            iso.SetProp("_Name", f"{name}{variant_suffix}")

            # 3D Embedding
            if config.get('do_3d', True):
                iso_h = Chem.AddHs(iso)
                
                # Use ETKDGv3 for better ring conformations
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                
                res = AllChem.EmbedMolecule(iso_h, params)
                
                if res == 0:
                    try:
                        # Energy Minimize using MMFF94
                        AllChem.MMFFOptimizeMolecule(iso_h)
                        prepared_mols.append(iso_h)
                    except:
                        pass # Minimization failed, skip this isomer
            else:
                # Keep 2D (Calculate 2D coords just in case)
                AllChem.Compute2DCoords(iso)
                prepared_mols.append(iso)

    return prepared_mols

# ---------------------------------------------------------
# Main Execution Flow
# ---------------------------------------------------------

def load_molecules(input_file):
    mols = []
    names = []
    print(f"Loading molecules from {input_file}...")
    
    if input_file.endswith('.sdf'):
        # Iterate over the file
        suppl = Chem.SDMolSupplier(input_file, removeHs=False)
        for i, m in enumerate(suppl):
            if m:
                # --- CRITICAL: READ CHIRALITY FROM 3D ---
                try:
                    if m.GetNumConformers() > 0:
                        Chem.AssignStereochemistryFrom3D(m)
                except:
                    pass # Fallback to graph topology if 3D assignment fails
                
                n = m.GetProp("_Name") if m.HasProp("_Name") else f"Ligand_{i}"
                mols.append(m)
                names.append(n)
                
    elif input_file.endswith('.smi') or input_file.endswith('.smiles'):
        with open(input_file, 'r') as f:
            for i, line in enumerate(f):
                if line.strip():
                    parts = line.split()
                    smiles = parts[0]
                    # Read name from second column if available
                    name_in_file = parts[1] if len(parts) > 1 else f"Ligand_{i}"
                    
                    m = Chem.MolFromSmiles(smiles)
                    if m:
                        mols.append(m)
                        names.append(name_in_file)
                        
    return list(zip(mols, names))

def remove_duplicates(mol_list):
    """
    Deduplicates molecules based on Canonical Isomeric SMILES.
    Using Isomeric SMILES ensures enantiomers (R vs S) are treated as different.
    """
    seen = set()
    unique_inputs = []
    
    print("Deduplicating...")
    for mol, name in mol_list:
        try:
            # canonical=True, isomericSmiles=True preserves chirality
            smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if smi not in seen:
                seen.add(smi)
                unique_inputs.append((mol, name))
        except:
            continue
            
    return unique_inputs

def save_mols(mol_list, output_file):
    if not mol_list:
        print("No molecules to save.")
        return

    w = Chem.SDWriter(output_file)
    for m in mol_list:
        if m: w.write(m)
    w.close()
    print(f"Saved {len(mol_list)} molecules to: {output_file}")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Molecule Preparation Pipeline: Filter, Standardize, Protonate, 3D Generate."
    )
    
    # Required Arguments
    parser.add_argument("input", help="Input file path (.sdf or .smi)")
    parser.add_argument("output", help="Output file path (.sdf)")

    # Modes
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--dedup-only", action="store_true", 
                       help="Only load, deduplicate, and save molecules.")
    group.add_argument("--filter-only", action="store_true", 
                       help="Load, deduplicate, apply Lipinski filter, and save (no 3D/protonation).")

    # Step Toggles (Flags to DISABLE specific steps)
    parser.add_argument("--skip-filter", action="store_true", help="Skip Lipinski rule check.")
    parser.add_argument("--skip-standardize", action="store_true", help="Skip MolVS standardization.")
    parser.add_argument("--skip-protonate", action="store_true", help="Skip pH 7.4 protonation.")
    parser.add_argument("--skip-stereo", action="store_true", help="Skip stereoisomer enumeration.")
    parser.add_argument("--skip-3d", action="store_true", help="Skip 3D embedding (save as 2D).")

    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # --- 1. Load Data ---
    raw_data = load_molecules(args.input)
    if not raw_data:
        print("Error: No molecules loaded.")
        sys.exit(1)
    
    count_initial = len(raw_data)
    print(f"-> Initial molecule count: {count_initial}")

    # --- 2. Deduplicate ---
    unique_inputs = remove_duplicates(raw_data)
    
    count_dedup = len(unique_inputs)
    print(f"-> Count after deduplication: {count_dedup}")
    print(f"   (Removed {count_initial - count_dedup} duplicates)")

    # --- Check for "Dedup Only" Mode ---
    if args.dedup_only:
        print("Mode: Deduplicate Only. Saving...")
        save_mols([m for m, n in unique_inputs], args.output)
        return

    # --- 3. Configure Processing ---
    # Build a config dictionary to pass to workers
    config = {
        'do_filter': not args.skip_filter,
        'do_standardize': not args.skip_standardize,
        'do_protonate': not args.skip_protonate,
        'do_stereo': not args.skip_stereo,
        'do_3d': not args.skip_3d,
        'stop_after_filter': args.filter_only
    }

    # Override config for "Filter Only" mode
    if args.filter_only:
        print("Mode: Filter Only (Lipinski). Skipping 3D/Protonation.")
        config['do_filter'] = True
        config['stop_after_filter'] = True
    else:
        print("Mode: Full Pipeline Active")
        print(f"   [x] Filter (Lipinski):  {'ON' if config['do_filter'] else 'OFF'}")
        print(f"   [x] Standardize:        {'ON' if config['do_standardize'] else 'OFF'}")
        print(f"   [x] Protonate (pH 7.4): {'ON' if config['do_protonate'] else 'OFF'}")
        print(f"   [x] Stereo Enumeration: {'ON' if config['do_stereo'] else 'OFF'}")
        print(f"   [x] 3D Generation:      {'ON' if config['do_3d'] else 'OFF'}")

    # --- 4. Prepare Tasks for Multiprocessing ---
    # We pack the config dictionary into the tuple for every molecule
    tasks = [(mol, name, config) for mol, name in unique_inputs]
    
    # Leave one CPU free to keep system responsive
    num_workers = max(1, cpu_count() - 1)
    print(f"\nStarting processing on {num_workers} cores...")

    all_results = []
    
    # Run Pool
    with Pool(processes=num_workers) as pool:
        # imap_unordered is faster as it yields results as soon as they finish
        results_iterator = pool.imap_unordered(process_single_molecule, tasks)
        
        # tqdm wraps the iterator to show the progress bar
        for result_list in tqdm(results_iterator, total=len(tasks), unit="mol", desc="Processing"):
            if result_list:
                all_results.extend(result_list)

    # --- 5. Final Report & Save ---
    count_final = len(all_results)
    print("\nProcessing Complete.")
    print("-" * 30)
    print(f"Initial raw inputs:      {count_initial}")
    print(f"Unique inputs:           {count_dedup}")
    print(f"Final output molecules:  {count_final}")
    print("-" * 30)

    save_mols(all_results, args.output)

if __name__ == "__main__":
    # Essential for Windows multiprocessing
    multiprocessing.freeze_support()
    main()