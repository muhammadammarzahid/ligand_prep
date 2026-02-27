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
# Helper Functions
# ---------------------------------------------------------

def check_lipinski(mol):
    if mol is None: return False
    try:
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
    Worker function.
    Collects ALL successful variants first, then names them.
    This ensures we don't add '_v1' if there is only one result.
    """
    mol, header_name, props, config = task_data
    
    # --- Step 1: Lipinski Filter ---
    if config.get('do_filter', True):
        if not check_lipinski(mol):
            return []

    # If "Filter Only", return immediately
    if config.get('stop_after_filter', False):
        return [(mol, header_name, props)]

    # --- Step 2: Standardization ---
    current_mol = mol
    if config.get('do_standardize', True):
        try:
            s = Standardizer()
            mol_std = s.standardize(mol)
            current_mol = s.charge_parent(mol_std, skip_standardize=True)
        except:
            return []

    # --- Step 3: Protonation ---
    protonated_mols = []
    if config.get('do_protonate', True):
        try:
            smi = Chem.MolToSmiles(current_mol, isomericSmiles=True)
            protonated_smiles_list = protonate_smiles(smi, ph_min=7.0, ph_max=7.4)
            protonated_mols = [Chem.MolFromSmiles(s) for s in protonated_smiles_list]
        except:
            return []
    else:
        protonated_mols = [current_mol]

    # --- Step 4: Collect 2D Candidates (Stereo) ---
    # We flatten the list of potential candidates here
    candidates_2d = []
    
    for p_mol in protonated_mols:
        if p_mol is None: continue
        
        if config.get('do_stereo', True):
            opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, maxIsomers=4)
            try:
                isomers = list(EnumerateStereoisomers(p_mol, options=opts))
                candidates_2d.extend(isomers)
            except:
                candidates_2d.append(p_mol)
        else:
            candidates_2d.append(p_mol)

    # --- Step 5: 3D Generation & Filtering ---
    successful_mols = []
    
    for candidate in candidates_2d:
        target_mol = candidate
        
        if config.get('do_3d', True):
            # Try 3D embedding
            iso_h = Chem.AddHs(candidate)
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            
            res = AllChem.EmbedMolecule(iso_h, params)
            if res == 0:
                try:
                    AllChem.MMFFOptimizeMolecule(iso_h)
                    target_mol = iso_h
                    successful_mols.append(target_mol)
                except:
                    # If optimization fails, keeps unoptimized 3D
                    target_mol = iso_h
                    successful_mols.append(target_mol)
            else:
                # Embedding failed, drop this candidate
                continue 
        else:
            # 2D Mode
            AllChem.Compute2DCoords(candidate)
            successful_mols.append(candidate)

    # --- Step 6: Final Naming Logic ---
    # Only append suffix if we actually produced multiple valid molecules
    results_to_return = []
    total_count = len(successful_mols)
    
    for i, final_mol in enumerate(successful_mols):
        if total_count > 1:
            # We have variants, so use suffixes (v1, v2...)
            final_name = f"{header_name}_v{i+1}"
        else:
            # Singleton: Keep original name exactly
            final_name = header_name
            
        results_to_return.append((final_mol, final_name, props))

    return results_to_return

# ---------------------------------------------------------
# Main Execution Flow
# ---------------------------------------------------------

def load_molecules(input_file):
    mols = []
    names = []
    all_props = [] 
    
    print(f"Loading molecules from {input_file}...")
    
    if input_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(input_file, removeHs=False)
        for i, m in enumerate(suppl):
            if m:
                # 1. Capture Properties
                props = {}
                try:
                    for key in m.GetPropNames():
                        props[key] = m.GetProp(key)
                except:
                    pass

                # 2. Capture Name
                n = ""
                if m.HasProp("_Name"):
                    n = m.GetProp("_Name").strip()
                
                if not n:
                    check_list = ["Name", "ID", "compound_id", "Title"]
                    for tag in check_list:
                        for key in props.keys():
                            if tag.lower() == key.lower():
                                n = props[key].strip()
                                break
                        if n: break
                
                if not n:
                    n = f"Ligand_{i+1}"
                
                mols.append(m)
                names.append(n)
                all_props.append(props)
                
    elif input_file.endswith('.smi') or input_file.endswith('.smiles'):
        with open(input_file, 'r') as f:
            for i, line in enumerate(f):
                if line.strip():
                    parts = line.strip().split(None, 1)
                    smiles = parts[0]
                    name_in_file = parts[1] if len(parts) > 1 else f"Ligand_{i+1}"
                    
                    m = Chem.MolFromSmiles(smiles)
                    if m:
                        mols.append(m)
                        names.append(name_in_file)
                        all_props.append({"Original_SMILES": smiles})
    
    # Visual update: Just return the list, main() handles printing counts
    return list(zip(mols, names, all_props))

def remove_duplicates(mol_data):
    seen = set()
    unique_inputs = []
    
    print("Deduplicating...")
    for mol, name, props in mol_data:
        try:
            smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            if smi not in seen:
                seen.add(smi)
                unique_inputs.append((mol, name, props))
        except:
            continue
    return unique_inputs

def main():
    args = parse_arguments() 
    
    # 1. Load (Visual Update: Print Initial Count)
    raw_data = load_molecules(args.input)
    if not raw_data: sys.exit(1)
    
    count_initial = len(raw_data)
    print(f"-> Initial molecule count: {count_initial}")

    # 2. Dedup (Visual Update: Print Deduplication stats)
    unique_inputs = remove_duplicates(raw_data)
    count_dedup = len(unique_inputs)
    
    print(f"-> Count after deduplication: {count_dedup}")
    print(f"   (Removed {count_initial - count_dedup} duplicates)")
    
    # 3. Config
    config = {
        'do_filter': not args.skip_filter,
        'do_standardize': not args.skip_standardize,
        'do_protonate': not args.skip_protonate,
        'do_stereo': not args.skip_stereo,
        'do_3d': not args.skip_3d,
        'stop_after_filter': args.filter_only
    }

    # --- VISUAL ELEMENT MERGE: Print Configuration Block ---
    if args.filter_only:
        print("Mode: Filter Only (Lipinski). Skipping 3D/Protonation.")
    else:
        print("Mode: Full Pipeline Active")
        print(f"   [x] Filter (Lipinski):  {'ON' if config['do_filter'] else 'OFF'}")
        print(f"   [x] Standardize:        {'ON' if config['do_standardize'] else 'OFF'}")
        print(f"   [x] Protonate (pH 7.4): {'ON' if config['do_protonate'] else 'OFF'}")
        print(f"   [x] Stereo Enumeration: {'ON' if config['do_stereo'] else 'OFF'}")
        print(f"   [x] 3D Generation:      {'ON' if config['do_3d'] else 'OFF'}")
    # -------------------------------------------------------

    # 4. Processing
    tasks = [(mol, name, props, config) for mol, name, props in unique_inputs]
    
    if args.use_cores is not None:
        num_workers = max(1, min(args.use_cores, cpu_count()))
    else:
        num_workers = max(1, cpu_count() - 1)
    
    print(f"Processing on {num_workers} cores...")
    
    w = Chem.SDWriter(args.output)
    count_saved = 0

    with Pool(processes=num_workers) as pool:
        results_iterator = pool.imap_unordered(process_single_molecule, tasks)
        
        for result_batch in tqdm(results_iterator, total=len(tasks), unit="mol"):
            if not result_batch: continue
            
            for (mol_obj, final_name, original_props) in result_batch:
                if mol_obj:
                    # 1. Set Header
                    mol_obj.SetProp("_Name", final_name)
                    
                    # 2. Set Data Tags
                    for k, v in original_props.items():
                        if k and str(k).strip(): 
                             mol_obj.SetProp(str(k), str(v))
                    
                    # 3. Ensure Name/ID exists
                    mol_obj.SetProp("Name", final_name)
                    mol_obj.SetProp("ID", final_name)

                    w.write(mol_obj)
                    count_saved += 1

    w.close()
    print(f"\nSuccess. Saved {count_saved} molecules to {args.output}")

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file path (.sdf or .smi)")
    parser.add_argument("output", help="Output file path (.sdf)")
    parser.add_argument("--use_cores", type=int, help="Limit the number of CPU cores used")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--dedup-only", action="store_true")
    group.add_argument("--filter-only", action="store_true")
    parser.add_argument("--skip-filter", action="store_true")
    parser.add_argument("--skip-standardize", action="store_true")
    parser.add_argument("--skip-protonate", action="store_true")
    parser.add_argument("--skip-stereo", action="store_true")
    parser.add_argument("--skip-3d", action="store_true")
    return parser.parse_args()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()