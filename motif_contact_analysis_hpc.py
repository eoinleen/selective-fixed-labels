#!/usr/bin/env python3
"""
Motif Contact Analysis for HPC Systems
Analyzes PDB files with two chains and identifies contact residues for fixing.
Outputs JSON format compatible with auto_config_fixed_script.py

This version is adapted for HPC systems (removed Colab dependencies)
and outputs structured data for automated script generation.

Usage:
python motif_contact_analysis_hpc.py --pdb input.pdb --output contacts.json
python motif_contact_analysis_hpc.py --pdb input.pdb --distance_cutoffs 3.5,4.0,4.5 --verbose
"""

# ============================================================================
# EASY CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# DEFAULT CHAIN ASSIGNMENTS
# Set these if you want to avoid specifying --motif_chain and --target_chain every time
DEFAULT_MOTIF_CHAIN = "A"      # Chain containing the motif residues
DEFAULT_TARGET_CHAIN = "B"     # Chain containing the target protein

# DEFAULT RESIDUE RANGES  
# Set these if you want to avoid specifying --motif_residues every time
# Examples:
#   Single motif: "66-93"
#   Dual motifs:  "66-93,109-127" 
#   Individual:   "70,75,80-85"
#   All residues: None or "all"
DEFAULT_MOTIF_RESIDUES = "66-93,109-127"   # ← CHANGE THIS FOR YOUR MOTIF
DEFAULT_TARGET_RESIDUES = "all"            # Usually "all" unless focusing on specific domain

# DEFAULT DISTANCE CUTOFFS
DEFAULT_DISTANCE_CUTOFFS = "3.0,3.5,4.0,4.5,5.0"

# ============================================================================
# COMMAND LINE EXAMPLES:
# 
# Using defaults (set above):
#   python motif_contact_analysis_hpc.py --pdb your_motif.pdb
#
# Override motif regions:
#   python motif_contact_analysis_hpc.py --pdb your_motif.pdb --motif_residues "66-93"
#
# Full specification:
#   python motif_contact_analysis_hpc.py --pdb your_motif.pdb \
#     --motif_chain A --motif_residues "66-93,109-127" \
#     --target_chain B --target_residues "all" --verbose
# ============================================================================

import json
import argparse
import math
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze motif contacts and output JSON for auto-config',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Configuration defaults (can be changed at top of script):
  Motif chain:     {DEFAULT_MOTIF_CHAIN}
  Target chain:    {DEFAULT_TARGET_CHAIN}
  Motif residues:  {DEFAULT_MOTIF_RESIDUES}
  Target residues: {DEFAULT_TARGET_RESIDUES}
  Distance cutoffs: {DEFAULT_DISTANCE_CUTOFFS}

Examples:
  # Use all defaults
  python {__file__} --pdb motif.pdb
  
  # Override motif region
  python {__file__} --pdb motif.pdb --motif_residues "66-93"
  
  # Full custom specification  
  python {__file__} --pdb motif.pdb --motif_chain A --motif_residues "66-93,109-127" --target_chain B
        """)
    
    parser.add_argument('--pdb', type=str, required=True, help='Input PDB file')
    parser.add_argument('--output', type=str, help='Output JSON file (default: auto-generated)')
    parser.add_argument('--distance_cutoffs', type=str, default=DEFAULT_DISTANCE_CUTOFFS, 
                       help=f'Distance cutoffs (default: {DEFAULT_DISTANCE_CUTOFFS})')
    parser.add_argument('--motif_chain', type=str, default=DEFAULT_MOTIF_CHAIN,
                       help=f'Motif chain ID (default: {DEFAULT_MOTIF_CHAIN})')
    parser.add_argument('--target_chain', type=str, default=DEFAULT_TARGET_CHAIN,
                       help=f'Target chain ID (default: {DEFAULT_TARGET_CHAIN})')
    parser.add_argument('--motif_residues', type=str, default=DEFAULT_MOTIF_RESIDUES,
                       help=f'Motif residue ranges (default: {DEFAULT_MOTIF_RESIDUES})')
    parser.add_argument('--target_residues', type=str, default=DEFAULT_TARGET_RESIDUES,
                       help=f'Target residue ranges (default: {DEFAULT_TARGET_RESIDUES})')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    return parser.parse_args()

def parse_residue_ranges(range_string):
    """
    Parse residue range string into list of residue numbers
    
    Examples:
    "66-93" -> [66, 67, 68, ..., 93]
    "66-93,109-127" -> [66, 67, ..., 93, 109, 110, ..., 127]
    "70,75,80-85" -> [70, 75, 80, 81, 82, 83, 84, 85]
    """
    if not range_string or range_string.lower() == 'all':
        return None  # Return None to indicate all residues
    
    residues = []
    
    for part in range_string.split(','):
        part = part.strip()
        if '-' in part:
            # Range format (e.g., "66-93")
            start, end = map(int, part.split('-'))
            residues.extend(range(start, end + 1))
        else:
            # Single residue
            residues.append(int(part))
    
    return sorted(list(set(residues)))

def filter_atoms_by_residues(atoms, allowed_residues):
    """
    Filter atoms to only include specified residues
    If allowed_residues is None, return all atoms
    """
    if allowed_residues is None:
        return atoms
    
    return [atom for atom in atoms if atom['res_num'] in allowed_residues]

def parse_pdb_simple(pdb_file):
    """
    Simple PDB parser that extracts coordinates and residue info
    Returns dict with chain information
    """
    chains = defaultdict(list)

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    atom_info = {
                        'res_num': res_num,
                        'res_name': res_name,
                        'atom_name': atom_name,
                        'coords': [x, y, z],
                        'chain': chain_id
                    }
                    chains[chain_id].append(atom_info)
    except Exception as e:
        raise ValueError(f"Error parsing PDB file: {e}")

    return chains
    """
    Simple PDB parser that extracts coordinates and residue info
    Returns dict with chain information
    """
    chains = defaultdict(list)

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    atom_info = {
                        'res_num': res_num,
                        'res_name': res_name,
                        'atom_name': atom_name,
                        'coords': [x, y, z],
                        'chain': chain_id
                    }
                    chains[chain_id].append(atom_info)
    except Exception as e:
        raise ValueError(f"Error parsing PDB file: {e}")

    return chains

def identify_chains(chains, motif_chain=None, target_chain=None, motif_residues=None, target_residues=None):
    """
    Identify motif (smaller) and target (larger) chains
    """
    chain_sizes = {chain: len(set(atom['res_num'] for atom in atoms))
                  for chain, atoms in chains.items()}

    print(f"Chain sizes: {chain_sizes}")

    if motif_chain and target_chain:
        print(f"Using specified chains: motif={motif_chain}, target={target_chain}")
    else:
        motif_chain = min(chain_sizes.keys(), key=lambda x: chain_sizes[x])
        target_chain = max(chain_sizes.keys(), key=lambda x: chain_sizes[x])

    # Parse residue ranges
    motif_res_list = parse_residue_ranges(motif_residues)
    target_res_list = parse_residue_ranges(target_residues)
    
    # Report residue filtering
    if motif_res_list:
        print(f"Motif chain: {motif_chain} (analyzing residues: {min(motif_res_list)}-{max(motif_res_list)}, {len(motif_res_list)} total)")
    else:
        print(f"Motif chain: {motif_chain} (analyzing all {chain_sizes[motif_chain]} residues)")
        
    if target_res_list:
        print(f"Target chain: {target_chain} (analyzing residues: {min(target_res_list)}-{max(target_res_list)}, {len(target_res_list)} total)")
    else:
        print(f"Target chain: {target_chain} (analyzing all {chain_sizes[target_chain]} residues)")

    return motif_chain, target_chain, motif_res_list, target_res_list

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates"""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def calculate_distances(motif_atoms, target_atoms, heavy_atoms_only=True):
    """
    Calculate minimum distances between all motif and target atoms
    Returns dict with residue-level minimum distances
    """
    motif_residues = defaultdict(list)
    target_coords = []

    # Group motif atoms by residue
    for atom in motif_atoms:
        if heavy_atoms_only and atom['atom_name'].startswith('H'):
            continue
        motif_residues[atom['res_num']].append(atom)

    # Get all target coordinates
    for atom in target_atoms:
        if heavy_atoms_only and atom['atom_name'].startswith('H'):
            continue
        target_coords.append(atom['coords'])

    if not target_coords:
        raise ValueError("No target coordinates found")

    # Calculate minimum distance for each motif residue
    residue_distances = {}

    for res_num, atoms in motif_residues.items():
        min_distance = float('inf')
        res_name = atoms[0]['res_name']

        for atom in atoms:
            for target_coord in target_coords:
                distance = calculate_distance(atom['coords'], target_coord)
                min_distance = min(min_distance, distance)

        residue_distances[res_num] = {
            'min_distance': min_distance,
            'res_name': res_name
        }

    return residue_distances

def analyze_contacts(pdb_file, distance_cutoffs, motif_chain=None, target_chain=None, 
                    motif_residues=None, target_residues=None, verbose=False):
    """
    Main analysis function with residue range filtering
    """
    if verbose:
        print(f"Analyzing PDB file: {pdb_file}")
        print("="*50)

    # Parse PDB
    chains = parse_pdb_simple(pdb_file)

    if len(chains) < 2:
        raise ValueError(f"Expected at least 2 chains, found {len(chains)}: {list(chains.keys())}")

    # Identify motif and target chains with residue filtering
    motif_chain, target_chain, motif_res_list, target_res_list = identify_chains(
        chains, motif_chain, target_chain, motif_residues, target_residues
    )

    if motif_chain not in chains or target_chain not in chains:
        raise ValueError(f"Specified chains not found in PDB")

    # Filter atoms by residue ranges
    motif_atoms = filter_atoms_by_residues(chains[motif_chain], motif_res_list)
    target_atoms = filter_atoms_by_residues(chains[target_chain], target_res_list)
    
    if not motif_atoms:
        raise ValueError(f"No motif atoms found in specified range")
    if not target_atoms:
        raise ValueError(f"No target atoms found in specified range")
    
    if verbose:
        motif_res_count = len(set(atom['res_num'] for atom in motif_atoms))
        target_res_count = len(set(atom['res_num'] for atom in target_atoms))
        print(f"Filtered to {motif_res_count} motif residues and {target_res_count} target residues")

    # Calculate distances
    if verbose:
        print("\nCalculating distances...")
    
    distances = calculate_distances(motif_atoms, target_atoms)

    # Analyze contacts at different cutoffs
    results = {}

    if verbose:
        print(f"\nContact Analysis for Chain {motif_chain}:")
        print("="*40)

    for cutoff in distance_cutoffs:
        contacts = {res: info for res, info in distances.items()
                   if info['min_distance'] <= cutoff}

        contact_list = [(res, info['res_name'], info['min_distance'])
                       for res, info in contacts.items()]
        contact_list.sort()

        results[cutoff] = contact_list

        if verbose:
            print(f"\nContacts within {cutoff}Å: {len(contact_list)} residues")
            for res, res_name, dist in contact_list:
                print(f"  {motif_chain}:{res}:{res_name} (min dist: {dist:.2f}Å)")

    return results, motif_chain, target_chain, distances

def generate_contact_json(results, motif_chain, target_chain, pdb_file, all_distances):
    """
    Generate JSON output compatible with auto_config_fixed_script.py
    """
    
    # Get all residues that appear in any cutoff
    all_residues = set()
    for contact_list in results.values():
        for res, _, _ in contact_list:
            all_residues.add(res)

    # Create residue details
    residue_details = {}
    for res in sorted(all_residues):
        res_name = None
        distances_dict = {}
        
        # Find residue name and distances at each cutoff
        for cutoff, contact_list in results.items():
            found = False
            for r, rn, d in contact_list:
                if r == res:
                    distances_dict[f"{cutoff}A"] = round(d, 2)
                    res_name = rn
                    found = True
                    break
            if not found:
                distances_dict[f"{cutoff}A"] = None
        
        residue_details[res] = {
            'residue_name': res_name,
            'distances': distances_dict,
            'min_distance': min([d for d in distances_dict.values() if d is not None])
        }

    # Generate contact lists for each cutoff
    contact_lists = {}
    for cutoff, contact_list in results.items():
        contact_lists[f"{cutoff}A"] = [res for res, _, _ in contact_list]

    # Identify structural residues
    structural_residues = []
    for res, details in residue_details.items():
        if details['residue_name'] in ['PRO', 'CYS']:
            structural_residues.append(res)

    # Create final JSON structure
    contact_data = {
        'analysis_info': {
            'pdb_file': pdb_file,
            'motif_chain': motif_chain,
            'target_chain': target_chain,
            'total_motif_residues': len(all_distances),
            'analysis_method': 'minimum_heavy_atom_distance'
        },
        'distance_cutoffs': list(results.keys()),
        'contact_lists': contact_lists,
        'residue_details': residue_details,
        'structural_residues': structural_residues,
        'recommendations': {
            'conservative_3.5A': contact_lists.get('3.5A', []),
            'liberal_4.0A': contact_lists.get('4.0A', []),
            'with_structural': sorted(set(contact_lists.get('3.5A', []) + structural_residues))
        }
    }

    return contact_data

def print_summary(contact_data, verbose=False):
    """Print human-readable summary"""
    
    info = contact_data['analysis_info']
    print("\n" + "="*60)
    print("MOTIF CONTACT ANALYSIS SUMMARY")
    print("="*60)
    print(f"PDB File: {info['pdb_file']}")
    print(f"Motif Chain: {info['motif_chain']} ({info['total_motif_residues']} residues)")
    print(f"Target Chain: {info['target_chain']}")
    
    print(f"\nContact counts by distance cutoff:")
    for cutoff_key, residues in contact_data['contact_lists'].items():
        print(f"  {cutoff_key}: {len(residues)} residues")
    
    if contact_data['structural_residues']:
        print(f"\nStructural residues found: {contact_data['structural_residues']}")
    
    print(f"\nRecommendations:")
    recs = contact_data['recommendations']
    print(f"  Conservative (3.5Å): {len(recs['conservative_3.5A'])} residues")
    print(f"  Liberal (4.0Å): {len(recs['liberal_4.0A'])} residues")
    print(f"  With structural: {len(recs['with_structural'])} residues")
    
    if verbose:
        print(f"\nDetailed residue information:")
        for res in sorted(contact_data['residue_details'].keys()):
            details = contact_data['residue_details'][res]
            distances_str = ", ".join([f"{k}:{v}" for k, v in details['distances'].items() if v is not None])
            print(f"  {info['motif_chain']}:{res}:{details['residue_name']} - {distances_str}")

def main():
    args = parse_args()
    
    # Parse distance cutoffs
    try:
        distance_cutoffs = [float(x.strip()) for x in args.distance_cutoffs.split(',')]
    except ValueError:
        print("Error: Invalid distance cutoffs. Use format: 3.5,4.0,4.5")
        return 1
    
    # Set output filename
    if not args.output:
        base_name = args.pdb.replace('.pdb', '').replace('.PDB', '')
        args.output = f"{base_name}_contacts.json"
    
    try:
        print("Motif Contact Analysis for HPC Systems")
        print("="*50)
        
        # Run analysis
        results, motif_chain, target_chain, all_distances = analyze_contacts(
            args.pdb, distance_cutoffs, args.motif_chain, args.target_chain, 
            args.motif_residues, args.target_residues, args.verbose
        )
        
        # Generate JSON output
        contact_data = generate_contact_json(results, motif_chain, target_chain, args.pdb, all_distances)
        
        # Save JSON
        with open(args.output, 'w') as f:
            json.dump(contact_data, f, indent=2)
        
        # Print summary
        print_summary(contact_data, args.verbose)
        
        print(f"\n✅ Contact analysis complete!")
        print(f"📄 Results saved to: {args.output}")
        print(f"\n🔧 Usage examples:")
        print(f"# Analyze specific motif regions:")
        print(f"python motif_contact_analysis_hpc.py --pdb motif.pdb --motif_residues '66-93,109-127' --output contacts.json")
        print(f"# Analyze single motif region:")
        print(f"python motif_contact_analysis_hpc.py --pdb motif.pdb --motif_chain A --motif_residues '66-93' --target_chain B")
        print(f"# Specify target region too:")
        print(f"python motif_contact_analysis_hpc.py --pdb motif.pdb --motif_residues '66-93' --target_residues '1-100'")
        print(f"\n🔧 Next step - Generate auto-config script:")
        print(f"python auto_config_fixed_script.py \\")
        print(f"  --trb_config your_config.json \\")
        print(f"  --contact_file {args.output} \\")
        print(f"  --distance_cutoff 4.0 \\")
        print(f"  --output selective_fixed_labels_auto.py")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
