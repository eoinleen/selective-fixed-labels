#!/usr/bin/env python
"""
Auto-Configuration Script for Selective FIXED Labels
Automatically generates selective_fixed_labels script from TRB config and contact analysis

This script takes:
1. TRB configuration (from parse_trb_config.py)
2. Contact analysis data (JSON format or manual input)
3. User-defined parameters (distance cutoff, extra residues)

And generates a customized selective_fixed_labels script.

Usage:
python auto_config_fixed_script.py --trb_config config.json --contact_file contacts.json --output selective_fixed_labels_auto.py
python auto_config_fixed_script.py --trb_config config.json --contact_residues "70,71,72,73" --distance_cutoff 4.0
"""

import json
import argparse
import os
import sys
from datetime import datetime

# ============================================================================
# USER CONFIGURATION SECTION - MODIFY THESE SETTINGS
# ============================================================================

# A) DISTANCE CUTOFF OPTIONS (choose one: 3.5, 4.0, or 4.5)
DEFAULT_DISTANCE_CUTOFF = 4.0

# B) EXTRA FIXED RESIDUES (regardless of distance cutoff)
# Add residue numbers that should always be FIXED, separated by commas
# Example: [74, 78, 120] to always fix A74, A78, A120
EXTRA_FIXED_RESIDUES = []

# ============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description='Auto-generate selective FIXED labels script')
    parser.add_argument('--trb_config', type=str, required=True, help='JSON config from parse_trb_config.py')
    parser.add_argument('--contact_file', type=str, help='Contact analysis file (JSON)')
    parser.add_argument('--contact_residues', type=str, help='Comma-separated contact residues (e.g., "70,71,72")')
    parser.add_argument('--distance_cutoff', type=float, choices=[3.5, 4.0, 4.5], 
                       default=DEFAULT_DISTANCE_CUTOFF, help='Distance cutoff for contacts')
    parser.add_argument('--extra_residues', type=str, help='Extra residues to fix (comma-separated)')
    parser.add_argument('--output', type=str, default='selective_fixed_labels_auto.py', 
                       help='Output script filename')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    return parser.parse_args()

def load_trb_config(config_path):
    """Load TRB configuration from JSON file"""
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        return config
    except Exception as e:
        raise ValueError(f"Error loading TRB config: {e}")

def parse_contact_file(contact_file, distance_cutoff):
    """Parse contact analysis file and extract residues within cutoff"""
    
    contact_residues = []
    
    try:
        # Check if it's JSON format (from motif_contact_analysis_hpc.py)
        if contact_file.endswith('.json'):
            with open(contact_file, 'r') as f:
                contact_data = json.load(f)
            
            # Extract residues from JSON format
            cutoff_key = f"{distance_cutoff}A"
            if cutoff_key in contact_data.get('contact_lists', {}):
                contact_residues = contact_data['contact_lists'][cutoff_key]
            else:
                # Find closest available cutoff
                available_cutoffs = []
                for key in contact_data.get('contact_lists', {}).keys():
                    try:
                        cutoff_val = float(key.replace('A', ''))
                        if cutoff_val <= distance_cutoff:
                            available_cutoffs.append((cutoff_val, key))
                    except ValueError:
                        continue
                
                if available_cutoffs:
                    # Use the highest cutoff that's <= requested cutoff
                    best_cutoff = max(available_cutoffs)[1]
                    contact_residues = contact_data['contact_lists'][best_cutoff]
                    print(f"Using {best_cutoff} cutoff (closest to {distance_cutoff}Å)")
        
        else:
            # Parse as TSV/CSV format (original functionality)
            with open(contact_file, 'r') as f:
                lines = f.readlines()
            
            # Skip header if present
            data_lines = [line for line in lines if not line.strip().startswith('#') and line.strip()]
            
            for line in data_lines:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                    
                try:
                    # Assuming format: residue, 3.0A, 3.5A, 4.0A, 4.5A, 5.0A, res_name
                    residue_str = parts[0]  # e.g., "A:70"
                    
                    # Extract residue number
                    if ':' in residue_str:
                        residue_num = int(residue_str.split(':')[1])
                    else:
                        residue_num = int(residue_str.replace('A', ''))
                    
                    # Check distance columns based on cutoff
                    distance_cols = {3.5: 2, 4.0: 3, 4.5: 4, 5.0: 5}
                    col_idx = distance_cols.get(distance_cutoff, 3)
                    
                    if col_idx < len(parts):
                        distance_val = parts[col_idx]
                        if distance_val != '-' and float(distance_val) <= distance_cutoff:
                            contact_residues.append(residue_num)
                            
                except (ValueError, IndexError):
                    continue
                
    except Exception as e:
        raise ValueError(f"Error parsing contact file: {e}")
    
    return sorted(list(set(contact_residues)))

def parse_contact_residues_string(contact_string):
    """Parse comma-separated contact residues string"""
    try:
        residues = [int(r.strip()) for r in contact_string.split(',') if r.strip()]
        return sorted(list(set(residues)))
    except ValueError as e:
        raise ValueError(f"Error parsing contact residues: {e}")

def validate_residues_against_motifs(contact_residues, motif_regions, extra_residues):
    """Validate that contact residues fall within motif regions"""
    
    # Get all valid motif residue ranges
    valid_ranges = []
    for region in motif_regions:
        start, end = region['original_range']
        valid_ranges.append((start, end))
    
    # Check contact residues
    invalid_contact = []
    for residue in contact_residues:
        valid = any(start <= residue <= end for start, end in valid_ranges)
        if not valid:
            invalid_contact.append(residue)
    
    # Check extra residues
    invalid_extra = []
    for residue in extra_residues:
        valid = any(start <= residue <= end for start, end in valid_ranges)
        if not valid:
            invalid_extra.append(residue)
    
    # Report conflicts
    if invalid_contact or invalid_extra:
        print("=" * 80)
        print("⚠️  RESIDUE VALIDATION ERRORS")
        print("=" * 80)
        
        if invalid_contact:
            print(f"❌ Contact residues NOT in motif regions: {invalid_contact}")
        
        if invalid_extra:
            print(f"❌ Extra residues NOT in motif regions: {invalid_extra}")
        
        print("\nValid motif ranges:")
        for i, (start, end) in enumerate(valid_ranges, 1):
            print(f"  Region {i}: A{start}-{end}")
        
        print("\nPlease check your contact analysis or TRB file.")
        print("=" * 80)
        return False
    
    return True

def assign_residues_to_regions(residues, motif_regions):
    """Assign residues to their corresponding motif regions"""
    
    region_assignments = {}
    
    for i, region in enumerate(motif_regions):
        start, end = region['original_range']
        region_residues = [r for r in residues if start <= r <= end]
        if region_residues:
            region_assignments[i] = region_residues
    
    return region_assignments

def generate_script_content(trb_config, contact_residues, extra_residues, distance_cutoff):
    """Generate the selective fixed labels script content"""
    
    # Combine contact and extra residues
    all_residues = sorted(list(set(contact_residues + extra_residues)))
    
    # Assign residues to regions
    region_assignments = assign_residues_to_regions(all_residues, trb_config['motif_regions'])
    
    # Generate timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Build contig string
    contig_string = trb_config['contig_string']
    
    # Create region info for documentation
    region_info = []
    for i, region in enumerate(trb_config['motif_regions']):
        start, end = region['original_range']
        length = region['length']
        offset = region['offset']
        region_info.append(f"Region {i+1}: A{start}-{end} (offset: {offset:+d})")
    
    # Generate contact lists for each region
    region_contact_lists = []
    
    for i, region in enumerate(trb_config['motif_regions']):
        region_residues = region_assignments.get(i, [])
        region_contact_lists.append(region_residues)
    
    # Determine if single or multiple motifs
    num_regions = len(trb_config['motif_regions'])
    is_single_motif = num_regions == 1
    
    # Start building the script content
    script_lines = []
    
    script_lines.append("#!/usr/bin/env python")
    script_lines.append('"""')
    script_lines.append("Auto-Generated Selective FIXED Labels Script")
    script_lines.append(f"Generated on: {timestamp}")
    script_lines.append(f"Distance cutoff: {distance_cutoff}Å")
    script_lines.append(f"Extra residues: {extra_residues}")
    script_lines.append("")
    script_lines.append("=== CURRENT CONFIGURATION ===")
    script_lines.append(f"TRB file: {trb_config.get('input_pdb', 'Unknown')}")
    script_lines.append(f"Contig mapping: {contig_string}")
    script_lines.append(f"Motif regions: {len(trb_config['motif_regions'])} region{'s' if not is_single_motif else ''}")
    for region_line in region_info:
        script_lines.append(region_line)
    script_lines.append(f"Contact analysis with {distance_cutoff}Å distance cutoff")
    script_lines.append("")
    script_lines.append("Usage:")
    script_lines.append("python selective_fixed_labels_auto.py --pdbdir . --verbose --dry_run")
    script_lines.append("python selective_fixed_labels_auto.py --pdbdir . --verbose")
    script_lines.append('"""')
    script_lines.append("")
    script_lines.append("import os")
    script_lines.append("import argparse")
    script_lines.append("from collections import defaultdict")
    script_lines.append("")
    
    # Parse args function
    script_lines.append("def parse_args():")
    script_lines.append("    parser = argparse.ArgumentParser(description='Apply selective FIXED labels based on auto-generated configuration')")
    script_lines.append("    parser.add_argument('--pdbdir', type=str, required=True, help='Directory containing PDB files with FIXED labels')")
    script_lines.append("    parser.add_argument('--chain', type=str, default='A', help='Chain to modify')")
    script_lines.append("    parser.add_argument('--verbose', action='store_true', help='Verbose output')")
    script_lines.append("    parser.add_argument('--dry_run', action='store_true', help='Show what would be changed without modifying files')")
    script_lines.append("    return parser.parse_args()")
    script_lines.append("")
    
    # Contact residues function
    script_lines.append("def get_critical_contact_residues():")
    script_lines.append('    """')
    script_lines.append("    Returns original motif residue numbers that should remain FIXED")
    script_lines.append(f"    Based on ≤{distance_cutoff}Å contact analysis + extra specified residues")
    script_lines.append("    ")
    script_lines.append(f"    Contact residues: {contact_residues}")
    script_lines.append(f"    Extra residues: {extra_residues}")
    script_lines.append(f"    Combined: {all_residues}")
    script_lines.append('    """')
    
    # Add region contact lists
    for i, region_residues in enumerate(region_contact_lists):
        start, end = trb_config['motif_regions'][i]['original_range']
        script_lines.append(f"    region{i+1}_contacts = {region_residues}  # Original A{start}-{end} region")
    
    script_lines.append("    ")
    return_stmt = ', '.join([f'region{i+1}_contacts' for i in range(len(region_contact_lists))])
    script_lines.append(f"    return {return_stmt}")
    script_lines.append("")
    
    # Find fixed regions function
    script_lines.append("def find_fixed_regions(pdb_file):")
    script_lines.append('    """')
    script_lines.append(f"    Extract FIXED positions and identify {len(trb_config['motif_regions'])} separate region{'s' if not is_single_motif else ''}")
    if is_single_motif:
        script_lines.append("    Returns list of FIXED positions")
    else:
        script_lines.append("    Returns tuple of region position lists")
    script_lines.append('    """')
    script_lines.append("    fixed_positions = []")
    script_lines.append("    with open(pdb_file, 'r') as f:")
    script_lines.append("        for line in f:")
    script_lines.append("            if line.startswith('REMARK PDBinfo-LABEL:') and 'FIXED' in line:")
    script_lines.append("                parts = line.split()")
    script_lines.append("                try:")
    script_lines.append("                    pos = int(parts[2])")
    script_lines.append("                    fixed_positions.append(pos)")
    script_lines.append("                except (ValueError, IndexError):")
    script_lines.append("                    continue")
    script_lines.append("    ")
    script_lines.append("    fixed_positions = sorted(fixed_positions)")
    script_lines.append("    ")
    script_lines.append("    if len(fixed_positions) == 0:")
    if is_single_motif:
        script_lines.append("        return []")
    else:
        script_lines.append(f"        return tuple([[] for _ in range({num_regions})])")
    
    if is_single_motif:
        script_lines.append("    # Single motif - return all FIXED positions as one region")
        script_lines.append("    return fixed_positions")
    else:
        script_lines.append(f"    # Find gaps between regions (expecting {len(trb_config['motif_regions'])} regions)")
        script_lines.append("    gaps = []")
        script_lines.append("    for i in range(1, len(fixed_positions)):")
        script_lines.append("        gap = fixed_positions[i] - fixed_positions[i-1]")
        script_lines.append("        if gap > 10:  # Significant gap indicates separate regions")
        script_lines.append("            gaps.append(i)")
        script_lines.append("    ")
        script_lines.append("    # Split into regions")
        script_lines.append("    regions = []")
        script_lines.append("    start_idx = 0")
        script_lines.append("    ")
        script_lines.append("    for gap_idx in gaps:")
        script_lines.append("        regions.append(fixed_positions[start_idx:gap_idx])")
        script_lines.append("        start_idx = gap_idx")
        script_lines.append("    ")
        script_lines.append("    # Add final region")
        script_lines.append("    regions.append(fixed_positions[start_idx:])")
        script_lines.append("    ")
        script_lines.append(f"    # Ensure we have the expected number of regions")
        script_lines.append(f"    while len(regions) < {len(trb_config['motif_regions'])}:")
        script_lines.append("        regions.append([])")
        script_lines.append("    ")
        script_lines.append(f"    return tuple(regions[:{len(trb_config['motif_regions'])}])")
    
    script_lines.append("")
    
    # Mapping function
    if is_single_motif:
        region = trb_config['motif_regions'][0]
        start, end = region['original_range']
        script_lines.append("def map_contacts_to_scaffold_positions(region1_pos, region1_contacts):")
        script_lines.append('    """')
        script_lines.append("    Map original motif contact residues to scaffold positions using TRB-derived offsets")
        script_lines.append("    Single motif version")
        script_lines.append('    """')
        script_lines.append("    critical_scaffold_positions = []")
        script_lines.append("    ")
        script_lines.append(f"    # Map region 1 contacts (A{start}-{end})")
        script_lines.append("    if region1_pos:")
        script_lines.append("        region1_start = region1_pos[0]")
        script_lines.append("        for original_pos in region1_contacts:")
        script_lines.append(f"            if {start} <= original_pos <= {end}:")
        script_lines.append(f"                scaffold_pos = original_pos - {start} + region1_start")
        script_lines.append("                if scaffold_pos in region1_pos:")
        script_lines.append("                    critical_scaffold_positions.append(scaffold_pos)")
        script_lines.append("    ")
        script_lines.append("    return sorted(critical_scaffold_positions)")
    else:
        region_params = ', '.join([f'region{i+1}_pos' for i in range(len(trb_config['motif_regions']))])
        contact_params = ', '.join([f'region{i+1}_contacts' for i in range(len(trb_config['motif_regions']))])
        script_lines.append(f"def map_contacts_to_scaffold_positions({region_params}, {contact_params}):")
        script_lines.append('    """')
        script_lines.append("    Map original motif contact residues to scaffold positions using TRB-derived offsets")
        script_lines.append('    """')
        script_lines.append("    critical_scaffold_positions = []")
        script_lines.append("    ")
        
        for i, region in enumerate(trb_config['motif_regions']):
            start, end = region['original_range']
            script_lines.append(f"    # Map region {i+1} contacts (A{start}-{end})")
            script_lines.append(f"    if region{i+1}_pos:")
            script_lines.append(f"        region{i+1}_start = region{i+1}_pos[0]")
            script_lines.append(f"        for original_pos in region{i+1}_contacts:")
            script_lines.append(f"            if {start} <= original_pos <= {end}:")
            script_lines.append(f"                scaffold_pos = original_pos - {start} + region{i+1}_start")
            script_lines.append(f"                if scaffold_pos in region{i+1}_pos:")
            script_lines.append("                    critical_scaffold_positions.append(scaffold_pos)")
            script_lines.append("    ")
        
        script_lines.append("    return sorted(critical_scaffold_positions)")
    
    script_lines.append("")
    
    # Remove function
    script_lines.append("def remove_non_critical_fixed_labels(pdb_file, critical_positions, dry_run=False):")
    script_lines.append('    """Remove FIXED labels from positions not in critical_positions list"""')
    script_lines.append("    if dry_run:")
    script_lines.append("        return True, 0, 0")
    script_lines.append("")
    script_lines.append("    with open(pdb_file, 'r') as f:")
    script_lines.append("        lines = f.readlines()")
    script_lines.append("")
    script_lines.append("    new_lines = []")
    script_lines.append("    removed_count = 0")
    script_lines.append("    kept_count = 0")
    script_lines.append("")
    script_lines.append("    for line in lines:")
    script_lines.append("        if line.startswith('REMARK PDBinfo-LABEL:') and 'FIXED' in line:")
    script_lines.append("            parts = line.split()")
    script_lines.append("            try:")
    script_lines.append("                pos = int(parts[2])")
    script_lines.append("                if pos in critical_positions:")
    script_lines.append("                    new_lines.append(line)")
    script_lines.append("                    kept_count += 1")
    script_lines.append("                else:")
    script_lines.append("                    removed_count += 1")
    script_lines.append("            except (ValueError, IndexError):")
    script_lines.append("                new_lines.append(line)")
    script_lines.append("        else:")
    script_lines.append("            new_lines.append(line)")
    script_lines.append("")
    script_lines.append("    if removed_count > 0 or kept_count > 0:")
    script_lines.append("        with open(pdb_file, 'w') as f:")
    script_lines.append("            f.writelines(new_lines)")
    script_lines.append("        return True, removed_count, kept_count")
    script_lines.append("")
    script_lines.append("    return False, 0, 0")
    script_lines.append("")
    
    # Process function
    contact_params = ', '.join([f'region{i+1}_contacts' for i in range(len(trb_config['motif_regions']))])
    script_lines.append(f"def process_scaffold(pdb_file, {contact_params}, verbose, dry_run):")
    script_lines.append('    """Process a single scaffold file"""')
    script_lines.append("    basename = os.path.basename(pdb_file)")
    script_lines.append("")
    script_lines.append("    if verbose:")
    script_lines.append("        print(f\"\\nProcessing: {basename}\")")
    script_lines.append("")
    script_lines.append("    try:")
    script_lines.append("        # Find FIXED regions")
    if is_single_motif:
        script_lines.append("        region1_pos = find_fixed_regions(pdb_file)")
    else:
        region_vars = ', '.join([f'region{i+1}_pos' for i in range(len(trb_config['motif_regions']))])
        script_lines.append(f"        {region_vars} = find_fixed_regions(pdb_file)")
    script_lines.append("        ")
    script_lines.append("        if verbose:")
    
    if is_single_motif:
        script_lines.append("            if region1_pos:")
        script_lines.append("                print(f\"  Region 1: {len(region1_pos)} positions ({region1_pos[0]}-{region1_pos[-1]})\")")
    else:
        for i in range(len(trb_config['motif_regions'])):
            script_lines.append(f"            if region{i+1}_pos:")
            script_lines.append(f"                print(f\"  Region {i+1}: {{len(region{i+1}_pos)}} positions ({{region{i+1}_pos[0]}}-{{region{i+1}_pos[-1]}})\")")
    
    script_lines.append("")
    script_lines.append("        # Map contact residues to scaffold positions")
    if is_single_motif:
        script_lines.append("        critical_scaffold_positions = map_contacts_to_scaffold_positions(")
        script_lines.append("            region1_pos, region1_contacts")
        script_lines.append("        )")
    else:
        region_params = ', '.join([f'region{i+1}_pos' for i in range(len(trb_config['motif_regions']))])
        contact_params = ', '.join([f'region{i+1}_contacts' for i in range(len(trb_config['motif_regions']))])
        script_lines.append("        critical_scaffold_positions = map_contacts_to_scaffold_positions(")
        script_lines.append(f"            {region_params}, ")
        script_lines.append(f"            {contact_params}")
        script_lines.append("        )")
    
    script_lines.append("")
    script_lines.append("        if verbose:")
    script_lines.append("            print(f\"  Critical contacts mapped to scaffold: {critical_scaffold_positions}\")")
    script_lines.append("")
    script_lines.append("        # Determine changes needed")
    if is_single_motif:
        script_lines.append("        all_fixed = list(region1_pos)")
    else:
        script_lines.append("        all_fixed = []")
        for i in range(len(trb_config['motif_regions'])):
            script_lines.append(f"        all_fixed.extend(list(region{i+1}_pos))")
    
    script_lines.append("        positions_to_remove = [pos for pos in all_fixed if pos not in critical_scaffold_positions]")
    script_lines.append("        positions_to_keep = [pos for pos in all_fixed if pos in critical_scaffold_positions]")
    script_lines.append("")
    script_lines.append("        if verbose:")
    script_lines.append("            print(f\"  Will remove FIXED from: {positions_to_remove}\")")
    script_lines.append("            print(f\"  Will keep FIXED on: {positions_to_keep}\")")
    script_lines.append("")
    script_lines.append("        if not positions_to_remove:")
    script_lines.append("            print(f\"  ✅ No changes needed - all FIXED positions are critical\")")
    script_lines.append("            return {'status': 'no_change', 'file': basename}")
    script_lines.append("")
    script_lines.append("        # Apply changes")
    script_lines.append("        if dry_run:")
    script_lines.append("            print(f\"  [DRY RUN] Would remove {len(positions_to_remove)} FIXED labels\")")
    script_lines.append("            print(f\"  [DRY RUN] Would keep {len(positions_to_keep)} FIXED labels\")")
    script_lines.append("            return {'status': 'dry_run', 'file': basename,")
    script_lines.append("                    'removed': len(positions_to_remove), 'kept': len(positions_to_keep)}")
    script_lines.append("        else:")
    script_lines.append("            success, removed_count, kept_count = remove_non_critical_fixed_labels(")
    script_lines.append("                pdb_file, critical_scaffold_positions, dry_run)")
    script_lines.append("")
    script_lines.append("            if success:")
    script_lines.append("                print(f\"  ✅ Removed {removed_count} FIXED labels, kept {kept_count}\")")
    script_lines.append("                return {'status': 'modified', 'file': basename,")
    script_lines.append("                        'removed': removed_count, 'kept': kept_count}")
    script_lines.append("            else:")
    script_lines.append("                print(f\"  ❌ Failed to modify FIXED labels\")")
    script_lines.append("                return {'status': 'failed', 'file': basename}")
    script_lines.append("")
    script_lines.append("    except Exception as e:")
    script_lines.append("        print(f\"  ❌ Error processing file: {e}\")")
    script_lines.append("        return {'status': 'error', 'file': basename, 'error': str(e)}")
    script_lines.append("")
    
    # Main function
    script_lines.append("def main():")
    script_lines.append("    args = parse_args()")
    script_lines.append("")
    script_lines.append("    # Get critical contact residues")
    return_vars = ', '.join([f'region{i+1}_contacts' for i in range(len(region_contact_lists))])
    script_lines.append(f"    {return_vars} = get_critical_contact_residues()")
    script_lines.append("")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("    print(\"AUTO-GENERATED SELECTIVE FIXED LABELS SCRIPT\")")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("    print(f\"Target directory: {args.pdbdir}\")")
    script_lines.append("    print(f\"Chain: {args.chain}\")")
    script_lines.append("    print(f\"Mode: {'DRY RUN' if args.dry_run else 'LIVE MODIFICATION'}\")")
    script_lines.append(f"    print(f\"Distance cutoff: {distance_cutoff}Å\")")
    script_lines.append(f"    print(f\"Contig: {contig_string}\")")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("")
    script_lines.append(f"    # Contact analysis summary")
    script_lines.append(f"    print(\"Contact analysis basis (≤{distance_cutoff}Å contacts + extra residues):\")")
    
    for i, region in enumerate(trb_config['motif_regions']):
        start, end = region['original_range']
        script_lines.append(f"    print(\"Region {i+1} (A{start}-{end}):\")")
        script_lines.append(f"    for orig_pos in region{i+1}_contacts:")
        script_lines.append("        print(f\"        A{orig_pos}\")")
    
    script_lines.append("    ")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("")
    script_lines.append("    # Find PDB files")
    script_lines.append("    pdb_files = [f for f in os.listdir(args.pdbdir) if f.endswith('.pdb')]")
    script_lines.append("    if not pdb_files:")
    script_lines.append("        print(\"No PDB files found in directory\")")
    script_lines.append("        return")
    script_lines.append("")
    script_lines.append("    print(f\"Found {len(pdb_files)} PDB files\")")
    script_lines.append("")
    script_lines.append("    # Process each file")
    script_lines.append("    results = []")
    script_lines.append("    for pdb_file in sorted(pdb_files):")
    script_lines.append("        pdb_path = os.path.join(args.pdbdir, pdb_file)")
    contact_args = ', '.join([f'region{i+1}_contacts' for i in range(len(region_contact_lists))])
    script_lines.append(f"        result = process_scaffold(pdb_path, {contact_args}, args.verbose, args.dry_run)")
    script_lines.append("        results.append(result)")
    script_lines.append("")
    script_lines.append("    # Summary")
    script_lines.append("    print(\"\\n\" + \"=\"*80)")
    script_lines.append("    print(\"SUMMARY\")")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("")
    script_lines.append("    status_counts = defaultdict(int)")
    script_lines.append("    total_removed = 0")
    script_lines.append("    total_kept = 0")
    script_lines.append("")
    script_lines.append("    for result in results:")
    script_lines.append("        status_counts[result['status']] += 1")
    script_lines.append("        if 'removed' in result:")
    script_lines.append("            total_removed += result['removed']")
    script_lines.append("        if 'kept' in result:")
    script_lines.append("            total_kept += result['kept']")
    script_lines.append("")
    script_lines.append("    for status, count in status_counts.items():")
    script_lines.append("        print(f\"{status}: {count} files\")")
    script_lines.append("")
    script_lines.append("    if total_removed > 0 or total_kept > 0:")
    script_lines.append("        print(f\"\\nTotal FIXED labels removed: {total_removed}\")")
    script_lines.append("        print(f\"Total FIXED labels kept: {total_kept}\")")
    script_lines.append("")
    script_lines.append("    print(\"\\nNext steps:\")")
    script_lines.append("    print(\"1. Run ProteinMPNN with these selective FIXED labels\")")
    script_lines.append("    print(\"2. Apply targeted mutations based on contact analysis\")")
    script_lines.append("    print(\"=\"*80)")
    script_lines.append("")
    script_lines.append("if __name__ == \"__main__\":")
    script_lines.append("    main()")
    
    # Join all lines with newlines
    script_content = '\n'.join(script_lines)
    
    return script_content

    return script_content

def main():
    args = parse_args()
    
    print("="*80)
    print("AUTO-CONFIGURATION SCRIPT FOR SELECTIVE FIXED LABELS")
    print("="*80)
    print(f"Distance cutoff: {args.distance_cutoff}Å")
    print(f"Extra residues: {args.extra_residues if args.extra_residues else 'None'}")
    print("="*80)
    
    try:
        # Load TRB configuration
        print(f"Loading TRB config: {args.trb_config}")
        trb_config = load_trb_config(args.trb_config)
        
        # Parse contact residues
        contact_residues = []
        
        if args.contact_file:
            print(f"Parsing contact file: {args.contact_file}")
            contact_residues = parse_contact_file(args.contact_file, args.distance_cutoff)
        elif args.contact_residues:
            print(f"Using provided contact residues: {args.contact_residues}")
            contact_residues = parse_contact_residues_string(args.contact_residues)
        else:
            print("❌ Error: Must provide either --contact_file or --contact_residues")
            return 1
        
        # Parse extra residues
        extra_residues = EXTRA_FIXED_RESIDUES.copy()
        if args.extra_residues:
            extra_residues.extend(parse_contact_residues_string(args.extra_residues))
        
        print(f"Contact residues (≤{args.distance_cutoff}Å): {contact_residues}")
        print(f"Extra residues: {extra_residues}")
        
        # Validate residues against motif regions
        if not validate_residues_against_motifs(contact_residues, trb_config['motif_regions'], extra_residues):
            return 1
        
        print("✅ All residues validated against motif regions")
        
        # Generate script
        print(f"Generating script: {args.output}")
        script_content = generate_script_content(trb_config, contact_residues, extra_residues, args.distance_cutoff)
        
        # Write script
        with open(args.output, 'w') as f:
            f.write(script_content)
        
        # Make executable
        import stat
        st = os.stat(args.output)
        os.chmod(args.output, st.st_mode | stat.S_IEXEC)
        
        print(f"✅ Script generated: {args.output}")
        print("\nNext steps:")
        print(f"1. Test: python {args.output} --pdbdir /path/to/pdbs --verbose --dry_run")
        print(f"2. Apply: python {args.output} --pdbdir /path/to/pdbs --verbose")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
