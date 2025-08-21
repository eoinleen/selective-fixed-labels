#!/usr/bin/env python
"""
TRB Configuration Parser
Extracts motif mapping information from RFdiffusion TRB files

This script reads a TRB file and extracts:
1. Original motif ranges from con_ref_pdb_idx
2. Scaffold position mappings from con_hal_pdb_idx  
3. Contig mapping string
4. Automatic offset calculations for motif regions

Usage:
python parse_trb_config.py --trb_file example.trb --output_json config.json
python parse_trb_config.py --trb_file example.trb --verbose
"""

import pickle
import argparse
import json
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Extract motif configuration from TRB file')
    parser.add_argument('--trb_file', type=str, required=True, help='Path to TRB file')
    parser.add_argument('--output_json', type=str, help='Output JSON file for configuration')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')
    return parser.parse_args()

def load_trb_file(trb_path):
    """Load TRB file and return the data dictionary"""
    try:
        with open(trb_path, 'rb') as f:
            data = pickle.load(f)
        return data
    except Exception as e:
        raise ValueError(f"Error loading TRB file: {e}")

def parse_motif_regions(con_ref_pdb_idx, con_hal_pdb_idx):
    """
    Parse motif regions from con_ref_pdb_idx and con_hal_pdb_idx
    
    Returns:
    - regions: List of (original_range, scaffold_range, offset) tuples
    - mapping_dict: Dictionary mapping original positions to scaffold positions
    """
    
    # Extract original positions (reference)
    original_positions = [pos for chain, pos in con_ref_pdb_idx]
    scaffold_positions = [pos for chain, pos in con_hal_pdb_idx]
    
    if len(original_positions) != len(scaffold_positions):
        raise ValueError("Mismatch between original and scaffold position counts")
    
    # Group consecutive positions into regions
    regions = []
    mapping_dict = {}
    
    # Create mapping dictionary
    for orig, scaffold in zip(original_positions, scaffold_positions):
        mapping_dict[orig] = scaffold
    
    # Find consecutive regions in original positions
    current_start = original_positions[0]
    current_end = original_positions[0]
    
    for i in range(1, len(original_positions)):
        if original_positions[i] == original_positions[i-1] + 1:
            # Consecutive position
            current_end = original_positions[i]
        else:
            # Gap found - end current region
            scaffold_start = mapping_dict[current_start]
            scaffold_end = mapping_dict[current_end]
            offset = scaffold_start - current_start
            
            regions.append({
                'original_range': (current_start, current_end),
                'scaffold_range': (scaffold_start, scaffold_end),
                'offset': offset,
                'length': current_end - current_start + 1
            })
            
            # Start new region
            current_start = original_positions[i]
            current_end = original_positions[i]
    
    # Add final region
    scaffold_start = mapping_dict[current_start]
    scaffold_end = mapping_dict[current_end]
    offset = scaffold_start - current_start
    
    regions.append({
        'original_range': (current_start, current_end),
        'scaffold_range': (scaffold_start, scaffold_end),
        'offset': offset,
        'length': current_end - current_start + 1
    })
    
    return regions, mapping_dict

def extract_contig_string(config):
    """Extract contig mapping string from config"""
    try:
        contigs = config['contigmap']['contigs']
        if isinstance(contigs, list):
            return contigs[0]
        return contigs
    except KeyError:
        return "Not found in config"

def create_configuration_dict(trb_data):
    """Create configuration dictionary from TRB data"""
    
    config = {
        'trb_info': {},
        'motif_regions': [],
        'mapping': {},
        'contig_string': '',
        'input_pdb': ''
    }
    
    # Extract basic info
    if 'config' in trb_data:
        config['contig_string'] = extract_contig_string(trb_data['config'])
        if 'inference' in trb_data['config']:
            config['input_pdb'] = trb_data['config']['inference'].get('input_pdb', '')
    
    # Parse motif regions
    if 'con_ref_pdb_idx' in trb_data and 'con_hal_pdb_idx' in trb_data:
        regions, mapping_dict = parse_motif_regions(
            trb_data['con_ref_pdb_idx'], 
            trb_data['con_hal_pdb_idx']
        )
        config['motif_regions'] = regions
        config['mapping'] = mapping_dict
    
    return config

def print_configuration(config, verbose=False):
    """Print configuration in human-readable format"""
    
    print("="*80)
    print("TRB CONFIGURATION SUMMARY")
    print("="*80)
    
    print(f"Input PDB: {config['input_pdb']}")
    print(f"Contig String: {config['contig_string']}")
    print(f"Number of Motif Regions: {len(config['motif_regions'])}")
    
    print("\nMotif Region Details:")
    print("-" * 60)
    
    for i, region in enumerate(config['motif_regions'], 1):
        orig_start, orig_end = region['original_range']
        scaffold_start, scaffold_end = region['scaffold_range']
        
        print(f"Region {i}:")
        print(f"  Original Range: A{orig_start}-{orig_end} ({region['length']} residues)")
        print(f"  Scaffold Range: A{scaffold_start}-{scaffold_end}")
        print(f"  Offset: {region['offset']:+d} (scaffold_pos = original_pos {region['offset']:+d})")
        print()
    
    if verbose:
        print("Complete Position Mapping:")
        print("-" * 40)
        for orig, scaffold in sorted(config['mapping'].items()):
            print(f"  A{orig} → A{scaffold}")
    
    print("="*80)

def save_configuration(config, output_path):
    """Save configuration to JSON file"""
    
    # Convert numpy integers to regular integers for JSON serialization
    json_config = {}
    for key, value in config.items():
        if key == 'mapping':
            # Convert mapping keys and values to regular integers
            json_config[key] = {str(k): int(v) for k, v in value.items()}
        elif key == 'motif_regions':
            # Convert motif regions
            json_regions = []
            for region in value:
                json_region = {
                    'original_range': [int(region['original_range'][0]), int(region['original_range'][1])],
                    'scaffold_range': [int(region['scaffold_range'][0]), int(region['scaffold_range'][1])],
                    'offset': int(region['offset']),
                    'length': int(region['length'])
                }
                json_regions.append(json_region)
            json_config[key] = json_regions
        else:
            json_config[key] = value
    
    with open(output_path, 'w') as f:
        json.dump(json_config, f, indent=2)
    
    print(f"Configuration saved to: {output_path}")

def main():
    args = parse_args()
    
    print(f"Loading TRB file: {args.trb_file}")
    
    try:
        # Load TRB file
        trb_data = load_trb_file(args.trb_file)
        
        # Extract configuration
        config = create_configuration_dict(trb_data)
        
        # Print configuration
        print_configuration(config, args.verbose)
        
        # Save to JSON if requested
        if args.output_json:
            save_configuration(config, args.output_json)
        
        print("\nNext steps:")
        print("1. Run contact analysis on your motif structure")
        print("2. Use auto_config_fixed_script.py with this TRB config and contact analysis")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
