#!/usr/bin/env python3
"""
Validation Script for Automated Selective FIXED Labels Pipeline

Validates that the selective FIXED label application worked correctly by:
1. Comparing original motif contact analysis results
2. Verifying TRB mapping (original → scaffold positions) 
3. Checking final FIXED labels in processed PDB files
4. Performing structural validation checks
5. Generating comprehensive accuracy and coverage reports

Usage:
    python validate_selective_fixed_pipeline.py \
        --pdb_dir /path/to/processed/pdbs \
        --trb_dir /path/to/trb/files \
        --contacts_file contacts.json \
        --motif_pdb original_motif.pdb \
        --distance_cutoff 4.0 \
        --output_report validation_report.txt \
        --verbose
"""

import os
import json
import argparse
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Set, Tuple, Optional

# Import required modules (adjust paths as needed)
try:
    from Bio import PDB
    from Bio.PDB import PDBParser, PDBIO, Select
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("Warning: BioPython not found. Basic validation will still work.")
    BIOPYTHON_AVAILABLE = False

class ValidationError(Exception):
    """Custom exception for validation errors"""
    pass

class FixedLabelValidator:
    """Validates the selective FIXED label pipeline results"""
    
    def __init__(self, distance_cutoff: float = 4.0, verbose: bool = False):
        self.distance_cutoff = distance_cutoff
        self.verbose = verbose
        
        # Initialize PDB parser if BioPython is available
        if BIOPYTHON_AVAILABLE:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
        
        # Results storage
        self.validation_results = {
            'summary': {},
            'file_results': {},
            'errors': [],
            'warnings': []
        }
        
    def log(self, message: str, level: str = "INFO"):
        """Log messages with optional verbosity"""
        if self.verbose or level in ["ERROR", "WARNING"]:
            prefix = f"[{level}]" if level != "INFO" else ""
            print(f"{prefix} {message}")
    
    def load_contacts_data(self, contacts_file: str) -> Dict:
        """Load the original motif contact analysis results"""
        try:
            with open(contacts_file, 'r') as f:
                contacts_data = json.load(f)
            self.log(f"Loaded contacts data from {contacts_file}")
            return contacts_data
        except Exception as e:
            raise ValidationError(f"Failed to load contacts file {contacts_file}: {e}")
    
    def parse_trb_file(self, trb_file: str) -> Dict:
        """Parse TRB file to get residue mapping"""
        mapping = {
            'con_ref_pdb_idx': [],
            'con_hal_pdb_idx': [],
            'motif_regions': [],
            'region_offsets': {}
        }
        
        try:
            # TRB files are typically pickled/binary files
            import pickle
            
            with open(trb_file, 'rb') as f:
                trb_data = pickle.load(f)
            
            # Extract mapping information from TRB data
            if isinstance(trb_data, dict):
                if 'con_ref_pdb_idx' in trb_data and 'con_hal_pdb_idx' in trb_data:
                    con_ref = trb_data['con_ref_pdb_idx']
                    con_hal = trb_data['con_hal_pdb_idx']
                    
                    # Extract residue numbers from tuples (chain, residue_number)
                    ref_positions = [int(item[1]) for item in con_ref]  # Extract residue numbers
                    hal_positions = [int(item[1]) for item in con_hal]  # Extract residue numbers
                    
                    mapping['con_ref_pdb_idx'] = ref_positions
                    mapping['con_hal_pdb_idx'] = hal_positions
                    
                    # Calculate motif regions and offsets (same logic as your pipeline)
                    if ref_positions and hal_positions:
                        # Group consecutive positions into regions
                        regions = []
                        current_region_start = ref_positions[0]
                        current_hal_start = hal_positions[0]
                        
                        for i in range(1, len(ref_positions)):
                            # Check for gap (indicates new motif region)
                            if ref_positions[i] != ref_positions[i-1] + 1:
                                # End current region
                                regions.append({
                                    'original_range': [current_region_start, ref_positions[i-1]],
                                    'scaffold_range': [current_hal_start, hal_positions[i-1]],
                                    'offset': current_hal_start - current_region_start
                                })
                                # Start new region
                                current_region_start = ref_positions[i]
                                current_hal_start = hal_positions[i]
                        
                        # Add final region
                        regions.append({
                            'original_range': [current_region_start, ref_positions[-1]],
                            'scaffold_range': [current_hal_start, hal_positions[-1]],
                            'offset': current_hal_start - current_region_start
                        })
                        
                        mapping['motif_regions'] = regions
                        
                        # Calculate offsets for each position
                        for i, ref_pos in enumerate(ref_positions):
                            hal_pos = hal_positions[i]
                            offset = hal_pos - ref_pos
                            mapping['region_offsets'][ref_pos] = offset
                        
                self.log(f"Parsed TRB mapping: {len(mapping.get('con_ref_pdb_idx', []))} residues, {len(mapping['motif_regions'])} regions")
                
                # Log the regions for debugging
                for i, region in enumerate(mapping['motif_regions']):
                    self.log(f"Region {i+1}: {region['original_range'][0]}-{region['original_range'][1]} -> {region['scaffold_range'][0]}-{region['scaffold_range'][1]} (offset: {region['offset']})")
                
                return mapping
            
        except Exception as e:
            self.log(f"Warning: Could not parse TRB file {trb_file}: {e}", "WARNING")
            
        return mapping
    
    def get_expected_fixed_residues(self, contacts_data: Dict, 
                                  distance_cutoff: float) -> Set[int]:
        """Get expected FIXED residues based on contact analysis"""
        expected_fixed = set()
        
        try:
            self.log(f"Looking for expected residues with cutoff {distance_cutoff}")
            
            # First try to get from recommendations (these are the processed lists)
            if 'recommendations' in contacts_data:
                recommendations = contacts_data['recommendations']
                self.log(f"Found recommendations: {list(recommendations.keys())}")
                
                # Try different recommendation keys based on cutoff
                cutoff_str = str(distance_cutoff)
                rec_keys = [
                    f"liberal_{cutoff_str}A",  # e.g., "liberal_4.0A"
                    f"conservative_{cutoff_str}A",  # e.g., "conservative_3.5A"
                    "with_structural"  # This includes structural residues
                ]
                
                found = False
                for key in rec_keys:
                    if key in recommendations:
                        contact_residues = recommendations[key]
                        if isinstance(contact_residues, list):
                            expected_fixed.update(contact_residues)
                            found = True
                            self.log(f"Found {len(contact_residues)} recommended residues for '{key}': {sorted(contact_residues)}")
                            break
                
                if found:
                    return expected_fixed
                else:
                    self.log(f"No matching recommendation found for cutoff {distance_cutoff}", "WARNING")
            else:
                self.log("No 'recommendations' section found", "WARNING")
            
            # Fallback: try contact_lists (raw contact data)
            if 'contact_lists' in contacts_data:
                contact_lists = contacts_data['contact_lists']
                self.log(f"Found contact_lists: {list(contact_lists.keys())}")
                
                # Try different possible keys for the distance cutoff
                cutoff_str = str(distance_cutoff)
                possible_keys = [
                    f"{cutoff_str}A",  # Format: "4.0A"
                    cutoff_str,        # Format: "4.0"
                    f"cutoff_{cutoff_str}",
                    f"distance_{cutoff_str}"
                ]
                
                found = False
                for key in possible_keys:
                    if key in contact_lists:
                        contact_residues = contact_lists[key]
                        if isinstance(contact_residues, list):
                            expected_fixed.update(contact_residues)
                            found = True
                            self.log(f"Found {len(contact_residues)} contact residues for cutoff {distance_cutoff}A")
                            break
                
                if not found:
                    # Print available keys for debugging
                    available_cutoffs = list(contact_lists.keys()) if contact_lists else []
                    self.log(f"Available cutoffs in contact_lists: {available_cutoffs}", "WARNING")
                    
            else:
                # Print available keys for debugging
                self.log(f"Available keys in contacts data: {list(contacts_data.keys())}", "WARNING")
                self.log(f"Warning: Neither 'recommendations' nor 'contact_lists' found in contacts data", "WARNING")
                
        except Exception as e:
            self.log(f"Error processing contacts data: {e}", "ERROR")
            
        self.log(f"Final expected_fixed set: {len(expected_fixed)} residues: {sorted(expected_fixed)}")
        return expected_fixed
    
    def get_fixed_residues_from_pdb(self, pdb_file: str) -> Set[int]:
        """Extract FIXED residues from PDB file"""
        fixed_residues = set()
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    # Look for REMARK lines with FIXED labels
                    if line.startswith('REMARK') and 'FIXED' in line:
                        # Format: REMARK PDBinfo-LABEL:   49 FIXED
                        parts = line.strip().split()
                        for i, part in enumerate(parts):
                            if part == 'FIXED' and i > 0:
                                # The residue number should be just before 'FIXED'
                                try:
                                    res_num = int(parts[i-1])
                                    fixed_residues.add(res_num)
                                    break
                                except (ValueError, IndexError):
                                    continue
                    
                    # Also check ATOM lines with FIXED in them (backup method)
                    elif line.startswith('ATOM') and 'FIXED' in line:
                        # Extract residue number from PDB ATOM line (columns 23-26)
                        try:
                            res_num = int(line[22:26].strip())
                            fixed_residues.add(res_num)
                        except (ValueError, IndexError):
                            continue
                        
        except Exception as e:
            self.log(f"Error reading PDB file {pdb_file}: {e}", "ERROR")
            
        return fixed_residues
    
    def load_original_motif_coordinates(self, motif_pdb_file: str) -> Dict:
        """Load coordinates from original motif PDB for structural validation"""
        motif_coords = {}
        
        if not motif_pdb_file or not os.path.exists(motif_pdb_file):
            return motif_coords
            
        try:
            with open(motif_pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21]
                        res_num = int(line[22:26].strip())
                        atom_name = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        if chain not in motif_coords:
                            motif_coords[chain] = {}
                        if res_num not in motif_coords[chain]:
                            motif_coords[chain][res_num] = {}
                            
                        motif_coords[chain][res_num][atom_name] = (x, y, z)
                        
            self.log(f"Loaded coordinates for {sum(len(chain_data) for chain_data in motif_coords.values())} residues from motif PDB")
            
        except Exception as e:
            self.log(f"Warning: Could not load motif coordinates from {motif_pdb_file}: {e}", "WARNING")
            
        return motif_coords
    
    def get_motif_regions_from_trb(self, trb_mapping: Dict) -> List[Tuple[int, int]]:
        """Extract motif region ranges from TRB mapping"""
        regions = []
        
        if 'motif_regions' in trb_mapping:
            for region in trb_mapping['motif_regions']:
                start, end = region['original_range']
                regions.append((start, end))
        elif 'con_ref_pdb_idx' in trb_mapping:
            # Fallback: derive regions from reference positions
            ref_positions = trb_mapping['con_ref_pdb_idx']
            if ref_positions:
                regions.append((min(ref_positions), max(ref_positions)))
                
        return regions
    
    def get_all_scaffold_residues(self, pdb_file: str) -> Set[int]:
        """Get all residue numbers present in scaffold PDB"""
        all_residues = set()
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        res_num = int(line[22:26].strip())
                        all_residues.add(res_num)
                        
        except Exception as e:
            self.log(f"Error reading scaffold residues from {pdb_file}: {e}", "ERROR")
            
        return all_residues
    
    def validate_structural_consistency(self, pdb_file: str, trb_mapping: Dict, 
                                      actual_fixed: Set[int], expected_fixed: Set[int],
                                      motif_coords: Dict) -> Dict:
        """Perform structural validation checks"""
        structural_results = {
            'motif_region_check': {'status': 'pass', 'details': []},
            'non_motif_freedom_check': {'status': 'pass', 'details': []},
            'coordinate_validation': {'status': 'pass', 'details': []}
        }
        
        try:
            # 1. Verify FIXED residues are in motif regions
            motif_regions = self.get_motif_regions_from_trb(trb_mapping)
            if motif_regions:
                # Map motif regions to scaffold positions
                scaffold_motif_ranges = []
                for orig_start, orig_end in motif_regions:
                    if 'region_offsets' in trb_mapping and orig_start in trb_mapping['region_offsets']:
                        offset = trb_mapping['region_offsets'][orig_start]
                        scaffold_start = orig_start + offset
                        scaffold_end = orig_end + offset
                        scaffold_motif_ranges.append((scaffold_start, scaffold_end))
                
                # Check if all FIXED residues are within motif regions
                fixed_outside_motif = set()
                for res in actual_fixed:
                    in_motif = any(start <= res <= end for start, end in scaffold_motif_ranges)
                    if not in_motif:
                        fixed_outside_motif.add(res)
                
                if fixed_outside_motif:
                    structural_results['motif_region_check']['status'] = 'warning'
                    structural_results['motif_region_check']['details'].append(
                        f"FIXED residues outside motif regions: {sorted(fixed_outside_motif)}"
                    )
                else:
                    structural_results['motif_region_check']['details'].append(
                        f"All {len(actual_fixed)} FIXED residues are within motif regions"
                    )
            
            # 2. Check that non-motif residues are properly freed
            all_scaffold_residues = self.get_all_scaffold_residues(pdb_file)
            if motif_regions and 'region_offsets' in trb_mapping:
                # Determine which scaffold residues correspond to motif
                scaffold_motif_residues = set()
                for orig_start, orig_end in motif_regions:
                    if orig_start in trb_mapping['region_offsets']:
                        offset = trb_mapping['region_offsets'][orig_start]
                        for orig_res in range(orig_start, orig_end + 1):
                            scaffold_res = orig_res + offset
                            scaffold_motif_residues.add(scaffold_res)
                
                # Check non-motif residues are not FIXED
                non_motif_residues = all_scaffold_residues - scaffold_motif_residues
                non_motif_fixed = actual_fixed.intersection(non_motif_residues)
                
                if non_motif_fixed:
                    structural_results['non_motif_freedom_check']['status'] = 'warning'
                    structural_results['non_motif_freedom_check']['details'].append(
                        f"Non-motif residues incorrectly FIXED: {sorted(non_motif_fixed)}"
                    )
                else:
                    structural_results['non_motif_freedom_check']['details'].append(
                        f"All {len(non_motif_residues)} non-motif residues are properly free"
                    )
            
            # 3. Coordinate validation (if motif PDB provided)
            if motif_coords and 'region_offsets' in trb_mapping:
                coord_issues = []
                
                # Load scaffold coordinates for comparison
                scaffold_coords = {}
                try:
                    with open(pdb_file, 'r') as f:
                        for line in f:
                            if line.startswith('ATOM'):
                                res_num = int(line[22:26].strip())
                                atom_name = line[12:16].strip()
                                x = float(line[30:38])
                                y = float(line[38:46])
                                z = float(line[46:54])
                                
                                if res_num not in scaffold_coords:
                                    scaffold_coords[res_num] = {}
                                scaffold_coords[res_num][atom_name] = (x, y, z)
                
                    # Compare key atoms (CA) for FIXED residues
                    for orig_res in expected_fixed:
                        if orig_res in trb_mapping['region_offsets']:
                            offset = trb_mapping['region_offsets'][orig_res]
                            scaffold_res = orig_res + offset
                            
                            # Check if coordinates exist in both
                            if ('A' in motif_coords and orig_res in motif_coords['A'] and 
                                'CA' in motif_coords['A'][orig_res] and
                                scaffold_res in scaffold_coords and 
                                'CA' in scaffold_coords[scaffold_res]):
                                
                                motif_ca = motif_coords['A'][orig_res]['CA']
                                scaffold_ca = scaffold_coords[scaffold_res]['CA']
                                
                                # Calculate distance (simple check for major displacement)
                                dist = ((motif_ca[0] - scaffold_ca[0])**2 + 
                                       (motif_ca[1] - scaffold_ca[1])**2 + 
                                       (motif_ca[2] - scaffold_ca[2])**2)**0.5
                                
                                if dist > 50.0:  # Large displacement threshold
                                    coord_issues.append(f"Residue {orig_res}->{scaffold_res}: CA displaced by {dist:.1f}Å")
                
                except Exception as e:
                    coord_issues.append(f"Could not compare coordinates: {e}")
                
                if coord_issues:
                    structural_results['coordinate_validation']['status'] = 'warning'
                    structural_results['coordinate_validation']['details'].extend(coord_issues)
                else:
                    structural_results['coordinate_validation']['details'].append(
                        "Coordinate validation passed - no major displacements detected"
                    )
                    
        except Exception as e:
            for check in structural_results.values():
                check['status'] = 'error'
                check['details'] = [f"Structural validation error: {e}"]
        
        return structural_results
    
    def validate_single_file(self, pdb_file: str, trb_file: str, 
                           expected_fixed: Set[int], motif_coords: Dict = None) -> Dict:
        """Validate a single PDB/TRB file pair"""
        result = {
            'pdb_file': pdb_file,
            'trb_file': trb_file,
            'status': 'success',
            'errors': [],
            'warnings': [],
            'metrics': {},
            'structural_validation': {}
        }
        
        try:
            # Get actual FIXED residues from PDB
            actual_fixed = self.get_fixed_residues_from_pdb(pdb_file)
            
            # Parse TRB mapping
            trb_mapping = self.parse_trb_file(trb_file) if trb_file else {}
            
            # Map expected residues to scaffold positions using TRB offsets
            mapped_expected = set()
            
            if trb_mapping and 'region_offsets' in trb_mapping:
                region_offsets = trb_mapping['region_offsets']
                
                for original_res in expected_fixed:
                    if original_res in region_offsets:
                        offset = region_offsets[original_res]
                        scaffold_res = original_res + offset
                        mapped_expected.add(scaffold_res)
                    else:
                        result['warnings'].append(f"Residue {original_res} not found in TRB mapping")
            else:
                # Fallback: if no TRB mapping, use original positions (for debugging)
                mapped_expected = expected_fixed.copy()
                result['warnings'].append("No TRB mapping available - using original residue numbers")
            
            # Calculate validation metrics
            correct_fixed = actual_fixed.intersection(mapped_expected)
            missing_fixed = mapped_expected - actual_fixed
            extra_fixed = actual_fixed - mapped_expected
            
            # Store metrics
            result['metrics'] = {
                'expected_count': len(mapped_expected),
                'actual_count': len(actual_fixed),
                'correct_count': len(correct_fixed),
                'missing_count': len(missing_fixed),
                'extra_count': len(extra_fixed),
                'accuracy': len(correct_fixed) / len(mapped_expected) if mapped_expected else 0,
                'precision': len(correct_fixed) / len(actual_fixed) if actual_fixed else 0,
                'recall': len(correct_fixed) / len(mapped_expected) if mapped_expected else 0,
                'expected_residues': sorted(mapped_expected),
                'actual_residues': sorted(actual_fixed),
                'correct_residues': sorted(correct_fixed),
                'missing_residues': sorted(missing_fixed),
                'extra_residues': sorted(extra_fixed)
            }
            
            # Perform structural validation
            result['structural_validation'] = self.validate_structural_consistency(
                pdb_file, trb_mapping, actual_fixed, mapped_expected, motif_coords or {}
            )
            
            # Check for issues
            if missing_fixed:
                result['warnings'].append(f"Missing {len(missing_fixed)} expected FIXED residues: {sorted(missing_fixed)}")
                
            if extra_fixed:
                result['warnings'].append(f"Found {len(extra_fixed)} unexpected FIXED residues: {sorted(extra_fixed)}")
                
            if result['metrics']['accuracy'] < 0.95 and mapped_expected:
                result['warnings'].append(f"Low accuracy: {result['metrics']['accuracy']:.2%}")
            
            # Add structural validation warnings
            for check_name, check_result in result['structural_validation'].items():
                if check_result['status'] in ['warning', 'error']:
                    for detail in check_result['details']:
                        result['warnings'].append(f"Structural {check_name}: {detail}")
                
        except Exception as e:
            result['status'] = 'error'
            result['errors'].append(str(e))
            self.log(f"Error validating {pdb_file}: {e}", "ERROR")
            
        return result
    
    def validate_pipeline(self, pdb_dir: str, trb_dir: str, 
                         contacts_file: str, motif_pdb_file: str = None) -> Dict:
        """Validate the entire pipeline results"""
        self.log("Starting pipeline validation...")
        
        # Load contacts data
        contacts_data = self.load_contacts_data(contacts_file)
        expected_fixed = self.get_expected_fixed_residues(contacts_data, self.distance_cutoff)
        
        # Load motif coordinates if provided
        motif_coords = {}
        if motif_pdb_file:
            motif_coords = self.load_original_motif_coordinates(motif_pdb_file)
        
        # Find matching PDB and TRB files
        pdb_files = list(Path(pdb_dir).glob("*.pdb"))
        trb_files = {f.stem: str(f) for f in Path(trb_dir).glob("*.trb")}
        
        self.log(f"Found {len(pdb_files)} PDB files and {len(trb_files)} TRB files")
        
        # Validate each file
        file_results = []
        summary_stats = defaultdict(int)
        structural_summary = defaultdict(int)
        
        for pdb_file in pdb_files:
            pdb_stem = pdb_file.stem
            trb_file = trb_files.get(pdb_stem)
            
            if not trb_file:
                self.log(f"Warning: No matching TRB file for {pdb_file.name}", "WARNING")
                
            result = self.validate_single_file(str(pdb_file), trb_file, expected_fixed, motif_coords)
            file_results.append(result)
            
            # Update summary statistics
            summary_stats['total_files'] += 1
            summary_stats[result['status']] += 1
            
            # Update structural validation summary
            if 'structural_validation' in result:
                for check_name, check_result in result['structural_validation'].items():
                    structural_summary[f'{check_name}_{check_result["status"]}'] += 1
            
        # Calculate averages for successful validations
        success_count = summary_stats['success']
        if success_count > 0:
            metric_totals = defaultdict(float)
            
            # Sum up all metrics from successful results
            for result in file_results:
                if result['status'] == 'success' and 'metrics' in result:
                    for metric_name, metric_value in result['metrics'].items():
                        if isinstance(metric_value, (int, float)):
                            metric_totals[f'avg_{metric_name}'] += metric_value
            
            # Calculate averages
            for key, total in metric_totals.items():
                summary_stats[key] = total / success_count
                    
        # Store results
        self.validation_results = {
            'summary': dict(summary_stats),
            'structural_summary': dict(structural_summary),
            'file_results': file_results,
            'pipeline_config': {
                'distance_cutoff': self.distance_cutoff,
                'expected_fixed_count': len(expected_fixed),
                'expected_fixed_residues': sorted(expected_fixed),
                'motif_pdb_provided': bool(motif_pdb_file),
                'structural_validation_enabled': bool(motif_coords)
            }
        }
        
        return self.validation_results
    
    def generate_report(self, output_file: str = None) -> str:
        """Generate a comprehensive validation report"""
        results = self.validation_results
        
        report_lines = [
            "=" * 80,
            "SELECTIVE FIXED LABELS PIPELINE VALIDATION REPORT",
            "=" * 80,
            "",
            "SUMMARY STATISTICS:",
            f"  Total files processed: {results['summary'].get('total_files', 0)}",
            f"  Successful validations: {results['summary'].get('success', 0)}",
            f"  Failed validations: {results['summary'].get('error', 0)}",
            f"  Distance cutoff used: {results['pipeline_config']['distance_cutoff']}",
            f"  Expected FIXED residues: {results['pipeline_config']['expected_fixed_count']}",
            "",
            "AVERAGE METRICS:",
            f"  Accuracy: {results['summary'].get('avg_accuracy', 0):.2%}",
            f"  Precision: {results['summary'].get('avg_precision', 0):.2%}",
            f"  Recall: {results['summary'].get('avg_recall', 0):.2%}",
            f"  Expected count: {results['summary'].get('avg_expected_count', 0):.1f}",
            f"  Actual count: {results['summary'].get('avg_actual_count', 0):.1f}",
        ]
        
        # Add structural validation summary if available
        if 'structural_summary' in results and results['structural_summary']:
            report_lines.extend([
                "",
                "STRUCTURAL VALIDATION SUMMARY:",
            ])
            
            structural_summary = results['structural_summary']
            total_files = results['summary'].get('success', 0)
            
            for check_type in ['motif_region_check', 'non_motif_freedom_check', 'coordinate_validation']:
                pass_count = structural_summary.get(f'{check_type}_pass', 0)
                warn_count = structural_summary.get(f'{check_type}_warning', 0)
                error_count = structural_summary.get(f'{check_type}_error', 0)
                
                check_name = check_type.replace('_', ' ').title()
                report_lines.append(f"  {check_name}:")
                report_lines.append(f"    Pass: {pass_count}/{total_files} ({pass_count/total_files*100:.1f}%)")
                if warn_count > 0:
                    report_lines.append(f"    Warnings: {warn_count}")
                if error_count > 0:
                    report_lines.append(f"    Errors: {error_count}")
        
        report_lines.extend([
            "",
            "DETAILED RESULTS:",
        ])
        
        # Add per-file details for files with issues
        issues_found = False
        for result in results['file_results']:
            if result['errors'] or result['warnings']:
                if not issues_found:
                    issues_found = True
                    report_lines.append("")
                    
                report_lines.extend([
                    f"\nFile: {Path(result['pdb_file']).name}",
                    f"  Status: {result['status']}",
                ])
                
                if result['errors']:
                    report_lines.extend([f"  ERROR: {error}" for error in result['errors']])
                    
                if result['warnings']:
                    # Group structural warnings separately
                    structural_warnings = [w for w in result['warnings'] if 'Structural' in w]
                    other_warnings = [w for w in result['warnings'] if 'Structural' not in w]
                    
                    for warning in other_warnings:
                        report_lines.append(f"  WARNING: {warning}")
                    
                    if structural_warnings:
                        report_lines.append("  STRUCTURAL ISSUES:")
                        for warning in structural_warnings:
                            report_lines.append(f"    {warning}")
                    
                if 'metrics' in result:
                    report_lines.extend([
                        f"  Accuracy: {result['metrics'].get('accuracy', 0):.2%}",
                        f"  Expected/Actual: {result['metrics'].get('expected_count', 0)}/{result['metrics'].get('actual_count', 0)}"
                    ])
        
        if not issues_found:
            report_lines.append("\nNo issues detected in any files - all validations passed!")
        
        report_lines.extend([
            "",
            "=" * 80,
            f"Expected FIXED residues: {results['pipeline_config']['expected_fixed_residues']}",
            "",
            "STRUCTURAL VALIDATION:",
            f"  Motif PDB provided: {'Yes' if results['pipeline_config'].get('motif_pdb_provided', False) else 'No'}",
            f"  Coordinate validation: {'Enabled' if results['pipeline_config'].get('structural_validation_enabled', False) else 'Disabled'}",
            "=" * 80
        ])
        
        report_text = "\n".join(report_lines)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
            self.log(f"Report saved to {output_file}")
            
        return report_text

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Validate Selective FIXED Labels Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument("--pdb_dir", required=True,
                       help="Directory containing processed PDB files")
    parser.add_argument("--trb_dir", required=True,
                       help="Directory containing TRB files")
    parser.add_argument("--contacts_file", required=True,
                       help="JSON file with motif contact analysis results")
    parser.add_argument("--motif_pdb",
                       help="Original motif PDB file for structural validation")
    parser.add_argument("--distance_cutoff", type=float, default=4.0,
                       help="Distance cutoff used in pipeline (default: 4.0)")
    parser.add_argument("--output_report", 
                       help="Output file for text validation report")
    parser.add_argument("--output_json",
                       help="Output JSON file for detailed results (auto-generated if not specified)")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Validate input paths
    for path, name in [(args.pdb_dir, "PDB directory"), 
                       (args.trb_dir, "TRB directory"),
                       (args.contacts_file, "Contacts file")]:
        if not os.path.exists(path):
            print(f"Error: {name} not found: {path}")
            return 1
    
    if args.motif_pdb and not os.path.exists(args.motif_pdb):
        print(f"Error: Motif PDB file not found: {args.motif_pdb}")
        return 1
    
    try:
        # Run validation
        validator = FixedLabelValidator(
            distance_cutoff=args.distance_cutoff,
            verbose=args.verbose
        )
        
        results = validator.validate_pipeline(
            args.pdb_dir, 
            args.trb_dir, 
            args.contacts_file,
            args.motif_pdb
        )
        
        # Generate report
        if args.output_report:
            if args.output_report.endswith('.json'):
                # User specified JSON extension for report, use as JSON output
                json_file = args.output_report
                report_file = args.output_report.replace('.json', '.txt')
            else:
                report_file = args.output_report
                json_file = args.output_json or args.output_report.replace('.txt', '.json')
        else:
            report_file = "validation_report.txt"
            json_file = args.output_json or "validation_results.json"
            
        report_text = validator.generate_report(report_file)
        
        if args.verbose:
            print("\n" + report_text)
        
        # Save JSON results
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Detailed results saved to {json_file}")
        print(f"Text report saved to {report_file}")
            
        # Print summary
        summary = results['summary']
        print(f"\nValidation Complete:")
        print(f"  Files processed: {summary.get('total_files', 0)}")
        print(f"  Success rate: {summary.get('success', 0)}/{summary.get('total_files', 0)}")
        print(f"  Average accuracy: {summary.get('avg_accuracy', 0):.2%}")
        
        return 0
        
    except Exception as e:
        print(f"Validation failed: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
