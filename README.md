# selective-fixed-labels
Automated pipeline for selective FIXED labels in ProteinMPNN motif design. Analyzes protein contacts, parses RFdiffusion TRB files, and generates scripts that fix only critical residues instead of entire motif regions. Improves design flexibility and success rates in complex scaffolding projects.

# Selective FIXED Labels Pipeline for ProteinMPNN

![Pipeline](https://img.shields.io/badge/Pipeline-ProteinMPNN-blue)
![Python](https://img.shields.io/badge/Python-3.8%2B-green)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Status](https://img.shields.io/badge/Status-Active-brightgreen)

A streamlined toolkit for applying selective FIXED labels in ProteinMPNN motif-based protein design workflows. This pipeline automates the critical step of identifying which motif residues should remain fixed during sequence design, based on contact analysis and structural mapping.

## 🎯 What This Does

Instead of fixing entire motif regions (which can be overly restrictive), this approach selectively fixes only the residues that form critical contacts with the target protein, allowing for more flexible and successful design outcomes. The toolkit consists of four Python scripts that work together seamlessly with existing RFdiffusion → ProteinMPNN workflows.

## 📋 Requirements

- Python 3.8+
- BioPython (optional, for enhanced validation)
- Standard libraries: argparse, json, pickle, collections, pathlib, datetime, math

## 🚀 Quick Start

```bash
# 1. Analyze motif contacts
python motif_contact_analysis_hpc.py --pdb your_motif.pdb --output contacts.json

# 2. Parse TRB configuration  
python parse_trb_config.py --trb_file your_file.trb --output_json config.json

# 3. Generate selective FIXED script
python auto_config_fixed_script.py --trb_config config.json --contact_file contacts.json --output selective_fixed.py

# 4. Apply selective FIXED labels
python selective_fixed.py --pdbdir /path/to/scaffolds --verbose

# 5. Validate results (optional)
python validate_selective_fixed_pipeline.py --pdb_dir /path/to/processed --trb_dir /path/to/trb --contacts_file contacts.json
```

## 📖 Script Documentation

### 1. `motif_contact_analysis_hpc.py`
**Purpose:** Analyzes contacts between motif and target chains in PDB files

**Basic Usage:**
```bash
python motif_contact_analysis_hpc.py --pdb motif.pdb
```

**Advanced Usage:**
```bash
python motif_contact_analysis_hpc.py --pdb motif.pdb \
  --motif_chain A --motif_residues "66-93,109-127" \
  --target_chain B --distance_cutoffs "3.5,4.0,4.5" \
  --output contacts.json --verbose
```

**Key Options:**
- `--pdb`: Input PDB file with motif and target chains
- `--motif_residues`: Specify motif regions (e.g., "66-93,109-127" for dual motif)
- `--distance_cutoffs`: Contact distance thresholds (default: 3.0,3.5,4.0,4.5,5.0)
- `--verbose`: Detailed output showing contact analysis

**Output:** JSON file with contact analysis results for different distance cutoffs

---

### 2. `parse_trb_config.py`
**Purpose:** Extracts motif mapping information from RFdiffusion TRB files

**Basic Usage:**
```bash
python parse_trb_config.py --trb_file example.trb --output_json config.json
```

**Advanced Usage:**
```bash
python parse_trb_config.py --trb_file example.trb --verbose
```

**Key Options:**
- `--trb_file`: Input TRB file from RFdiffusion
- `--output_json`: Output JSON configuration file
- `--verbose`: Show detailed mapping information

**Output:** JSON file containing:
- Original motif residue ranges
- Scaffold position mappings  
- Offset calculations for each region
- Contig mapping strings

---

### 3. `auto_config_fixed_script.py`
**Purpose:** Generates customized selective FIXED label scripts

**Basic Usage:**
```bash
python auto_config_fixed_script.py --trb_config config.json --contact_file contacts.json
```

**Manual Contact Specification:**
```bash
python auto_config_fixed_script.py --trb_config config.json \
  --contact_residues "70,71,72,73" --distance_cutoff 4.0
```

**Key Options:**
- `--trb_config`: JSON config from parse_trb_config.py
- `--contact_file`: Contact analysis JSON file
- `--contact_residues`: Manual specification of contact residues
- `--distance_cutoff`: Distance threshold (3.5, 4.0, or 4.5 Å)
- `--extra_residues`: Additional residues to always fix
- `--output`: Output script filename

**Output:** Executable Python script that applies selective FIXED labels

---

### 4. `validate_selective_fixed_pipeline.py`
**Purpose:** Validates that selective FIXED labeling worked correctly

**Basic Usage:**
```bash
python validate_selective_fixed_pipeline.py \
  --pdb_dir /path/to/processed/pdbs \
  --trb_dir /path/to/trb/files \
  --contacts_file contacts.json
```

**Full Validation:**
```bash
python validate_selective_fixed_pipeline.py \
  --pdb_dir /path/to/processed/pdbs \
  --trb_dir /path/to/trb/files \
  --contacts_file contacts.json \
  --motif_pdb original_motif.pdb \
  --distance_cutoff 4.0 \
  --output_report validation_report.txt \
  --verbose
```

**Key Options:**
- `--pdb_dir`: Directory with processed PDB files
- `--trb_dir`: Directory with TRB files
- `--contacts_file`: Original contact analysis JSON
- `--motif_pdb`: Original motif PDB for structural validation
- `--distance_cutoff`: Distance cutoff used in pipeline
- `--output_report`: Text report file
- `--verbose`: Detailed validation output

**Output:** Validation report with accuracy metrics and structural checks

## 📁 File Flow

```
Motif PDB ──→ [Script 1] ──→ contacts.json
                                     ↓
TRB file ───→ [Script 2] ──→ config.json ──→ [Script 3] ──→ selective_fixed.py
                                                                     ↓
Scaffold PDBs ←─────────────────────────────────────────────────────┘
     ↓
[Script 4] ──→ validation_report.txt
```

## ⚙️ Configuration

### Script 1 - Contact Analysis
Edit the top of `motif_contact_analysis_hpc.py`:
```python
DEFAULT_MOTIF_CHAIN = "A"
DEFAULT_TARGET_CHAIN = "B"  
DEFAULT_MOTIF_RESIDUES = "66-93,109-127"  # Change for your motif
DEFAULT_DISTANCE_CUTOFFS = "3.0,3.5,4.0,4.5,5.0"
```

### Script 3 - Script Generator
Edit the top of `auto_config_fixed_script.py`:
```python
DEFAULT_DISTANCE_CUTOFF = 4.0
EXTRA_FIXED_RESIDUES = []  # Add residue numbers to always fix
```

## 🔧 Typical Workflow

1. **Prepare your motif PDB** with chains A (motif) and B (target)
2. **Run contact analysis** to identify critical residues
3. **Parse your TRB file** from RFdiffusion scaffolding
4. **Generate selective FIXED script** combining contact + TRB data
5. **Apply to scaffold PDBs** using the generated script
6. **Validate results** to ensure accuracy
7. **Run ProteinMPNN** with selectively fixed scaffolds

## 🐛 Troubleshooting

**No contacts found:** Check motif/target chain assignments and try higher distance cutoffs

**TRB parsing fails:** Ensure TRB file is from RFdiffusion pipeline

**Script generation errors:** Verify contact residues fall within motif regions defined in TRB

**Validation warnings:** Review structural consistency and position mappings

## 📧 Contact

e.leen@leeds.ac.uk

## 📄 License

MIT License - see LICENSE file for details
