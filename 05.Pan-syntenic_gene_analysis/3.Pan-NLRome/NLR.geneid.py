import os
import re
import shutil
from collections import defaultdict

def setup_directories(base_dir):
    """Create all required directory structures"""
    dirs = {
        'nlr_ann': os.path.join(base_dir, "0.NLR_Ann"),
        'gff3': os.path.join(base_dir, "0.gff3"),
        'nlr_id': os.path.join(base_dir, "1.NLR_id"),
        'nlr_id_sort': os.path.join(base_dir, "1.NLR_id_sort"),
    }
    
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
        print(f"Ensuring directory exists: {dir_path}")
    
    return dirs

def log_step(step, message):
    """Log processing steps"""
    print(f"\n{'='*50}")
    print(f"Step {step}: {message}")
    print(f"{'='*50}")

def parse_nlr_gff(nlr_gff_file):
    """Parse NLR GFF file"""
    nlr_data = []
    try:
        with open(nlr_gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                match = re.search(r'nlrClass=([^;]+)', attributes)
                nlr_class = match.group(1) if match else "Unknown"
                nlr_data.append((chrom, start, end, nlr_class))
        print(f"Successfully parsed {len(nlr_data)} records from {nlr_gff_file}")
    except Exception as e:
        print(f"Error parsing {nlr_gff_file}: {str(e)}")
    return nlr_data

def parse_species_gff(species_gff_file):
    """Parse species GFF3 file"""
    species_data = {}
    try:
        with open(species_gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene':
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                match = re.search(r'ID=(?:gene:)?([^;]+)', attributes)
                gene_id = match.group(1) if match else "Unknown"
                if chrom not in species_data:
                    species_data[chrom] = []
                species_data[chrom].append((start, end, gene_id))
        print(f"Successfully parsed {sum(len(v) for v in species_data.values())} gene records from {species_gff_file}")
    except Exception as e:
        print(f"Error parsing {species_gff_file}: {str(e)}")
    return species_data

def find_matches(nlr_data, species_data):
    """Find matching genes"""
    results = []
    for chrom, start, end, nlr_class in nlr_data:
        if chrom in species_data:
            for mrna_start, mrna_end, gene_id in species_data[chrom]:
                if not (end < mrna_start or start > mrna_end):
                    results.append((gene_id, start, end, nlr_class))
    print(f"Found {len(results)} matching genes")
    return results

def process_gene_matching(dirs, species_list):
    """Process gene matching step"""
    log_step(1, "Processing gene matching")
    
    for species in species_list:
        print(f"\nProcessing species: {species}")
        nlr_gff = os.path.join(dirs['nlr_ann'], f"{species}.NLR.gff")
        species_gff = os.path.join(dirs['gff3'], f"{species}.HC.gff3")
        output_file = os.path.join(dirs['nlr_id'], f"{species}.idgene.txt")
        
        if not os.path.exists(nlr_gff):
            print(f"Warning: NLR GFF file not found: {nlr_gff}")
            continue
        if not os.path.exists(species_gff):
            print(f"Warning: Species GFF file not found: {species_gff}")
            continue
        
        nlr_data = parse_nlr_gff(nlr_gff)
        species_data = parse_species_gff(species_gff)
        matches = find_matches(nlr_data, species_data)
        
        with open(output_file, 'w') as f:
            f.write("Gene_ID\tStart\tEnd\tNLR_Class\n")
            for gene_id, start, end, nlr_class in matches:
                f.write(f"{gene_id}\t{start}\t{end}\t{nlr_class}\n")
        print(f"Match results written to: {output_file}")

def main():
    # Base directory
    base_dir = "./"
    
    # Setup directory structure
    dirs = setup_directories(base_dir)
    
    # Species list
    species_list = ["sample1","sample2"]
    
    # Execute processing pipeline
    process_gene_matching(dirs, species_list)
    
    print("\nAll processing steps completed!")

if __name__ == "__main__":
    main()
