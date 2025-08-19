#!/bin/bash


# Step 1: Extract OG IDs where column 2 == 13 and column 3 == 13
awk -F'\t' '$2 == 13 && $3 == 13 {print $1}' sample.list.SG.pan > single-copy_og_list.txt

# Step 2: Merge all FASTA files with species prefix (e.g., >SpeciesA|gene1)
> all_proteins.fasta
for file in fasta/*; do
    [[ -f "$file" ]] || continue
    species=$(basename "$file" | sed 's/\.[^.]*$//')  # Remove extension
    cat "$file" | sed "s/^>/>${species}|/g" >> all_proteins.fasta
done

# Step 3: Create directories for outputs
mkdir -p OG_alignments alignments

# Step 4: Extract sequences for each OG into separate FASTA files
while IFS= read -r og_id; do
    # Get gene list from Orthogroups.tsv (skip first column: OG ID)
    genes=$(grep "^${og_id} " Orthogroups.tsv | cut -f2- | tr '\t,' '\n' | grep -v "^$")
    
    if [ -z "$genes" ]; then
        echo "Warning: $og_id not found in Orthogroups.tsv"
        continue
    fi

    # Output unaligned file
    output_fasta="OG_alignments/${og_id}.fa"
    > "$output_fasta"

    # Extract each gene from merged FASTA
    echo "$genes" | while read gene; do
        sed -n "/^>${gene}\$/,/^\s*>/p" all_proteins.fasta | \
        sed '/^\s*>/d; /^$/d'  # Remove next header and blank lines
    done >> "$output_fasta"

    # Step 5: Run MAFFT alignment (default parameters) on this OG
    if [ -s "$output_fasta" ]; then
        mafft --auto "$output_fasta" > "alignments/${og_id}.aln.fa" 2>/dev/null
        echo "Aligned: ${og_id}"
    else
        echo "Warning: No sequences for $og_id, skipping alignment"
    fi

done < single-copy_og_list.txt

# Step 6: (Optional) Concatenate all aligned files into one
cat alignments/*.aln.fa > alignments/All_SingleCopy_Aligned.fasta
