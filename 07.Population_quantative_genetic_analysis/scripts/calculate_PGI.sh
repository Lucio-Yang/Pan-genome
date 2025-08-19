#!/bin/bash
# calculate_pgi.sh
# Compute PGI = sum((end-start)*fd) / genome_length for all CSV files

INPUT_DIR="output_2_csv_cat"     # Input folder with CSV files
OUTPUT_FILE="pgi_summary.csv"   # Output summary file
GENOME_LENGTH=10422265814        # Wheat genome effective length (bp)

# Write header
echo "File,P1,P2,Target,PGI" > "$OUTPUT_FILE"

# Process each CSV file
for file in "$INPUT_DIR"/*.csv; do
    [[ ! -f "$file" ]] && continue
    filename=$(basename "$file")

    # Extract P1_P2_Target from filename (e.g., DEW_WEW_Ethiopia)
    if [[ "$filename" =~ \.([^.]+_.*_.*?)\.w500m50\.csv ]]; then
        IFS='_' read -r P1 P2 Target <<< "${BASH_REMATCH[1]}"
    else
        P1="Unknown"; P2="Unknown"; Target="Unknown"
    fi

    # Calculate PGI using awk
    pgi=$(awk -F',' '
    NR==1 {
        # Get column indices
        for(i=1; i<=NF; i++) {
            gsub(/[^a-zA-Z0-9_]/, "", $i)
            if($i == "start") s = i
            if($i == "end")   e = i
            if($i == "fd")    f = i
        }
        next
    }
    {
        # Skip invalid rows
        if ($s == "" || $e == "" || $f == "") next
        if ($f < 0 || $f > 1) next
        if ($e <= $s) next
        total += ($e - $s) * $f
    }
    END {
        printf "%.6f", total / ENVIRON["GENOME_LENGTH"]
    }' "GENOME_LENGTH=$GENOME_LENGTH" "$file")

    # Save result
    echo "$filename,$P1,$P2,$Target,$pgi" >> "$OUTPUT_FILE"
    echo "Processed: $filename -> PGI=$pgi"
done

echo "PGI calculation done. Output: $OUTPUT_FILE"
