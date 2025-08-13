awk '
BEGIN {
    FS = "\t";          # Set input field separator to tab
    OFS = "\t";         # Set output field separator to tab

    # Load all gene IDs from all.id.txt into associative array 'ids'
    while (getline < "all.id.txt") { #all.id.txt is the NLR gene id file.
        ids[$1] = 1;    # Use first column (gene ID) as key
    }
    close("all.id.txt"); # Close file after reading
}
{
    if (FNR == 1) {
        # Process header line (first line of the file)

        # Print first 4 columns unchanged
        for (i = 1; i <= 4; i++) {
            printf "%s%s", $i, OFS;
        }

        # For columns 5 to end: check if any gene in the comma-separated list
        # appears in the whitelist (ids). If yes, keep the original column name;
        # otherwise, output empty string.
        for (i = 5; i <= NF; i++) {
            split($i, genes, ",");   # Split gene list by comma
            found = 0;
            for (j in genes) {
                if (genes[j] in ids) {
                    found = 1;        # Mark as found if any gene matches
                    break;
                }
            }
            printf "%s", found ? $i : "";  # Keep original header if matched
            if (i < NF) printf OFS;        # Add OFS unless it's the last field
        }
        printf "\n";
    } else {
        # Process data lines (non-header rows)

        # Print first 4 columns unchanged
        for (i = 1; i <= 4; i++) {
            printf "%s%s", $i, OFS;
        }

        # For columns 5 to end: filter gene lists to keep only those in the whitelist
        for (i = 5; i <= NF; i++) {
            split($i, genes, ",");   # Split gene list by comma
            output = "";              # Initialize filtered gene list
            for (j in genes) {
                if (genes[j] in ids) {
                    if (output != "") output = output ",";
                    output = output genes[j];  # Append matching gene
                }
            }
            printf "%s", output;      # Output filtered gene list (or empty)
            if (i < NF) printf OFS;  # Add separator if not last field
        }
        printf "\n";
    }
}' sample.list.SG.pan > NLR.sample.list.SG.pan	#sample.list.SG.pan file is the whole pan-genome list of all proteins of samples, NLR.sample.list.SG.pan file is the result file of NLR-pangenome list.
