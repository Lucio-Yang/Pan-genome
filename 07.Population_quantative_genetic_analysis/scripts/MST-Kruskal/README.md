# MST from Phylogenetic Trees

Generate Minimum Spanning Trees (MST) from phylogenetic trees using Kruskal's algorithm.

## Requirements

pip install biopython networkx

## Run

python MST-Kruskal.py

## Directory Structure

```text
project/
├── tree/					# Input folder containing tree files (.out format)
│   ├── tree1.out			
│   ├── tree2.out
│   └── ...
├── MST-Kruskal.py			# Main script to process tree files and extract MST edges
├── mst_edges.txt			# Output file storing the extracted minimum spanning tree edges
└── README.md				# Project documentation and instructions
```

