#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : zhangyan
# @Email         : zhangyan1220zj@163.com
# @FileName      : MST-Kruskal.py
# @Time          : 2024-09-25 20:30:23
# @description   : Analysis of Minimum Spanning Tree using Kruskal's Algorithm
"""

import os
import glob
from Bio import Phylo  # Biopython module for handling phylogenetic trees
import networkx as nx  # Library for creating and manipulating complex networks/graphs


def generate_mst_from_tree(input_file, output_file):

    tree = Phylo.read(input_file, "newick")
    leaves = tree.get_terminals()
    leaf_names = [leaf.name for leaf in leaves]
    
    G = nx.Graph()
    
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            distance = tree.distance(leaves[i], leaves[j])
            G.add_edge(leaf_names[i], leaf_names[j], weight=distance)
    
    # Compute the Minimum Spanning Tree using Kruskal's algorithm
    mst = nx.minimum_spanning_tree(G, algorithm='kruskal')
    
    mst_edges = list(mst.edges(data=True))
    
    
    with open(output_file, "a") as f:
        f.write(f"{os.path.basename(input_file)}\t")
        f.write(f"{mst_edges}\n")


if __name__ == "__main__":
    input_folder = "tree"  # Modify: change to your input folder name

    # Define the output file to store MST edges from all trees
    output_file = "mst_edges.txt"  # Modify: change to your desired output filename

    # Find all tree files with extension '.out' in the input folder
    tree_files = glob.glob(os.path.join(input_folder, "*.out"))

    for tree_file in tree_files:
        generate_mst_from_tree(tree_file, output_file)
    print(f"All MST edges have been written to {output_file}")
