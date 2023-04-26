#!python3

import argparse
import sys
import pandas as pd
import bte

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "--mat_pb", help="MAT protobuf file", type=str, 
    default="mutation_annotated_tree.pb")
parser.add_argument(
    "--introductions", help="Introductions file", type=str,
    default="introductions.tsv")
args = parser.parse_args()

# Load data
mat = bte.MATree(args.mat_pb)
introductions = pd.read_table(args.introductions)

# Get introduction nodes
introduction_nodes = set(introductions['introduction_node'].to_list())

# Get the leaves for each introduction node in tree
leaf_data = pd.DataFrame(columns=['introduction_node', 'leaf'])

for introduction_node in introduction_nodes:
    node_id = "node_" + introduction_node.split("_node_", 1)[1]
    leaf_list = mat.get_leaves_ids(nid=node_id)
    print(str(introduction_node) + " has " + str(len(leaf_list)) + " leaves")
    list_of_tuples = list(zip([introduction_node] * len(leaf_list), leaf_list))
    leaf_data_tmp = pd.DataFrame(list_of_tuples, columns = ['introduction_node', 'leaf'])
    leaf_data = pd.concat([leaf_data, leaf_data_tmp], ignore_index=True)

leaf_data.to_csv("node_to_leaves.csv", index=False)