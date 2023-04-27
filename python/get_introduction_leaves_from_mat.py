# Local: conda activate cluster-detection
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import toytree as ttplt
import bte

# Load data
data_dir = "../nextflow/results/urmc_klebsiella_2023-04-21T10:03:44.801358-07:00/"
mat = bte.MATree(data_dir + "matutils/mutation_annotated_tree.pb")
introductions = pd.read_table(data_dir + "clustertracker/introductions.tsv")

# Get introduction nodes
introduction_nodes = []
for index, row in introductions.iterrows():
    introduction_nodes.append(row['introduction_node'])

# Get the leaves for each introduction node in tree
leaf_data = pd.DataFrame.from_dict({
    'introduction_node': [],
    'leaf': []
})

for introduction_node in introduction_nodes:
    node_id = introduction_node.split("_", 1)[1]
    leaf_list = mat.get_leaves_ids(nid=node_id)
    leaf_data_tmp = pd.DataFrame.from_dict({'introduction_node':introduction_node, 'leaf':leaf_list})
    leaf_data = pd.concat([leaf_data, leaf_data_tmp], ignore_index=True)

leaf_data.to_csv(data_dir + "clustertracker/node_to_leaves.csv", index=False)