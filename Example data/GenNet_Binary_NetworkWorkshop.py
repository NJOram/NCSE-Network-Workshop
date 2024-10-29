import numpy as np
import networkx as nx
import pandas as pd
from NEMtropy import UndirectedGraph, matrix_generator
from NEMtropy.network_functions import build_adjacency_from_edgelist



df = pd.read_csv('Conventional network edgelist.csv')


G = nx.from_pandas_edgelist(df,source='V1',target='V2')
G = G.to_undirected()

# # Transform into adjacency matrix
adj_bin = nx.adjacency_matrix(G).toarray()

dseq= np.sum(adj_bin, axis=1)


graph = UndirectedGraph(adj_bin)

graph.solve_tool(model="cm_exp",
                 method="newton",
                 initial_guess="random")



graph.ensemble_sampler(999, cpu_n=1, output_dir="Conventional_null_ensemble/")

