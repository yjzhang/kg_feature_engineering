import kgfe

df = kgfe.load_graph('kegg_pathway_data.csv')
graph = kgfe.df_to_networkx(df)
results = kgfe.explanations.hypgergeom_test(graph, [7124, 3667, 3630], 'Gene')
print(results)
