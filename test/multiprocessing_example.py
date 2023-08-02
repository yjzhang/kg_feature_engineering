# basic usage example for kgfe
import os
import time

import kgfe

# timing ran on intel i7-9700 CPU @ 3.00GHz with 64gb memory

df = kgfe.load_graph('reactome_genes_chems.csv.gz')
t = time.time()
graph = kgfe.df_to_graph(df)
graph.simplify()
print('igraph graph from df time:', time.time() - t)
# timing: 12s
topic_ids = ['NCBIGene::5972',
        'NCBIGene::958',
        'NCBIGene::100', 'NCBIGene::8797', 'NCBIGene::26762']

# TODO: null model

import multiprocessing as mp

def null_stats_process(n_ids, n_samples, q):
    output = kgfe.explanations.null_graph_stats(graph, 'Gene', n_ids, n_samples)
    q.put(output)

def run_mp(n_procs=0, n_ids=10, n_samples=100):
    if n_procs == 0:
        n_procs = os.cpu_count()
    print(n_procs)
    n_samples_per_proc = int(n_samples/n_procs)
    print(n_samples_per_proc)
    q = mp.Queue()
    processes = []
    #pool = mp.Pool(processes=n_procs)
    #pool.imap()
    for _ in range(n_procs):
        p = mp.Process(target=null_stats_process, args=(n_ids, n_samples_per_proc, q,))
        processes.append(p)
        p.start()
    results = []
    for p in processes:
        p.join()
        results.extend(q.get())
    return results


t = time.time()
result = run_mp(n_ids=20, n_samples=100)
print('time for multiprocessing null graph stats:', time.time() - t)
print('num samples:', len(result))

t = time.time()
null_stats = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100)
print('time for single-process null graph stats:', time.time() - t)

t = time.time()
null_stats_p = kgfe.explanations.null_graph_stats(graph, 'Gene', 20, 100, parallel=True)
print('time for multiprocessing pool null graph stats:', time.time() - t)
print(len(null_stats_p))

