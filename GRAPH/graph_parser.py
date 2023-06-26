'''
takes in an edge list, creates a graph, runs dfs, outputs relevant info:
1. acyclic?
2. # ccs
3. size of each cc
4. ...
'''

import argparse
from graph import *
import sys

sys.setrecursionlimit(10**6)

parser = argparse.ArgumentParser(
    description="""
    parse edgelist....
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

args = parser.parse_args()
edges = read_input(args.input)

graphs = construct_graph(edges)
g = graphs[0]
u = graphs[1]
u = DFS(u)
g = DFS(g)

#print_reads(args.output, reads, args.kval)
print_reset(args.output)
print_info(args.output, g, u)
print_cc_dummycheck(args.output, u)

