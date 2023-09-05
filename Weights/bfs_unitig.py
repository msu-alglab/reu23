import argparse
from bfs_uni import *

parser = argparse.ArgumentParser(
    description="""
    get our edgeweights!!
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-q', '--seq', type=str, help='seq file', required=True)
requiredNamed.add_argument('-g', '--seg', type=str, help='seg file', required=True)
requiredNamed.add_argument('-m', '--sam', type=str, help='sam file', required=True)
requiredNamed.add_argument('-k', '--kval', type=int, help='k value', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

args = parser.parse_args()

vertices = readseg(args.seg)
edges, source, sink, longestpath = readseq(args.seq, vertices, args.kval)
readsam(args.sam, args.kval, edges, vertices)

print_counts(args.output, edges, vertices, source, sink)
#counts = edgewt_bfs(edges, vertices, source, args.kval)
#output_edges(args.output, counts, vertices)
#output_inexact_edges(args.output, counts, vertices)