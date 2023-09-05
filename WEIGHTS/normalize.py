import argparse
import time
from norm import *

parser = argparse.ArgumentParser(
    description="""
    get our edgeweights!!
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-s', '--sam', type=str, help='sam file', required=True)
requiredNamed.add_argument('-e', '--edges', type=str, help='edges', required=True)
requiredNamed.add_argument('-k', '--kval', type=int, help='k value', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

args = parser.parse_args()

kmers, reads = readsam(args.sam, args.kval)
vertices = readvertices(args.edges, args.kval)
iterate(reads, vertices, args.kval)
bestfit(vertices)
normalize(vertices, kmers)
output_edges(args.output, vertices, args.kval)