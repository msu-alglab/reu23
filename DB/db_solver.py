import argparse
import time
from db import *

parser = argparse.ArgumentParser(
    description="""
    Build a DeBruijn graph from FASTQ data.
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )

parser.add_argument('-k', '--kval', type=int, default=10, help='val of k for k-mer:\n')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--input', type=str, help='Input filename', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)

args = parser.parse_args()
start = time.time()
#reads = read_reads(args.input)
reads = read_map(args.input, args.kval)

#g = construct_graph(reads, args.kval)
g = map_constructor(reads, args.kval)
#g = source_sink(g)
end = time.time()

#print_reads(args.output, reads, args.kval)
print_reset(args.output)
print_runtime(args.output, (end - start))
print_graph(args.output, g)

