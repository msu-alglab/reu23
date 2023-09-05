import math
#takes in a seq file, a seg file, sam file, kval
#want to return maximal unitigs, with weights
#preferably return the ids and not the actual letters

class Read:
    def __init__(self, lab, seq):
        self.label = lab
        self.seq = seq

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.indegree = 0
        self.outdegree = 0
        self.outcount = 0
        self.neighbors = []

class Edge:
    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2
        self.weight = 0
        self.counts = 0

def reversecomplement(s):
    output = ''
    for i in range(1, len(s) + 1):
        index = -i
        if s[index] == 'A':
            output += ('T')
        if s[index] == 'T':
            output += ('A')
        if s[index] == 'C':
            output += ('G')
        if s[index] == 'G':
            output += ('C')
    return output

#reads in list of nodes from seg file
#returns vertex map
def readseg(fname):
    lines = open(fname, 'r').read().split('\n')
    lines.pop()
    vertices = dict()
    for line in lines:
        s = line.split()
        plus = s[0] + '+'
        minus = s[0] + '-'
        vertices[plus] = Vertex(s[1])
        vertices[minus] = Vertex(reversecomplement(s[1]))
    return vertices


#reads in all present edges from the seq file
#returns edge map
def readseq(fname, vertices, k):
    edges = dict()
    lines = open(fname, 'r').read().split('\n')
    lines.pop()
    longestpath = 0

    for j in range(0, len(lines)):
        path = lines[j].split()
        if len(path) > longestpath:
            longestpath = len(path)

        for i in range(1,len(path) - 1):
            if i == 1 and j == 0:
                source = path[i]
            if i == len(path) - 2 and j == 0:
                sink = path[i+1]
            seq1 = vertices[path[i]].seq
            seq2 = vertices[path[i+1]].seq
            combine = seq1[-(k):] + seq2[k-1]
            if combine in edges.keys():
                continue
            else:
                edges[combine] = Edge(path[i], path[i+1])
                vertices[path[i]].neighbors.append(path[i+1])
                vertices[path[i]].outdegree += 1
                vertices[path[i]].indegree += 1
    return edges, source, sink, longestpath

#reads kmers in....
def readsam(fname, k, edges, vertices):
    k = k+1
    lines = open(fname, 'r').read().split('\n')

    for z in range(0, len(lines)):
        line = lines[z]
        if len(line) == 0:
            continue
        elif line[0] == '@':
            continue
        elif len(line.split()) > 5:
            i = line.split()
            seq = i[9]
            
            j = 0
            while(j + k <= len(seq)):
                edge = seq[j:j+k]
                if edge in edges.keys():
                    edges.get(edge).weight += 1
                    u = edges.get(edge).v1
                    vertices.get(u).outcount += 1
                j += 1

def read_fasta(fname, k, edges, vertices):
    k = k+1
    lines = open(fname, 'r').read().split('\n')

    for z in range(1, len(lines), 2):
        seq = lines[z]
        j = 0
        while(j + k <= len(seq)):
            edge = seq[j:j+k]
            if edge in edges.keys():
                edges.get(edge).weight += 1
                u = edges.get(edge).v1
                vertices.get(u).outcount += 1
            j += 1

def edgewt_bfs(edges, vertices, src, k):
    queue = [src]
    inwt = {key: 0 for key in vertices}
    visited = {key: False for key in vertices}
    inwt[src] = 100
    visited[src] = True
    counts = dict()

    while queue:
        u = queue.pop(0)
        for v in vertices.get(u).neighbors:
            e = vertices[u].seq[-(k):] + vertices[v].seq[k-1]
            counts[e] = Edge(u, v)
            if vertices.get(u).outcount == 0:
                counts[e].weight = 0
                inwt[v] += inwt[u] / len(vertices.get(u).neighbors)
            else:
                counts[e].weight = (edges.get(e).weight / vertices.get(u).outcount) * inwt.get(u)
                inwt[v] += counts[e].weight

            if not visited[v]:
                queue.append(v)
                visited[v] = True
    return counts

def print_counts(fname, edges, vertices, source, sink):
    outputfile = open(fname, 'w')
    outputfile.write(str(int(len(vertices.keys())/2)) + '\n' + str(source) + '\n' + str(sink) + '\n')
    for e in edges.keys():
        outputfile.write(str(edges.get(e).v1) + ' ' + str(edges.get(e).v2) + ' ' + str(edges.get(e).weight) + ' ' +  '\n')
    outputfile.close()


def output_edges(outputfilename, edges, vertices):
    """ takes in outputfilename, edge map, kwrites normalized edge weights
        edges = map of k+1mer to edge with v1, v2, weight """
    outputfile = open(outputfilename, 'w')
    outputfile.write('#graph 1\n')
    outputfile.write(str(int(len(vertices.keys())/2)) + '\n')
    for e in edges.keys():
        w = round(edges.get(e).weight)
        if w == 0:
            w = -1
        outputfile.write(str(edges.get(e).v1) + ' ' + str(edges.get(e).v2) + ' ' + str(w) + ' ' +  '\n')
    outputfile.close()

def output_inexact_edges(outputfilename, edges, vertices):
    """ takes in outputfilename, edge map, kwrites normalized edge weights
        edges = map of k+1mer to edge with v1, v2, weight """
    outputfile = open(outputfilename, 'w')
    outputfile.write('#graph 1\n')
    outputfile.write(str(int(len(vertices.keys())/2)) + '\n')
    for e in edges.keys():
        w = round(edges.get(e).weight)
        if w == 0:
            lo = -1
            hi = -1
        else:
            lo = w - 3
            if lo < 0:
                lo = 0
            hi = w + 3
        outputfile.write(str(edges.get(e).v1) + ' ' + str(edges.get(e).v2) + ' ' + str(lo) + ' ' + str(hi) + ' ' + '\n')
    outputfile.close()