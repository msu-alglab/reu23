class Node:
    """ Class Node to represent a vertex in the de bruijn graph """
    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0

def read_reads(fname, k):
    lines = open(fname, 'r').read().split('\n')
    output = []
    for i in range(1, len(lines), 2):
        output.append[lines[i]]
    return lines

def read_map(fname, k):
    lines = open(fname, 'r').read().split('\n')
    LH = dict()
    RH = dict()
    for i in range(1, len(lines), 2):
    #for i in range(0, len(lines)):
        l = lines[i][:(k-1)]
        r = lines[i][-(k - 1):]

        if l not in LH.keys():
            LH[l] = [lines[i]]
        else: 
            LH[l].append(lines[i])

        if r not in RH.keys():
            RH[r] = [lines[i]]
        else:
            RH[r].append(lines[i])

    return (LH, RH)

def db_original(reads, k):
    vertices = dict()
    # vertices maps a vertex to its node object
    # type string -> node
    edgeweights = dict()
    # edgeweights maps an edge to a weight
    # type string -> int

    for read in reads:
        i = 0

        while i+k < len(read):
            v1 = read[i:i+k]
            v2 = read[i+1:i+k+1]
            e12 = v1 + ' ' + v2

            if v1 in vertices.keys():
                vertices[v1].outdegree += 1
                if e12 in edgeweights.keys():
                    edgeweights[e12] = edgeweights.get(e12) + 1
                else: 
                    edgeweights[e12] = 1                 
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edgeweights[e12] = 1
                            
            if v2 in vertices.keys():
                vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1

            i += 1

def map_constructor(x, k):
    prefix = x[0]
    suffix = x[1]

    edgewts = dict()
    vertices = dict()

    #iterate thru all the prefixes
    for p in prefix.keys():

        #for each prefix, create corresponding vertices
        for v1 in prefix[p]:
            if v1 not in vertices.keys():
                vertices[v1] = Node(v1)

            #check if there is corresponding suffix
            if p in suffix.keys():
                #create vertices as necessary
                for v2 in suffix[p]:
                    if v2 not in vertices.keys():
                        vertices[v2] = Node(v2)

                    #add edges
                    e12 = v2 + ' ' + v1
                    #e12 = v1 + '>' + v2
                    edgewts[e12] = 1

                    #update in/outdegree
                    vertices[v2].outdegree += 1
                    vertices[v1].indegree += 1

    return (vertices, edgewts)

def basic_iterator(reads, k):
    """ Construct de bruijn graph from sets of short reads with k length word"""
    readlen = len(reads)
    vertices = dict()
    # vertices maps a vertex to its node object
    # type string -> node
    edgeweights = dict()
    # edgeweights maps an edge to a weight
    # type string -> int

    #iterate thru vertices
    for i in range(0, readlen - 2):
        v1 = reads[i]
        if v1 not in vertices.keys():
            vertices[v1] = Node(v1)

        #compare to each vertex occuring after
        for j in range(0, readlen - 1):
            v2 = reads[j]
            if v2 not in vertices.keys():
                vertices[v2] = Node(v2)
            
            #if there exists an overlap, create edge
            if v1[-(k - 1):] == v2[:(k-1)]:
                    e12 = v1 + ' ' + v2

                    vertices[v1].outdegree += 1
                    vertices[v2].indegree += 1
                    if e12 in edgeweights.keys():
                        edgeweights[e12] = edgeweights.get(e12) + 1
                    else: 
                        edgeweights[e12] = 1                 

    return (vertices, edgeweights)

def print_reads(outputfilename, reads):
    #for testing purposes
    outputfile = open(outputfilename, 'w')
    for read in reads:
        outputfile.write(read + '\n')
    outputfile.close

def source_sink(g):
    V = g[0]
    E = g[1]

    sources = []
    sinks = []

    for v in V.keys():
        if V[v].indegree == 0:
            sources.append(v)
        if V[v].outdegree == 0:
            sinks.append(v)

    if len(sources) > 1:
        V['xxx'] = Node('xxx')
        for s in sources:
            ess = 'xxx ' + s
            E[ess] = V[s].outdegree

    if len(sinks) > 1:
        V['yyy'] = Node('yyy')
        for i in sinks:
            ess = i + ' yyy'
            E[ess] = V[i].indegree

    return(V, E)

def print_reads(outputfilename, reads, k):
    #for testing purposes
    outputfile = open(outputfilename, 'w')
    for read in reads:
        outputfile.write(read[-(k - 1):] + '\n')
    outputfile.close()
    
def print_reset(outputfilename):
    outputfile = open(outputfilename, 'w')
    outputfile.close()

def print_runtime(outputfilename, time):
     outputfile = open(outputfilename, 'a')
     outputfile.write('# Runtime: ' + str(time) + '\n')
     outputfile.close()

def print_graph(outputfilename, g):
    """ Print the information in the graph to be (somewhat) presentable """
    outputfile = open(outputfilename, 'a')

    V = g[0]
    E = g[1]
    outputfile.write(str(len(V.keys())) + '\n')

    for a in E.keys():
        outputfile.write(a + ' ' + str(E[a]) + '\n')
        #outputfile.write(a + '\n')
    outputfile.close()