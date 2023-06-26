class Graph:
    def __init__(self):
        self.graph = dict()
        self.acyclic = True
        self.cc = dict()
        self.vertices = 0
        self.edges = 0
        self.sources = 0
        self.sinks = 0

    # edges are modeled as dictionary, v1 -> [(v2, weight)]
    def addEdge(self, node1, node2, cost):
        if node1 not in self.graph:
            self.graph[node1] = []
        if node2 not in self.graph:
            self.graph[node2] = []

        if node1 == 'xxx':
            self.sources += 1
        if node2 == 'yyy':
            self.sinks += 1
        if node1 != 'xxx' and node2 != 'yyy':
            self.edges += 1
        self.graph[node1].append((node2, int(cost)))

def read_input(fname):
    lines = open(fname, 'r').read().split('\n')

    return (lines[1], lines[2:len(lines) - 1])

def construct_graph(info):
    g = Graph()
    u = Graph()
    u.vertices = info[0]
    g.vertices = info[0]
    edges = info[1]

    for edge in edges:
        e = edge.split()
        g.addEdge(e[0], e[1], e[2])
        u.addEdge(e[0], e[1], e[2])
        u.addEdge(e[1], e[0], e[2])

    return (g, u)

def DFS_visit(g, v, visited, c):
    # takes in graph g, v to visit, visited list
    # mark vertex as visited
    visited[v] = 'gray'
    # get adjacency list
    adj = g.graph[v]

    # for each object in the adjacency list
    for a in adj:
        # the first of the tuple is the adj node
        v1 = a[0]
        
        if visited[v1] == 'gray':
            g.acyclic = False
        elif visited[v1] == 'white':
            # add to cc
            g.cc[c].append(v1)
            # visit the node
            DFS_visit(g, v1, visited, c)


    visited[v] = 'black'

def DFS(g):
    # create a visited list, where each v is mappped to False
    visited = dict()
    for v in g.graph.keys():
        visited[v] = 'white'
    c = 0

    # then for each vertex in the vertex list
    for v in g.graph.keys():
        # if not already, visit, create new cc
        if visited[v] == 'white':
            c += 1
            g.cc[c] = [v]
            DFS_visit(g, v, visited, c)

    return g

def print_reset(outputfilename):
    outputfile = open(outputfilename, 'w')
    outputfile.close()

def print_info(outputfilename, g, u):
    outputfile = open(outputfilename, 'a')
    outputfile.write('# vertices: ' + str(g.vertices) + '\n')
    outputfile.write('# non supersource/sink edges: ' + str(g.edges) + '\n')
    outputfile.write('# sources: ' + str(g.sources) + '\n')
    outputfile.write('# sinks: ' + str(g.sinks) + '\n')
    outputfile.write('graph is acyclic? ' + str(g.acyclic) + '\n')
    outputfile.write('# cc: ' + str(len(u.cc.keys())) + '\n')
    outputfile.close()

def print_cc_dummycheck(outputfilename, u):
    outputfile = open(outputfilename, 'a')
    v = 0
    for c in u.cc.keys():
        if (len(u.cc[c]) > 400):
            outputfile.write('CC #' + str(c) + ' size: ' + str(len(u.cc[c])) + '\n')
            for x in u.cc[c]:
                outputfile.write(x + '\n')
        v += len(u.cc[c])
    outputfile.write('total vertices counted: ' + str(v) + '\n')
    outputfile.write('total original vertices: ' + str(u.vertices))
    outputfile.close()
