class Read:
    """ Class Node to represent a vertex in the de bruijn graph """

    def __init__(self, lab, start, end, seq):
        self.label = lab
        self.start = start
        self.end = end
        self.seq = seq

class Coord:
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.weight = 1

class Kmer:
    def __init__(self, lab):
        self.label = lab
        self.weight = 0
        self.numcoords = 0
        self.coords = dict()
        self.originalmap = ''
        self.start = 0
        self.end = 0
    
def readvertices(fname ,k):
    lines = open(fname, 'r').read().split('\n')
    lines.pop()
    vertices = {line: Kmer(line) for line in lines}
    return vertices

def readsam(fname, k):
    k = k+1
    lines = open(fname, 'r').read().split('\n')

    reads = dict()
    kmers = [0]*30000
    repeatlist = set()

    for z in range(0, len(lines)):
        line = lines[z]
        if len(line) == 0:
            continue
        elif line[0] == '@':
            continue
        elif len(line.split()) > 5:
            i = line.split()
            name = i[0]
            if name in reads.keys():
                del reads[name]
                repeatlist.add(name)
            elif name not in repeatlist:
                start = int(i[3])
                seq = i[9]
                cig = i[5]
                
                lastindex = 0
                lastseq = 0
                letters = 0
                
                # parse cigar code, modify data as necessary
                for c in range(0, len(cig)):
                    if cig[c] == 'M':
                        #match... count it in the seq, move onto next
                        lastseq += int(cig[lastindex:c])
                        lastindex = c + 1
                        letters += 1
                    if cig[c] == 'H':
                        #hard cuts... not included in seq, ignore
                        lastindex = c + 1
                        letters += 1
                    if cig[c] == 'S':
                        slength = int(cig[lastindex:c]) 
                        seq = seq[0:lastindex] + seq[lastindex + slength: len(seq)]
                        lastindex = c + 1
                        letters += 1
                    if cig[c] == 'I':
                        #insertion from reference
                        ilength = int(cig[lastindex:c])
                        newseq = seq[lastseq:len(seq)]
                        seq = seq[0:lastseq + k - 1] 
                        newcig = cig[c+1: len(cig)]
                        newstart = start + lastseq - ilength
                        cig = cig[0:c + 1]

                        lastseq += 1
                        lastindex = c + 1
                        letters += 1

                        lines.append(str(z) + ' 1 2 ' + str(newstart) + 
                                    ' 4 ' + newcig + ' 6 7 8 ' + newseq)

                        break
                    if cig[c] == 'D':
                        #deletion from reference.
                        dlength = int(cig[lastindex:c])
                        newseq = seq[lastseq:len(seq)]
                        seq = seq[0:lastseq + k - 1] 
                        newcig = cig[c+1: len(cig)]
                        dlength = int(cig[lastindex:c])
                        newstart = start + lastseq + dlength
                        cig = cig[0:c + 1]

                        lastseq += 1
                        lastindex = c + 1
                        letters += 1

                        lines.append(str(z) + ' 1 2 ' + str(newstart) + 
                                    ' 4 ' + newcig + ' 6 7 8 ' + newseq)
                        break
                
                end = start + len(seq) - 1
                reads[name] = Read(name, start, end, seq)

                for i in range(start, end - k + 2):
                    kmers[i] += 1

    return kmers, reads

def iterate(reads, vertices, k):
    """ takes in vertex map, read map, kval, iterates thru cleaned reads, return updated vertex map """
    k = k+1
    for r in reads.keys():
        read = reads.get(r)
        i = 0
        while(i + k <= len(read.seq)):
            node = read.seq[i:i+k]
            nstart = read.start + i
            nend =read.start + i + k -1
            if node in vertices.keys():
                vertices.get(node).weight += 1
                if nstart in vertices.get(node).coords.keys():
                    vertices.get(node).coords[nstart].weight += 1
                else:
                    vertices.get(node).coords[nstart] = Coord(nstart, nend)
                    vertices.get(node).numcoords += 1
                    vertices.get(node).originalmap = read.label
            i += 1
    return vertices


def bestfit(vertices):
    for v in vertices.keys():
        maxwt = 0
        for c in vertices.get(v).coords.keys():
            if vertices.get(v).coords.get(c).weight > maxwt:
                mxwt = vertices.get(v).coords.get(c).weight
                vertices.get(v).start = vertices.get(v).coords.get(c).start
                vertices.get(v).end = vertices.get(v).coords.get(c).end

def normalize(vertices, kmers):
    """ takes in vertex map, kmer array, returns normalized vertex map """
    for v in vertices.keys():
        if kmers[(vertices.get(v).start)] == 0:
            continue
        else:
            w = vertices.get(v).weight / kmers[vertices.get(v).start]
            vertices.get(v).weight = w

def output_edges(outputfilename, vertices, k):
    """ takes in outputfilename, vertex map, kwrites normalized edge weights """
    outputfile = open(outputfilename, 'w')
    outputfile.write('# counted occurences in reads' + '\n')
    outputfile.write('however many nodes lol...' + '\n')
    for v in vertices.keys():
        if vertices.get(v).weight > 0:
            outputfile.write(str(v[:k]) + ' ' + str(v[-k:]) + ' ' + str(vertices.get(v).weight) + '\n')
    outputfile.close()