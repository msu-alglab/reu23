import argparse

parser = argparse.ArgumentParser(
    description="""
    get our edgeweights!!
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-p', '--paths', type=str, help='seq file', required=True)
requiredNamed.add_argument('-n', '--nodes', type=str, help='seg file', required=True)
requiredNamed.add_argument('-o', '--output', type=str, help='Output filename', required=True)




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

def readseg(fname):
    lines = open(fname, 'r').read().split('\n')
    lines.pop()
    vertices = dict()
    for line in lines:
        s = line.split()
        plus = s[0] + '+'
        minus = s[0] + '-'
        vertices[plus] = s[1]
        vertices[minus] = reversecomplement(s[1])
    return vertices

args = parser.parse_args()

pathlines = open(args.paths, 'r').read().split('\n')
pathlines.pop()
paths = dict()
for i in range (0, len(pathlines)):
    pathinfo = pathlines[i].split()
    key = pathinfo[0]
    nodelist = [pathinfo[1]]
    #nodelist = [pathinfo[1][2:len(pathinfo[1]) - 2]]
    for j in range(2, len(pathinfo)):
        #nodelist.append(pathinfo[j][1:len(pathinfo[j]) - 2])
        nodelist.append(pathinfo[j])
    paths[key] = nodelist
    
nodes = readseg(args.nodes)
outputfile = open(args.output, 'w')

finalpaths = dict()
for p in paths.keys():
    pathlist = paths.get(p)
    seq = nodes.get(pathlist[0])
    for i in range(1, len(pathlist)):
        currnode = nodes.get(pathlist[i])
        seq += currnode[26:len(currnode)]
    finalpaths[p] = seq
    outputfile.write('>' + p + '\n')
    outputfile.write(seq)
    outputfile.write('\n')
    #print('>output_10_' + p)
    #print(seq + '\n')

outputfile.close()