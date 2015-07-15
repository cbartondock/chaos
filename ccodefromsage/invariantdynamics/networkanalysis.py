#!/usr/local/bin/python
from graph_tool.all import *
import numpy as np
import matplotlib.cm as cm
import matplotlib.pylab as plt
import time
import cPickle as pickle

start = time.time()

infile = open("outputs/save.p","rb")
nodes = pickle.load(infile)
print "1"
keys = sorted(nodes.keys())
p = np.loadtxt("outputs/parameters.txt")
grid = p[0]
imatrix = np.loadtxt("outputs/imatrix.txt").astype('uint')
print "2"
graph = Graph()
graph.add_vertex(len(keys))
summable = graph.new_vertex_property("int")
sequence_position = graph.new_vertex_property("int")
matrix_position = graph.new_vertex_property("vector<int>")
selfloop = graph.new_vertex_property("int")
graph.vertex_properties["seqpos"] = sequence_position
graph.vertex_properties["matpos"] = matrix_position
graph.vertex_properties["summable"] = summable
graph.vertex_properties["selfloop"] = selfloop
print "4"
for i in range(0,graph.num_vertices()):
    graph.vp["summable"][graph.vertex(i)] = 1
    graph.vp["seqpos"][graph.vertex(i)] = keys[i]
    graph.vp["matpos"][graph.vertex(i)] = (keys[i]/grid, keys[i]%grid)

print "5"

access_dict={}
for i in range(0,len(keys)):
        access_dict[keys[i]]=i
for dom in keys:
    for im in nodes[dom]:
        graph.add_edge(access_dict[dom],access_dict[im])

for v in graph.vertices():
    if graph.edge(i,i) is None:
        graph.vp["selfloop"][graph.vertex(i)] = 0
    else:
        graph.vp["selfloop"][graph.vertex(i)] = 1
        print "found self loop!"

print "6"
#pos = random_layout(graph)
#graph_draw(graph, pos=pos,output ="outputs/graph.ps")
print "7"
components, histogram = label_components(graph)

cg, community, vcount, ecount, va, ea = condensation_graph(graph, components, avprops = [graph.vp["summable"], graph.vp["matpos"],graph.vp["selfloop"]])




numprop = va[0]
sumposprop = va[1]
avposprop = graph.new_vertex_property("string")
selfloopprop = va[2]
recurrenceprop = graph.new_vertex_property("bool")
for cv in cg.vertices():
    if selfloopprop[cv] + numprop[cv] > 1:
        recurrenceprop[cv] = True
    else:
        recurrenceprop[cv] = False

print "8"
cg = transitive_closure(cg)
cg.set_vertex_filter(recurrenceprop)
cg.purge_vertices()
vc = cg.vc_property_map("")
print "9"
cpos = sfdp_layout(cg)
graph_draw(cg,pos=cpos,vertex_fill_color=deg,output="outputs/condensation_graph.ps",geometry=[1000,1000])


for vert in graph.vertices():
    avposprop[vert] = str([float(coord)/float(numprop[vert]*grid) for coord in sumposprop[vert]])

def callback(widget, g, keyval, picked, pos, vprops, eprops):
    if picked!=None and keyval==104L:
        plt.clf()
        q = [[0 for i in range(0,int(p[0]))] for j in range(0,int(p[0]))]
        print("\n")
        for v in graph.vertices():
            if community[picked]==components[v]:
                (i, j) = tuple(graph.vp["matpos"][v])
                q[i][j]=-.3
        plt.imshow(imatrix+q,vmin=0,vmax=1,
        interpolation="nearest",
        cmap= cm.bwr,
        extent=[p[1],p[1]+p[0]*p[3],p[2],p[2]+p[0]*p[4]])
        plt.show()
    elif keyval==105:
        plt.clf()
        plt.imshow(imatrix,vmin=0,vmax=1,interpolation='nearest',
        cmap=cm.bwr,
        extent=[p[1],p[1]+p[0]*p[3],p[2],p[2]+p[0]*p[4]])
        plt.show()
print "10"
interactive_window(cg,pos=cpos, geometry=(1000,1000),key_press_callback=callback,display_props = [numprop, avposprop,selfloopprop,recurrenceprop],display_props_size=14,vertex_size=15,edge_pen_width=2,self_loops=True)
infile.close()
end= time.time()
print(end-start)
