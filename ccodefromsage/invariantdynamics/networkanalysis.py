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

graph.vertex_properties["seqpos"] = graph.new_vertex_property("int")
graph.vertex_properties["matpos"] = graph.new_vertex_property("vector<int>")
graph.vertex_properties["summable"] = graph.new_vertex_property("int")
graph.vertex_properties["selfloop"] = graph.new_vertex_property("int")
print "3"
for i in range(0,graph.num_vertices()):
    graph.vp["summable"][graph.vertex(i)] = 1
    graph.vp["seqpos"][graph.vertex(i)] = keys[i]
    graph.vp["matpos"][graph.vertex(i)] = list((keys[i]/grid, keys[i]%grid))

print "4"

access_dict={}
for i in range(0,len(keys)):
        access_dict[keys[i]]=i
for dom in keys:
    for im in nodes[dom]:
        graph.add_edge(access_dict[dom],access_dict[im])

for v in graph.vertices():
    if graph.edge(v,v) is not None:
        graph.vp["selfloop"][v] = 1
        print "found self loop!"

print "5"
if graph.num_vertices()<150:
    posg = sfdp_layout(graph)
    graph_draw(graph, posg,output ="outputs/graph.ps")
    interactive_window(graph,posg, geometry=(1000,1000),display_props=[graph.vp["matpos"],graph.vp["selfloop"]],display_props_size=16,vertex_size=15,edge_pen_width=2)
print "6"


for vert in graph.vertices():
    for neighbor in [out_edge.target() for out_edge in vert.out_edges()]:
        if (graph.edge(neighbor, vert) is None and
            abs(graph.vp["matpos"][neighbor][0]-graph.vp["matpos"][vert][0]) <= 1 and
                abs(graph.vp["matpos"][neighbor][1]-graph.vp["matpos"][vert][1]) <= 1):
                    graph.add_edge(neighbor, vert)


print "7"
components = label_components(graph)[0]

cg, community, vcount, ecount, va, ea = condensation_graph(graph, components, avprops = [graph.vp["summable"], graph.vp["matpos"],graph.vp["selfloop"]])




cg = transitive_closure(cg)

cg.vertex_properties["num"] = va[0]
cg.vertex_properties["numloops"] = va[2]
cg.vertex_properties["avpos"] = cg.new_vertex_property("string")
cg.vertex_properties["rec"]=cg.new_vertex_property("bool")

for cv in cg.vertices():
    cg.vp["avpos"][cv] = str([float(coord)/float(cg.vp["num"][cv]*grid) for coord in va[1][cv]])
    if cg.vp["numloops"][cv] + cg.vp["num"][cv] > 1:
        cg.vp["rec"][cv] = 1
    else:
        cg.vp["rec"][cv] = 0

print(cg.vp["num"])
print(cg.vp["numloops"])
print(cg.vp["rec"])

print "8"
print(cg.vp["rec"])
cg.set_vertex_filter(cg.vp["rec"])
cg.purge_vertices()
print "9"
cpos = sfdp_layout(cg)
graph_draw(cg,pos=cpos,output="outputs/condensation_graph.ps",geometry=[1000,1000])


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
interactive_window(cg,pos=cpos, geometry=(1000,1000),key_press_callback=callback,display_props = [cg.vp["avpos"],cg.vp["num"], cg.vp["numloops"],cg.vp["rec"]],display_props_size=14,vertex_size=15,edge_pen_width=2)
infile.close()
end= time.time()
print(end-start)
