# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 09:23:03 2023

@author: malia
"""

import networkx as nx
import time


#prim , kruskal , borouvka with difffrent data structure data structure 
def get_line():
     with open("ban5000w-0.01-adjlist.txt") as file:
         for i in file:
             yield i

myset = set()
G = nx.Graph()
file_name= "ban5000w-0.01-adjlist.txt"

lines_required = 10000000
gen = get_line()
chunk = [i for i, j in zip(gen, range(lines_required))]
for i in chunk:
    tmp = i.split()
    a = int(tmp[0])
    b = int(tmp[1].replace('\n','').replace('q',''))
    c=float(tmp[2].replace('\n','').replace('q',''))
    G.add_edge(a, b, weight=c )
    lines_required=lines_required-1
    if lines_required==0:
        break
print ("done reading the file")
numberOfNodes=100
G=G.subgraph(range(1,numberOfNodes+1))
#print('number of node = %d'% numberOfNodes)

def getLargestNode():
    mx=-1
    for node in G.nodes:
        if node>mx:
            mx=node
    return mx

class GraphList:
    def __init__(self):
        self.graph=[]
        edgs=G.edges()
        for edge in edgs:
            u=edge[0]
            v=edge[1]
            w=G.get_edge_data(u, v)['weight']
            self.graph.append((u,v,w))

    def viewGraph(self,lst):
        ch=input('want to view tree??(y/n)')
        if ch!='y':
            return
        for edge in lst:
            print ("%d -- %d (w:%f)" % (edge[0],edge[1],edge[2]))

    def nodeHasEdge(self,n, e):
        return e[0] == n or e[1] == n

    def hasPath(self,lst,u,v):
        queue=[]
        ns=[item for item in lst]
        visited=[]
        def getNeibors(node):
            neibors=[]
            for item in ns:
                if self.nodeHasEdge(node,item):
                    if item[0]==node:
                        neibors.append(item[1])
                    else:
                        neibors.append(item[0])
            return neibors
        queue.append(u)
        visited.append(u)
        while len(queue)!=0:
            current=queue[0]
            del queue[0]
            if current==v:
                return True
            neibors=getNeibors(current)
            for nei in neibors:
                if nei not  in visited:
                    queue.append(nei)
                    visited.append(nei)
        return False

    def prem(self):
        tree=[]
        edges=[elem for elem in self.graph]
        possibleEdge=[]
        visitedNodes=[edges[0][0]]
        start_time=time.time()
        def nodeHasEdge(n,e):
            return e[0]==n or e[1]==n
        def appendEdgeForNode(node):
            i=0
            n=len(edges)
            while i<n:
                if nodeHasEdge(node,edges[i]):
                    possibleEdge.append(edges[i])
                    del edges[i]
                    i=i-1
                    n=n-1
                i=i+1

        def getNextEdge():

            nxtEdg=None
            bestWeight=1000
            n=len(possibleEdge)
            i=0
            while i<n:
                if possibleEdge[i][1] in visitedNodes and possibleEdge[i][0] in visitedNodes:
                    del possibleEdge[i]
                    i=i-1
                    n=n-1
                else:
                    if possibleEdge[i][2]<bestWeight:
                        bestWeight=possibleEdge[i][2]
                        nxtEdg=possibleEdge[i]
                i=i+1
            return nxtEdg

        appendEdgeForNode(visitedNodes[0])
        while True:
            nxtEdg=getNextEdge()

            if nxtEdg is None:
                break
            if nxtEdg[0] not in visitedNodes:
                visitedNodes.append(nxtEdg[0])
                appendEdgeForNode(nxtEdg[0])
            if nxtEdg[1] not in visitedNodes:
                visitedNodes.append(nxtEdg[1])
                appendEdgeForNode(nxtEdg[1])
            possibleEdge.remove(nxtEdg)
            tree.append(nxtEdg)
        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)


def kruskal(self):
        tree = [] # store MST
        edges = [elem for elem in self.graph]#STORE ALL EFGES
        start_time = time.time()
        def getNextEdge():# return the edge eith lowest weight
            nxtEdg=None
            bestWeight=1000 #initialize best weight
            n=len(edges)
            i=0
            while i<n:
                if self.hasPath(tree,edges[i][0],edges[i][1]):#check if there is a cycle 
                    del edges[i]
                    i=i-1
                    n=n-1
                else:
                    if edges[i][2]<bestWeight:
                        bestWeight=edges[i][2]
                        nxtEdg=edges[i]
                i=i+1
            return nxtEdg
        while True:
            nxtEdg=getNextEdge()

            if nxtEdg is None:
                break
            tree.append(nxtEdg)
        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)


def boruvka(self):
        nodes=[]# create empty list so store nodes
        tree = [] # store edges of MST
        edges = [elem for elem in self.graph] # store all edges in the graph
        start_time = time.time()
        for edg in edges:#Loop through each edge in the list edges.
            if edg[0] not  in nodes:#Check if the first node of the current edge is not in the list nodes.
                nodes.append(edg[0])#If the first node is not in nodes, append it to the list.
            elif edg[1] not in nodes:#Check if the second node of the current edge is not in the list nodes.
                nodes.append(edg[1])
        def getBestEdg(node):#Define a function getBestEdg that takes in a node as an argument and returns the best edge for that node.
            bestweight=1000
            bestEdge=None
            for edg in edges:
                if self.nodeHasEdge(node,edg):
                    if edg[0]==node:
                        v=edg[1]
                    else:
                        v=edg[0]
                    if self.hasPath(tree,node,v):#Check if there is already a path between node and v in the current tree. If there is, 
                                                 #continue to the next iteration of the loop.
                        continue;
                    if edg[2]<bestweight:
                        bestEdge=edg#bestEdge=edg - Assign the current edge to bestEdge and
                        bestweight=edg[2]# Update the value of bestweight to the weight of the current edge.
            return bestEdge
        for node in nodes:
            bestEdg=getBestEdg(node)
            if bestEdg is not None:
                tree.append(getBestEdg(node))
        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)



class GraphMatrix:
    def __init__(self):
        n = getLargestNode() + 1
        self.graph = [[] for i in range(n)]
        edgs = G.edges()
        for edge in edgs:
            u = edge[0]
            v = edge[1]
            w = G.get_edge_data(u, v)['weight']
            self.graph[u].append((v, w))
            self.graph[v].append((u,w))

    def viewGraph(self,mtrx):
        ch = input('want to view tree??(y/n)')
        if ch != 'y':
            return
        for nodeIndex in range(len(mtrx)):
            for edge in mtrx[nodeIndex]:
                print("%d -- %d (w:%f)" % (nodeIndex, edge[0], edge[1]))
                for edg in mtrx[edge[0]]:
                    if edg[0]==nodeIndex:
                        mtrx[edge[0]].remove(edg)
                        break

    def hasPath(self, mtx, u, v):
        queue = []
        visited = []

        def getNeibors(node):
            return [v for (v,_) in mtx[node]]

        queue.append(u)
        visited.append(u)
        while len(queue) != 0:
            current = queue[0]
            del queue[0]
            if current == v:
                return True
            neibors = getNeibors(current)
            for nei in neibors:
                if nei not in visited:
                    queue.append(nei)
                    visited.append(nei)
        return False

    def hasEdge(self,g,u,v):
        for edg in g[u]:
            if edg[0]==v:
                return True
        return False

    def prem(self):
        start_time = time.time()
        tree=[[] for i in range(len(self.graph))]#is a list of lists where each sublist represents a node in the graph, and will hold its edges in the minimum spanning tree.
        visitedNodes=[1]#is initialized with node 1, which is the starting node for the Prim's algorithm.
        def getNextEdge(): #getNextEdge is defined to find the next edge that should be added to the minimum spanning tree
            nxtEdg=None
            bestWeight=1000
            for node in visitedNodes:#The function loops over the nodes in visitedNodes to find the next edge.
                for (v,w) in self.graph[node]:
                    if (not self.hasEdge(tree,node,v)) and v not in visitedNodes:##function checks if the edge (node,v) has already been added to the tree
                        if w<bestWeight:
                            bestWeight=w
                            nxtEdg=(node,v,w)
            return nxtEdg;

        while True:
            nxtEdg=getNextEdge()#nxtEdg=getNextEdge(): The value of nxtEdg is assigned the result of the function getNextEdge.

            if nxtEdg is None:
                break
            if nxtEdg[0] not in visitedNodes:#If the starting node of the next edge is not in visitedNodes, it is added to the list.
                visitedNodes.append(nxtEdg[0])
            if nxtEdg[1] not in visitedNodes:
                visitedNodes.append(nxtEdg[1])
            tree[nxtEdg[0]].append((nxtEdg[1],nxtEdg[2]))
            tree[nxtEdg[1]].append((nxtEdg[0], nxtEdg[2]))

        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)


    def kruskal(self):
        start_time = time.time()
        tree = [[] for i in range(len(self.graph))]

        def getNextEdge():# this function returns the next edge to be added to the minimum spanning tree.
            nxtEdg=None
            bestWeight=1000
            for node in range(len(self.graph)):
                for (v,w) in self.graph[node]:
                    if not self.hasPath(tree,node,v):#f the edge (v, w) does not form a cycle in the minimum
                                                     #spanning tree (checked using the "hasPath" function), and the weight of the edge (w) is less than the current best weight, the edge is updated as the next edge to be added.
                        if w<bestWeight:
                            bestWeight=w
                            nxtEdg=(node,v,w)
            return nxtEdg
        while True:
            nxtEdg=getNextEdge()
            if nxtEdg is None:
                break
            tree[nxtEdg[0]].append((nxtEdg[1], nxtEdg[2]))
            tree[nxtEdg[1]].append((nxtEdg[0], nxtEdg[2]))
        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)


    def boruvka(self):
        start_time = time.time()
        tree = [[] for i in range(len(self.graph))]

        def getBestEdg(node):#he function defines an inner function "getBestEdg()",
                              # which takes a node as an input and returns the best edge starting from that node. 
            bestweight=1000
            bestEdge=None
            for (v,w) in self.graph[node]:
                if not self.hasPath(tree,node,v): #If the edge (v, w) does not form a cycle and its                                       
        # weight is less than the current best weight, the edge is updated as the best edge.
                    if w<bestweight:
                        bestEdge=(node,v,w)
                        bestweight=w
            return bestEdge
        for node in range(len(self.graph)):
            bestEdg=getBestEdg(node)#The main loop then loops over all nodes in the graph and calls the "getBestEdg()"
                                    # function for each node.
            if bestEdg is not None:
                tree[bestEdg[0]].append((bestEdg[1], bestEdg[2])) #is added to the minimum spanning tree 
                                                        #by adding the edge to both the nodes in the "tree" list.
                tree[bestEdg[1]].append((bestEdg[0], bestEdg[2]))
        print("--- %s seconds ---" % (time.time() - start_time))
        self.viewGraph(tree)
        
mat=GraphMatrix()
li=GraphList()
while True:
    inp=input('''Enter data structure\n(1) list\n(2) matrix\n(3)exit\n''')
    if inp=='1':
        inp2=input('''Enter algorithm\n(1) prem\n(2) kruskal\n(3) boruvka\n''')
        if inp2=='1':
            print('<<LIST:PRIM>>')
            li.prem()
        elif inp2=='2':
            print('<<LIST:KRUSKAL>>')
            li.kruskal()
        elif inp2=='3':
            print('<<LIST:BORUVKA>>')
            li.boruvka()
    elif inp=='2':
        inp2 = input('''Enter algorithm\n(1) prem\n(2) kruskal\n(3) boruvka''')
        if inp2 == '1':
            print('MATRIX:PREM')
            mat.prem()
        elif inp2 == '2':
            print('MATRIX:KRUSKAL')
            mat.kruskal()
        elif inp2 == '3':
            print('MATRIX:BORUVKA')
            mat.boruvka()
    elif inp=='3':
        break



        