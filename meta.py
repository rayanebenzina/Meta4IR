import numpy as np
import math
import copy
import random
import time
import functools
X = [int(1000/2**i) for i in range(0,int(math.log2(1000))+1)]

class Graph:
    def __init__(self,idSafeNode):
        self.maxlvl, self.lvls, self.nodes = 0, [[]], {}
        self.tree = Node(self,idSafeNode)
        self.evacNodes, self.order = [],[]
    def end(self):
        return self.tree.end()

class Node:
    def __init__(self,g,id,parent = None, c = math.inf, l = 0):
        self.g, self.id, self.parent, self.c, self.l = g, id, parent, c, l
        self.successors, self.states, self.lvl, self.g.nodes[id], self.t = {}, [], 0,self,0
        if parent:
            self.lvl = parent.lvl + 1
            self.parent.successors[id] = self
            self.c = min(self.c,parent.c)
            self.t = parent.t + self.l
        if(self.lvl > g.maxlvl):
            g.maxlvl+=1
            g.lvls+=[[self]]
        else:
            g.lvls[self.lvl] += [self]

    def fusionL(self,lastStates = None, last = None):
        states = []
        for s in self.successors.values():
            sc = copy.deepcopy(lastStates) if s==last else copy.deepcopy(s.states)
            for state in sc:
                state[0]+=s.l
            states = fusion(states,sc)
        return states
    def fusion(self):
        self.states = Node.fusionL(self)
    def testC(self,s = None):
        return functools.reduce(lambda a,b: not a is False and a+b[1]<=self.c and b[1]+a, s or self.states, 0) is not False
    def end(self):
        return not len(self.states) == 0 and int(self.states[-1][0]) + self.t

class EvacNode(Node):
    def __init__(self, g, id, parent,c,l,population, max_rate):
        Node.__init__(self,g,id,parent,c,l)
        self.m = min(max_rate,self.c)
        self.p = population
        g.evacNodes+=[self]
        self.statescopy = []
    def fusion(self):
        Node.fusion(self)
        self.states = fusion(self.states,copy.deepcopy(self.statescopy))
    def fusionL(self,lastStates,last):
        return lastStates if not last else fusion(Node.fusionL(self,lastStates,last),copy.deepcopy(self.statescopy))
    def save(self):
        self.statescopy = copy.deepcopy(self.states)
    def getsave(self):
        return self.statescopy
    def set(self,a,b = None,c = None,d = None):
        self.states = [[a,b],[c,d] if c and d else [a+self.p/b,-b]] if b else a
        self.save()
    def sub(self,x):
        s = self.statescopy
        return s[0][0]-x >=0 and [[s[0][0]-x,s[0][1]],[s[1][0]-x,s[1][1]]]
    def subL(self,x):
        self.set(self.sub(x))

def readGraph(filename):
    f = open(filename, "r",encoding='cp437');
    f.readline()
    nb, id = [int(i) for i in f.readline().split()]
    g = Graph(id)
    paths = [[int(i) for i in f.readline().split()] for j in range(nb)]
    g.order = evacNodes = [line[0] for line in paths]
    f.readline()
    nbArcs=int(f.readline().split()[1])
    arcs = {(n1,n2):(l,c) for n1,n2,_,l,c in [[int(j) for j in f.readline().split()] for i in range(nbArcs)]}
    for idEvac,p,m,k,*path,idSafe in paths:
        addPath(g,arcs,evacNodes,g.tree,path[::-1]+[idEvac],p,m)
    g.lvls = g.lvls[::-1]
    return g

def readSolution(g,s):
    f = open(s,"r");f.readline()
    [g.nodes[int(id)].set(int(start),int(rate)) for i in range(int(f.readline())) for id,rate,start in [f.readline().split(" ")]]

def writeSolution(g,d,s,t,m):
    f = open(d+"solution_"+s,"w")
    f.write(s.split('.')[0] + "\n"+str(len(g.evacNodes))+ "\n")
    for st in ((g.nodes[id].getsave(),id) for id in g.order):
        f.write(str(st[1]) + " " + str(st[0][0][1]) + " " + str(st[0][0][0])+ "\n")
    f.write(("valid\n" + str(g.end()) if check(g) else "invalid\ninf") + "\n"+ str(t) + "\n" + m + "\n" + "Khalil et Rayane")

def addPath(g,a,e,n,path,p,m):
    if(path == []):
        return
    id = path[0]
    if not id in g.nodes.keys():
        x = (id,n.id) if id < n.id else (n.id,id)
        if(id in e):
            node = EvacNode(g,id,n,a[x][1],a[x][0],p,m)
        else:
            node = Node(g,id,n,a[x][1],a[x][0])
    else:
        node = g.nodes[id]
    return addPath(g,a,e,node,path[1:],p,m)

def inf(g):
    [n.set(0,n.c) for n in g.evacNodes]

def sup(g):
    [n.set(s and s.end()+1,n.c) for n,s in zip(g.evacNodes,[0]+g.evacNodes)]

def rand(g):
    rates = random.shuffle(g.evacNodes) or [random.randint(1,int(n.c)) for n in g.evacNodes]
    [n.set(s and s.end()+1,r) for n,s,r in zip(g.evacNodes,[0]+g.evacNodes,rates)]

def fusion(t1,t2):
    d1,d2 = {x:y for x,y in t1},{x:y for x,y in t2}
    return sorted([[i,(i in d1 and d1[i])+(i in d2 and d2[i])] for i in {*d1,*d2}])

def process(g):
    [n.fusion() for level in g.lvls for n in level]

def check(g):
    for level in g.lvls:
        for n in level:
            if n.fusion() or not n.testC():return False
    return True

def checkVoisin(g,n,s,l = None):
    return not n or not (not s or l and not l.testC(s)) and checkVoisin(g,n.parent,n.fusionL(s,l),n)

def descenteTemps(g):
    for n in sorted(g.evacNodes,key=lambda node:node.states[0][0]):
        for x in X:
            while(checkVoisin(g,n,n.sub(x))):
                n.subL(x)
        process(g)

def bestVoisins(g):
    minimum = g.tree.states[-1][0]
    best = None
    for n in g.evacNodes[:-1]:
        cg = copy.deepcopy(g)
        cg.nodes[n.id].subL(-round(minimum))
        process(cg)
        descenteTemps(cg)
        if(cg.tree.states[-1][0] < minimum):
            minimum = cg.tree.states[-1][0]
            best = cg
    return best

def multiStart(g,nb = 3):
    s = 40000
    ref = s
    cg = None
    for i in range(nb):
        while(s>=ref):
            cg = copy.deepcopy(g)
            rand(cg)
            process(cg)
            descenteTemps(cg)
            while 1:
                ccg = bestVoisins(cg)
                if(ccg==None):
                    break
                cg=ccg
            s=cg.tree.states[-1][0]
        print(int(s))
        ref = s
    return cg

#g = readGraph("sparse_10_30_3_8_I.full")
#start = time.time()
#g = sup(g)
#writeSolution(g,"sparse_10_30_3_8_I.full",time.time()-start,"Descente")