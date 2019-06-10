import meta as m
import time
import copy
file = "sparse_10_30_3_8_I.full"

g = m.readGraph(file)
cg = copy.deepcopy(g)



start = time.time()
inf(cg)
m.process(cg)
writeSolution(cg,"borneInf_"+file,time.time()-start,"Borne inf")

start = time.time()
m.sup(g)
m.check(g)
writeSolution(cg,"borneSup_"+file,time.time()-start,"Borne inf")