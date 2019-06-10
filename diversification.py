import meta as m
import time

file = "sparse_10_30_3_8_I.full"
g = m.readGraph(file)
start = time.time()
m.multiStart(g)
writeSolution(g,"Diversification_"+file,time.time()-start,"Diversification")