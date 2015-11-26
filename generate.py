#!/usr/bin/env python3

import sys
import random

if len(sys.argv)<5:
    print("usage: ./generate.py size_of_network size_of_community rate_of_inner_edge rate_of_outer_edge ")
    exit(1)

n=int(sys.argv[1]) # size of network
c=int(sys.argv[2]) # number of comunity
innerp=float(sys.argv[3]) # rate of inner edge
outerp=float(sys.argv[4]) # rate of outer edge
m=0

f=open("data/network.dat","w")
for u in range(n):
    for v in range(u+1,n):
        if (u<c and v<c)  and random.random()<innerp:
            f.write(str(u)+" "+str(v)+"\n")
            m+=1
        if (u>=c or v>=c) and random.random()<outerp:
            f.write(str(u)+" "+str(v)+"\n")
            m+=1
f.close()

f=open("data/group.dat","w")
for u in range(c):
    f.write(str(u)+"\n")
f.close()

print("n="+str(n)+", m="+str(m))

