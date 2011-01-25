#!/usr/bin/python

import fileinput

q = []
s = 5

for line in fileinput.input():
    if line:
        q.insert(0,float(line))
        if len(q) == s:
            #print str(q)+' -> '+str(sum(q)/float(s))
            print sum(q)/float(s)
            q.pop()
