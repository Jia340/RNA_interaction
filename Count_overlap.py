##input files are paired fragment

import sys

input1=open(sys.argv[1],"r")
input2=open(sys.argv[2],"r")
output=open("Intersections.txt","w")
sf=open("Single.txt","w")
set1=set()
set2=set()

def SingleFrag(line):
  return(line[0]==line[11] and line[7]==line[18] and len(line[4])==90 and len(line[15])==100 and line[9]==line[20]) 

while True:
  line=input1.readline()
  if line=="":break
  line=line.strip().split("\t")
  if SingleFrag(line):
    print >>sf,'\t'.join(line)
    continue
  set1.add("--".join(sorted([":".join([line[0],line[7],line[9]]),":".join([line[11],line[18],line[20]])])))

while True:
  line=input2.readline()
  if line=="":break
  line=line.strip().split("\t")
  if SingleFrag(line):
    print >>sf,'\t'.join(line)
    continue
  set2.add("--".join(sorted([":".join([line[0],line[7],line[9]]),":".join([line[11],line[18],line[20]])])))

print "Unique interaction in library 1: %d"%(len(set1))
print "Unique interaction in library 2: %d"%(len(set2))
print "Intersection between two libraries: %d"%(len(set1.intersection(set2)))
for item in set1.intersection(set2):
  print >>output,item
  
  
