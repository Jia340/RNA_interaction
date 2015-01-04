import sys,argparse
from data_structure import *
RE={}

def ParseArg():
    p=argparse.ArgumentParser(description="Remove repeat elements")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is the output file of Stitch-seq-Aligner.py")
    p.add_argument("-r","--repeat_masker",type=str,help="RepeatMasker file")
    p.add_argument('-o','--output',type=str,help="specify output file")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

def RepeatElement(chr,start,end):

  if "mmu-miR" in chr:
    return 1 
  if "ENSMUSG" in chr:
    chr=chr.split('_')[1]
  if chr=="chrMT" or chr=="chrNT":
    return 0
  S=0
  E=len(RE[chr])-1
  i=E/2
  res=RE[chr]
 
  while True:
#    print S,E,i
    
    if not (res[i][1]<start or end<res[i][0]): #overlap
      return 0

    if end<res[i][0]:
      E=i
      i=(i+S)/2
    elif res[i][1]<start:
      S=i
      i=(i+E)/2

    if E-S==1:
      if not (res[i+1][1]<start or end<res[i+1][0]): #overlap
        return 0
   
      return 1

def Main():
  args=ParseArg()
  file=open(args.repeat_masker,"r")
  frags_file=open(args.input,"r")
  output=open(args.output,"w")
  file.readline()
  global RE
  k=0
  print >>sys.stderr,"Input repeat masker"
  while True:
    line=file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if line[0] not in RE:
      RE[line[0]]=[]
    RE[line[0]].append((int(line[1]),int(line[2])))

  Types = ["snoRNA","protein_coding","lincRNA","tRNA","misc_RNA","pseudogene","miRNA","antisense","sense_intronic","non_coding","processed_transcript","sense_overlapping","protein_coding-noncoding"]
  print >>sys.stderr,"Input paired fragments"  
  while True:
    line=frags_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    
#    p1=annotated_bed(line[0:5]+line[6:9],proper=line[9],source=line[5],repeat=0)
#    p2=annotated_bed(line[11:16]+line[17:20],proper=line[20],source=line[16],repeat=0)
    '''
    if isinstance(p1.start, list):
      p1.start=int(p1.start[0])
      p1.end=int(p1.end[-1])
    if isinstance(p2.start, list):
      p2.start=int(p2.start[0])
      p2.end=int(p2.end[-1])
    '''
    if "dead" in line[0] or "dead" in line[11]:
      continue

    if line[6] not in Types or line[17] not in Types:
      continue
    if "," not in line[1]:
      if RepeatElement(line[0],int(line[1]),int(line[2]))==0:
        continue
    else:
      for i in range(len(line[1].split(","))):
        if RepeatElement(line[0],int(line[1][i]),int(line[2][i]))==0:continue

    if "," not in line[12]:
      if RepeatElement(line[11],int(line[12]),int(line[13]))==0:
        continue
    else:
      for i in range(len(line[12].split(","))):
        if RepeatElement(line[11],int(line[12][i]),int(line[13][i]))==0:continue

    print >>output,'\t'.join(line)
    k+=1
    if k%10000==0:
      print >>sys.stderr,"Fragments:%d\r"%(k),

if __name__=="__main__":
  Main()
    

  
  
