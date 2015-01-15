import argparse,sys
Mapped_trans_count={}
trans_len={}

def ParseArg():
  p=argparse.ArgumentParser(description="Calculate fpkm")
  p.add_argument("-i","--input",type=str,required=True,help="paired fragments.")
  p.add_argument("--RNA1",type=str,required=True,help="Original RNA1 file")
  p.add_argument("--RNA2",type=str,required=True,help="Orignial RNA2 file")
  p.add_argument("-s","--lib_size",action="store_true",help="Whether normalized by library size.")
  p.add_argument("-A","--annotation",type=str,required=True,help="Transcriptome annotation file. MUST BE sorted by transcript ID")
  p.add_argument("-o","--output",type=str,required=True,help="Output file name")
  if len(sys.argv)==1:
    print >>sys.stderr,p.print_help()
    sys.exit(0)
  return p.parse_args()

def Single_Frags(line):
  return(line[0]==line[11] and line[7]==line[18] and len(line[4])==90 and len(line[15])==100 and line[9]==line[20])  

def FPKM(lib_size,output):
  for t in Mapped_trans_count:
    if "mmu-miR" in t:
      fpkm=1e9*Mapped_trans_count[t]/(lib_size*22)
    else:
      if t not in trans_len:
        print >>sys.stderr,"Error:unfound transcript:%s"%(t)
        exit(0)
      fpkm=1e9*Mapped_trans_count[t]/(lib_size*trans_len[t])
#      print t
#      print Mapped_trans_count[t],trans_len[t]
    print >>output,'\t'.join([t,str(fpkm)])

def Trans_length(name,trans):
  if name in trans_len:
    print >>sys.stderr,"Error:Transcripts are not sorted"
    print >>sys.stderr,name
    exit(0)
  trans_len[name]=0
  for item in trans:
    trans_len[name]+=item[1]-item[0]+1

def Main():
  args=ParseArg()
  global Mapped_trans_count
  global trans_len
  RNA1={}
  RNA2={}
  output=open(args.output,"w")  
  sf=open("SingleF.txt","w")  
  print >>sys.stderr,"Input original files"
  rna1_file=open(args.RNA1,"r")
  rna2_file=open(args.RNA2,"r")
  while True:
    line=rna1_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if line[5]=="genome" or line[5]=="miRNA" or line[9]=="NonProperStrand":continue
    if line[10] not in RNA1:
      RNA1[line[10]]=[]
    RNA1[line[10]].append(line[0])
  
  while True:
    line=rna2_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if line[5]=="genome" or line[5]=="miRNA" or line[9]=="NonProperStrand":continue
    if line[10] not in RNA2:
      RNA2[line[10]]=[]
    RNA2[line[10]].append(line[0])
    
  anno_file=open(args.annotation,"r")
  anno_file.readline()
  
  c=0
  prev=""
  print >>sys.stderr,"reading annotation file"
  while True:
    line=anno_file.readline()
    if line=="":
      Trans_length(prev,trans)
      break 
    line=line.strip().split("\t")    
    t_id=line[0]
    if t_id!=prev:
      if c!=0:
        Trans_length(prev,trans)
      trans=[(int(line[6]),int(line[7]))]
      prev=t_id
    else:
      trans.append((int(line[6]),int(line[7])))
    c+=1
    
 # print >>sys.stderr,trans_len
  
  input_file=open(args.input,"r")
  total_t=0
  while True:
    line=input_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if Single_Frags(line):
      print >>sf,"\t".join(line)
      continue

    linker=line[10]
    if line[5]=="mm9_rna.fa" and line[9]=="ProperStrand":
      t_mapped=len(RNA1[line[10]])
      for t in RNA1[line[10]]:
        if t not in Mapped_trans_count:
          Mapped_trans_count[t]=0.0
        Mapped_trans_count[t]+=1.0/(t_mapped)
      total_t+=1
    elif line[5]=="miRNA" and line[9]=="ProperStrand":
      if line[0] not in Mapped_trans_count:
        Mapped_trans_count[line[0]]=0.0
      Mapped_trans_count[line[0]]+=1.0
      total_t+=1

    if line[16]=="mm9_rna.fa" and line[20]=="ProperStrand":
      t_mapped=len(RNA2[line[10]])
      for t in RNA2[line[10]]:
        if t not in Mapped_trans_count:
          Mapped_trans_count[t]=0.0
        Mapped_trans_count[t]+=1.0/(t_mapped)
      total_t+=1
    elif line[16]=="miRNA" and line[20]=="ProperStrand":
      if line[11] not in Mapped_trans_count:
        Mapped_trans_count[line[11]]=0.0
      Mapped_trans_count[line[11]]+=1.0
      total_t+=1
      if total_t%5000==0:
        print >>sys.stderr,"Processed paired fragments: %d"%(total_t)
  print >>sys.stderr,"Library size is: %d"%(total_t)
  if args.lib_size:
    FPKM(total_t,output)
  elif not args.lib_size:
    FPKM(1,output)
  
if __name__=='__main__':
  Main()
