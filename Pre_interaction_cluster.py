import argparse,sys
from data_structure import *
from copy import deepcopy
Ensembl_strand={}
miRNA_anno={}
Rev_Strand={"+":"-","-":"+"}

class mature_miRNA():
  def __init__(self,x,**kwargs):
    self.chr="chr"+x[0]
    self.start=int(x[3])
    self.end=int(x[4])
    self.strand=x[6]
    self.name=x[8].strip().split(";")[2][5:]
    if self.name.split("-")[-1]=="5p" or self.name.split("-")[-1]=="3p":
      self.txname="Mir"+"-".join(self.name.split("-")[2:-1])
    else:
      self.txname="Mir"+"-".join(self.name.split("-")[2:])
    for key in kwargs.keys():
      setattr(self,key,kwargs[key])

  def Genomic_position(self,s,e):
    if self.strand=="+":
      return (self.start+s-1,self.end+e-1)
    elif self.strand=="-":
      return (self.end-e+1,self.end-s+1)      

def ParseArg():
  p=argparse.ArgumentParser(description="Preprocessing for Select_strongInteraction_pp")
  p.add_argument("-i","--input",type=str,help="Paired fragment file")
  p.add_argument("-A","--annotation",type=str,default="/data2/jialu/RNA_RNA/new_mapping/annotation/mm9_Ensembl_Exon.txt",help="Ensembl annotation file for transcript")
  p.add_argument("-M","--miRNA",type=str,default="/data2/jialu/RNA_RNA/new_mapping/annotation/mmu_mm9.gff3",help="Genomic position of miRNA(miRbase)")
  p.add_argument("-o","--output",type=str,help="Output file name")
  if len(sys.argv)==1:
    print >>sys.stderr,p.print_help()
    sys.exit(0)
  return p.parse_args()

def Output_annotatedBed(bed):
  return '\t'.join([bed.chr,str(bed.start),str(bed.end),bed.strand,bed.seq,bed.source,bed.type,bed.name,bed.subtype,bed.proper])

def Convert_miRNA(bed):
  name=bed.name
  if name not in miRNA_anno:
    print name
    return 0
  new_pos=miRNA_anno[name].Genomic_position(bed.start,bed.end)
  bed.start=min(new_pos)
  bed.end=max(new_pos)
  bed.chr=miRNA_anno[name].chr
  if bed.proper=="ProperStrand":
    bed.strand=miRNA_anno[name].strand 
  else:
    bed.strand=Rev_Strand[miRNA_anno[name].strand]
  bed.name=miRNA_anno[name].txname
  return bed  

def Convert(half):
  if half.source=="miRNA":
    half=Convert_miRNA(half)
    return [half]
  elif half.source=="mm9_rna.fa":
    halfs=[]
    txID=half.chr.split("_")[0]
#    print txID
    half.chr=half.chr.split("_")[1]
#    print half.chr.split("_")
    if half.proper=="ProperStrand":
      half.strand=Ensembl_strand[txID]
    else:
      half.strand=Rev_Strand[Ensembl_strand[txID]]
    if isinstance(half.start, list):
      for i in range(len(half.start)):
        h=deepcopy(half)
        h.start=half.start[i]
        h.end=half.end[i]
        halfs.append(h)
      return halfs
    else:
      return [half]
  elif half.source=="genome":
    if half.type=="snoRNA" and "ENSMUS" in half.name:
      half.name=half.name.split("_")[1]
    
    if half.part==2:
      half.strand=Rev_Strand[half.strand] 
    return [half]


def Main():
  args=ParseArg()
  output=open(args.output,"w")    
  global Ensembl_strand
  global miRNA_anno

  Strand={"1":"+","-1":"-"}
  txanno_file=open(args.annotation,"r")
  txanno_file.readline()
  while True:
    line=txanno_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    if line[1] not in Ensembl_strand:
      Ensembl_strand[line[1]]=Strand[line[3]]
 
  miRNAanno_file=open(args.miRNA,"r")
  while True:
    line=miRNAanno_file.readline()
    if line=="":break
    line=line.strip().split("\t")
    p=mature_miRNA(line)
    if p.name not in miRNA_anno:
      miRNA_anno[p.name]=p

  Input=open(args.input,"r")
  while True:
    line=Input.readline()
    if line=="":
      break
    line=line.strip().split("\t")
    linker=line[10]
    p1=annotated_bed(line[0:10],part=1)
    p2=annotated_bed(line[11:],part=2)
    p1_converted=Convert(p1)
    p2_converted=Convert(p2)
    for i1 in p1_converted:
      for i2 in p2_converted:
        if i1!=0 and i2!=0:
          print >>output, Output_annotatedBed(i1)+"\t"+linker+"\t"+Output_annotatedBed(i2)

 
if __name__=="__main__":
  Main()
