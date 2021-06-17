#!/usr/bin/env python
import os
# from os import listdir
import csv

   
seq_dict = {}
seqf=open("seq_prefixes.txt","r")
for seq in seqf:
   print(seq.rstrip())
   seq_dict[seq.rstrip()]=[0,0,0]

os.chdir("/production/at820/exons")
exonfiles=os.listdir(".")

print("exons")
for file in exonfiles:
   print(file)
   f=open(file,"r")
   sampleid = f.readline().rstrip()[1:]
   while True:
      line = f.readline().rstrip()
      if not line: break
      elif line[0]==">":
         sampleid = line[1:]
      else:
         seq = line.rstrip()
         seq_dict[sampleid][0]+=len(seq)
   f.close()

os.chdir("/production/at820/introns")
intronfiles=os.listdir(".")

print("introns")
for file in intronfiles:
   print(file)
   f=open(file,"r")
   sampleid = f.readline().rstrip()[1:-5]
   while True:
      line = f.readline().rstrip()
      if not line: break
      elif line[0]==">":
         sampleid = line[1:-5]
      else:
         seq = line.rstrip()
         seq_dict[sampleid][1]+=len(seq)
   f.close()

os.chdir("/production/at820/supercontig")
superfiles=os.listdir(".")

print("supercontigs")
for file in superfiles:
   print(file)
   f=open(file,"r")
   sampleid = f.readline().rstrip()[1:-5]
   while True:
      line = f.readline().rstrip()
      if not line: break
      elif line[0]==">":
         sampleid = line[1:-5]
      else:
         seq = line.rstrip()
         seq_dict[sampleid][2]+=len(seq)
   f.close()



# write kbp counts for each sample to csv
with open("species_bp_output.csv", "w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(("Sample","Exon bp","Intron bp","Supercontig bp"))
    for key, value in seq_dict.items():
        writer.writerow((key, value[0], value[1], value[2]))


## test in interactive python
# f.close()
# seq_dict={"3-19-E":[0,0,0],"3-19-F":[0,0,0]}
# f=open("exons/4471.FNA","r")
# sampleid = f.readline().rstrip()[1:]
# while True:
#    line = f.readline().rstrip()
#    if not line: break
#    elif line[0]==">":
#       print(sampleid)
#       sampleid = line[1:]
#    else:
#       seq = line.rstrip()
#       seq_dict[sampleid][0]+=len(seq)

# with open("testoutput.csv", "w") as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(("Sample","Exon bp","Intron bp","Supercontig bp"))
#     for key, value in seq_dict.items():
#         writer.writerow((key, value[0], value[1], value[2]))

