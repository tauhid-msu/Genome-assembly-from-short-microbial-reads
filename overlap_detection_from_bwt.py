import numpy as np
import argparse
import time
import math
import sys
from Bio import SeqIO
#import random
#import sa

def read_seqfile(input_file):
    sys.stderr.write("\tReading Sequence File {0} ... ".format(input_file))
  
    records = list(SeqIO.parse(input_file, "fasta"));
    nseq = len(records);
    sys.stderr.write(" Number of sequences: %d\n"%nseq)

    d_arr = np.empty(nseq, dtype='object')
    seq_labels = np.empty(nseq, dtype='object')
    i = 0;
    for record in SeqIO.parse(input_file, "fasta"):
        d_arr[i] = record.seq;
        seq_labels[i] = record.id
        i += 1;
    return (d_arr,seq_labels)

def updateBackward(Occ,C,c_index,l,u):
    
    l = C[c_index] + Occ[l-1,c_index];
    u = C[c_index] + Occ[u,c_index]-1;
    return (l,u)

def findDbseqid(l,isqprefix):
    nseq = seqarray.shape[0];
    for i in range(nseq):
        if isqprefix==1:
            if (rdbseq_offset[i]<=l and l<rdbseq_offset[i+1]):
                return i;
        elif isqprefix==0:
            if (dbseq_offset[i]<=l and l<dbseq_offset[i+1]):
               return i;  

def outputOverlap(i,l,j,isqprefix,quseqarray):
    #print(i,l,j) 
    db_i = findDbseqid(l,isqprefix);
    curOlap = [];
    curOlap.append(i);
    curOlap.append(db_i);
    olaplen = len(quseqarray[i])-j-1;
    curOlap.append(olaplen);
    curOlap.append(isqprefix);  
    allolaplists.append(curOlap);
    
def FMquery(L,Occ,C,ssa,T,th,isqprefix,quseqarray):
    
    dictionary = {"$":0,"A":1,"C":2,"G":3,"T":4}; 
    nqseq = quseqarray.shape[0];
    for i in range(nqseq):
        qlen = len(quseqarray[i])-1;
        j = qlen;
        curQseq = quseqarray[i];
        if isqprefix==1:     
            curQseq = ''.join(reversed(curQseq));    
            
        c = curQseq[j];
        c_index = dictionary[c];
        l = C[c_index];
        u = C[c_index+1] - 1;
               
        j -= 1;
        while (l<=u and j>=0):  
            if (qlen-j+1)>=th:
                l1,u1 = updateBackward(Occ,C,0,l,u);
                if l1<=u1:
                    for rr in range(l,u+1,1):
                        if ssa[rr]>0:
                            if T[ssa[rr]-1]=='$':
                                outputOverlap(i,ssa[rr],j,isqprefix,quseqarray);
                        else:
                            outputOverlap(i,0,j,isqprefix,quseqarray);
                    
            c = curQseq[j]; 
            c_index = dictionary[c];
            l,u = updateBackward(Occ,C,c_index,l,u);
            j -= 1;


def seq_db_text(textOrder):
    text = "";
    seq_offset = [];
    for i in range(seqarray.shape[0]):
        seq_offset.append(len(text));
        if textOrder==0:
            text +=seqarray[i];
        elif textOrder==1:
            rseq = ''.join(reversed(seqarray[i]));
            text +=rseq;
        if i<(seqarray.shape[0]-1):
            text += "$";
    seq_offset.append(len(text)-1);
    T = text.strip() 
    
    return T,seq_offset;

def createSuffixArray(T):
    myT = []
    for chr in T:
        myT.append(ord(chr))
    sa = ksa(myT)
    return sa;
   
def addNewQseqs(stPos):

    oldcount = len(allqlabels);
    #print('oldcount',oldcount);
    #print(stPos);
    #print(len(allolaplists));
    
    redQseqs = [];
    redQlabels = [];

    for i in range(stPos,len(allolaplists)):
        if len(allolaplists[i])>0:
            #print(allolaplists[i]);
            rqPos = allolaplists[i][1];
            if seqlabels[rqPos] not in allqlabels:
                allqlabels.append(seqlabels[rqPos]);
                allqseqs.append(seqarray[rqPos]);
            if seqlabels[rqPos] not in redQlabels:
                redQlabels.append(seqlabels[rqPos]);
                redQseqs.append(seqarray[rqPos]); 
    
    newcount = len(allqlabels);
    #print('newcount',newcount);

    redQseqarray = np.empty(len(redQlabels), dtype='object');
    #redQlabelsarr = np.empty(len(redQlabels), dtype='object');    
 
    for i in range(len(redQlabels)):
        redQseqarray[i] = redQseqs[i];
        #redQlabelsarr[i] = redQlabels[i];
    return redQseqarray;

def populate_bwt(bwtfile,lenT):
    F = [];
    L = [];
    rF = [];
    rL = [];

    i = 0;
    with open(bwtfile) as f:
        for line in f:
            line_content = line.strip("\n").split(" ");
            if i<lenT:
                F.append(int(line_content[0]));
                L.append(line_content[1]);
            else:
                rF.append(int(line_content[0]));
                rL.append(line_content[1]);
            i+=1;
   
    f.close();
   
    return (F,L,rF,rL);

def populate_FM_index(F,L):
    dictionary = {"$":0,"A":1,"C":2,"G":3,"T":4};
    salen = len(F);
    Occ = np.zeros(shape=(salen,5),dtype=int);
    C = [];

    for i in range(salen):
        if i>0:
            Occ[i,:] = Occ[i-1,:];
        ch = dictionary[L[i]];
        Occ[i,ch]+=1;

    for i in range(Occ.shape[1]+1):
        if i<=0:
            C.append(int(0));
        elif i>0:
            C.append(Occ[salen-1][i-1]+C[i-1]);

    return (Occ,C)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bwtfile', help='Input BWT index file')
    parser.add_argument('seqfile', help='Inupt filename')
    parser.add_argument('qseqfile', help='Query filename')
    parser.add_argument('nOverlap', help='Length of overlap', type=int)
    parser.add_argument('outpath', help='Output file to store overlaps')
    args = parser.parse_args()

    global seqarray
    global seqlabels
    seqarray,seqlabels = read_seqfile(args.seqfile);
        
    qseqarray,qseqlabels = read_seqfile(args.qseqfile);

    global allqseqs;
    global allqlabels;
    allqseqs = [];
    allqlabels = [];
    for i in range(qseqlabels.shape[0]):
        allqseqs.append(qseqarray[i]);
        allqlabels.append(qseqlabels[i]);

    global dbseq_offset;
    T,dbseq_offset = seq_db_text(0);
    #sa = createSuffixArray(T);  
    
    global allolaplists
    allolaplists = [[]]   
           
    global rdbseq_offset;
    rT,rdbseq_offset = seq_db_text(1);
    #rsa = createSuffixArray(rT);

    print('T',len(T));
    print('rT',len(rT));

    F,L,rF,rL = populate_bwt(args.bwtfile,(len(T)+1));
    Occ,C = populate_FM_index(F,L);    
    rOcc,rC = populate_FM_index(rF,rL);

    
    stPos = len(allolaplists);
    #F,L,Occ,C = buildFMindex(sa,T);
    FMquery(L,Occ,C,F,T,int(args.nOverlap),0,qseqarray);
    
    print('olapcount',stPos);  
    print('olapcount',len(allolaplists));    
    
    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs)); 
 
    stPos = len(allolaplists);
    FMquery(L,Occ,C,F,T,int(args.nOverlap),0,newQseqs);
    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs));
    print(len(allqlabels));
    print(len(allolaplists));
    
    #Now work with reverse reads    
    #rF,rL,rOcc,rC = buildFMindex(rsa,rT);
    stPos = len(allolaplists); 
    FMquery(rL,rOcc,rC,rF,rT,int(args.nOverlap),1,qseqarray);     
    
    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs));
  
    stPos = len(allolaplists);
    FMquery(rL,rOcc,rC,rF,rT,int(args.nOverlap),1,newQseqs);
    newQseqs = addNewQseqs(stPos);

    print('Final Olaplist length',len(allolaplists));
    print('Final read count',len(allqlabels))
    
    '''
    file_allolap = open(args.outpath+"/overlap_list","w");
    for i in range(1,len(allolaplists),1):
        if allolaplists[i][3]==0:
            file_allolap.write("%s %s %d\n"%(qseqlabels[allolaplists[i][0]],seqlabels[allolaplists[i][1]], allolaplists[i][2] ) );
        elif allolaplists[i][3]==1:
            file_allolap.write("%s %s %d\n"%(seqlabels[allolaplists[i][1]],qseqlabels[allolaplists[i][0]], allolaplists[i][2] ) );
    file_allolap.close();
    '''    
         
    seed_file = open(args.outpath+"/overlap_reads.fa","w")
    for i in range (len(allqlabels)):
        if len(allqlabels[i])>0:
            seed_file.write(">%s\n"%(allqlabels[i]));
            seed_file.write("%s\n"%(allqseqs[i]))
    seed_file.close(); 	
    
if __name__ == '__main__':
    main()

