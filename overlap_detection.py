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

class Triple(object):
    """Represent each sortable character in R with three integers"""
    #todo: input validation, errors
    def __init__(self, T, idx, length):
        t_i = lambda i: T[i] if i < length else 0
        self._triple = (t_i(idx), t_i(idx + 1), t_i(idx + 2))
        self._index = idx
        self._rank = None
        self._rpos = None

    @property
    def triple(self):
        """Character for R_k strings"""
        return self._triple

    @property
    def index(self):
        """Position of R_k character in source string"""
        return self._index

    @property
    def rpos(self):
        """Sorted order of R_k charcter"""
        return self._rpos

    @rpos.setter
    def rpos(self, pos):
        self._rpos = pos

    @property
    def rank(self):
        """Sorted order of R_k charcter"""
        return self._rank

    @rank.setter
    def rank(self, pos):
        self._rank = pos

    def __repr__(self):
        return "Triple({0}, {1}, {2})".format(self.triple, self.index, self.rank)


class NonsamplePair(object):
    #todo: property decorators for validation
    def __init__(self, T, idx, S_i_ranks):
        self.index = idx
        self.pair = None
        max_index = len(T)
        if idx < max_index:
            self.pair = (T[self.index], S_i_ranks[self.index + 1])
        else:
            self.pair = (0, S_i_ranks[self.index + 1]) #defined to be 0 by KS algorithm 


# Recursive Karkkainen-Sanders implementation
#   Input: list of integers (representing characters)
#   Returns suffix array for list
def ksa(T):
    length = len(T) # n
    # B_k = { i \in [0,n] | i mod 3 = k }
    B_0, B_1, B_2 = xrange(0, length+1, 3), xrange(1, length+1, 3), xrange(2, length+1, 3)

    #karkkainen-sanders step 1: sort sample suffixes
    R_0 = [ Triple(T, idx, length) for idx in B_0 ]
    R_1 = [ Triple(T, idx, length) for idx in B_1 ]
    R_2 = [ Triple(T, idx, length) for idx in B_2 ]

    R = R_1 + R_2
    #enable reverse-lookup of characters in R from a list of sorted characters from R
    for i, r_char in enumerate(R):
        r_char.rpos = i
    sorted_suffixes_R = sorted(R, key=lambda suffix_char: suffix_char.triple)

    #Enables 0 as unique terminating character by starting ranks at 1
    def rank_suffixes(suffixes, rank=1):
        for i, suffix in enumerate(suffixes):
            if i > 0 and suffix.triple != suffixes[i-1].triple:
                rank += 1
            suffix.rank = rank
        return rank
    rank = rank_suffixes(sorted_suffixes_R)
    R_prime = [suffix.rank for suffix in R]

    #recursive call
    if (rank < len(R)):  #we had repeats of characters of R, make a recursive call to sort
        R_prime_suffix_array = ksa(R_prime)
    else:
        #directly form suffix array
        R_prime_suffix_array = [len(R)] + [suffix.rpos for suffix in sorted_suffixes_R]
    rank_Si = [None] * (length + 3) #why plus 3? -> additiionally define rank(S_(n+1) = rank(S_(n+2)) = 0
    rank_Si[-2] = rank_Si[-1] = 0

    #build rank(S_i) lookup array
    for i, SAi in enumerate(R_prime_suffix_array):
        if SAi < len(R): #ignore the index pointing to the terminating character of R_prime
            rank_Si[R[SAi].index] = i

    sorted_suffixes_R = [R[i] for i in R_prime_suffix_array[1:]]

    #karkkainen-sanders step 2: sort nonsample suffixes
    nonsample_suffix_pairs = [NonsamplePair(T, idx, rank_Si) for idx in B_0]
    sorted_nonsample_suffix_pairs = sorted(nonsample_suffix_pairs, key=lambda p: p.pair)

    #karkkainen-sanders step 3: merge
    cur_Sc, cur_Sb0 = 0, 0
    objs_SA = []
    def getT(idx):
        if idx < len(T):
            return T[idx]
        return 0
    while cur_Sc < len(sorted_suffixes_R) and cur_Sb0 < len(sorted_nonsample_suffix_pairs):
        i = sorted_suffixes_R[cur_Sc].index
        j = sorted_nonsample_suffix_pairs[cur_Sb0].index
        if i % 3 == 1: #i in B_1
            #S_i =< S_j iff (T[i], rank(S_t+1) =< (t_j, rank(s_j+1))
            if (getT(i), rank_Si[i+1]) < (getT(j), rank_Si[j+1]):
                objs_SA.append(sorted_suffixes_R[cur_Sc])
                cur_Sc += 1
            else:
                objs_SA.append(sorted_nonsample_suffix_pairs[cur_Sb0])
                cur_Sb0 += 1
        else: #i in B_2
            if (getT(i), getT(i+1),  rank_Si[i+2]) < (getT(j), getT(j+1), rank_Si[j+2]):
                objs_SA.append(sorted_suffixes_R[cur_Sc])
                cur_Sc += 1
            else:
                objs_SA.append(sorted_nonsample_suffix_pairs[cur_Sb0])
                cur_Sb0 += 1

    objs_SA += sorted_suffixes_R[cur_Sc:]
    objs_SA += sorted_nonsample_suffix_pairs[cur_Sb0:]
    SA = [suffix_object.index for suffix_object in objs_SA]
    return SA


def whichElt(Lelt):

    if Lelt=='$':
        return 0;
    elif Lelt=='A':
        return 1;
    elif Lelt=='C':
        return 2;
    elif Lelt=='G':
        return 3;
    elif Lelt=='T':
        return 4;

def buildFMindex(sa,T):

    Flist = [];
    Llist = [];
    salen = len(sa);
    Occ = np.zeros(shape=(salen,5),dtype=int); 
    C = [];
    
    for i in range(salen):
        if sa[i]==0:
            Llist.append('$');
        else:
            Llist.append(T[sa[i]-1]); 

        if i==0:
            Flist.append('$');
        else: 
            Flist.append(T[sa[i]]);
    
    for i in range(salen):
        if i>0:
            Occ[i,:] = Occ[i-1,:];
        Occ[i,whichElt(Llist[i])] += 1;
    
    for i in range(Occ.shape[1]+1):
        if i<=0:
            C.append(int(0));
        elif i>0:
            C.append(Occ[salen-1][i-1]+C[i-1]);    
    
    #for i in range(1,Occ.shape[1]+1,1):
    #    C[i] -= 1;        
    return (Flist,Llist,Occ,C)

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
    
def qseqOverlap(l,olaplen,isqprefix,qseqOlaps):
    db_i = findDbseqid(l,isqprefix);
    curOlap = [];
    curOlap.append(-1);
    curOlap.append(db_i);
    curOlap.append(olaplen);
    curOlap.append(isqprefix);
    qseqOlaps.append(curOlap);

def FMqueryForEachSeq(Occ,C,ssa,T,curQseq,th,isqprefix):
    qseqOlaps = [[]];
    dictionary = {"$":0,"A":1,"C":2,"G":3,"T":4};
    qlen = len(curQseq)-1;
    j = qlen;
    #curQseq = qseqarray[i];
    if isqprefix==1:
        curQseq = ''.join(reversed(curQseq));
    c = curQseq[j];
    c_index = dictionary[c];
    l = C[c_index];
    u = C[c_index+1] - 1;
    j -= 1; 
    
    while (l<=u and j>=0): 
        if (qlen-j+1)>=th:
            
            #print(th);
            l1,u1 = updateBackward(Occ,C,0,l,u);
            if l1<=u1:
                #print('lu$',l1,u1,'lu',l,u);
                for rr in range(l,u+1,1):
                    if ssa[rr]>0:
                        #print(rr,ssa[rr],T[ssa[rr]-1]);
                        olaplen = len(curQseq)-j-1;
                        if T[ssa[rr]-1]=='$':
                            qseqOverlap(ssa[rr],olaplen,isqprefix,qseqOlaps);
                        else:
                            #print(rr,ssa[rr],'$');
                            qseqOverlap(0,olaplen,isqprefix,qseqOlaps);
            c = curQseq[j];
            c_index = dictionary[c];
            l,u = updateBackward(Occ,C,c_index,l,u);
            #print(c,c_index,'lu',l,'->',u);
            
        j -= 1;
    
    return qseqOlaps;       
    
def FMquery(F,L,Occ,C,ssa,T,th,isqprefix,quseqarray):
    
    dictionary = {"$":0,"A":1,"C":2,"G":3,"T":4}; 
    nqseq = quseqarray.shape[0];
    for i in range(nqseq):
        qlen = len(quseqarray[i])-1;
        j = qlen;
        curQseq = quseqarray[i];
        if isqprefix==1:     
            #print('>',curQseq);
            curQseq = ''.join(reversed(curQseq));    
            #print('<',curQseq);

        c = curQseq[j];
        #print(c);
        #print(dictionary[c]);
        c_index = dictionary[c];
        l = C[c_index];
        u = C[c_index+1] - 1;
        #print('lu',l,'->',u);       
        j -= 1;
        while (l<=u and j>=0):  
            if (qlen-j+1)>=th:
                #print(th);
                l1,u1 = updateBackward(Occ,C,0,l,u);
                if l1<=u1:
                    #print('lu$',l1,u1,'lu',l,u);
                    for rr in range(l,u+1,1):
                        if ssa[rr]>0:
                            #print(rr,ssa[rr],T[ssa[rr]-1]);
                            if T[ssa[rr]-1]=='$':
                                outputOverlap(i,ssa[rr],j,isqprefix,quseqarray);
                        else:
                            #print(rr,ssa[rr],'$');
                            outputOverlap(i,0,j,isqprefix,quseqarray);
                    #outputOverlap(i,l,j,isqprefix);
                    #break;

            c = curQseq[j]; 
            c_index = dictionary[c];
            l,u = updateBackward(Occ,C,c_index,l,u);
            #print(c,c_index,'lu',l,'->',u);
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile', help='Inupt filename')
    parser.add_argument('qseqfile', help='Query filename')
    #parser.add_argument('nSeq', help='Number of sequences', default=100, type=int)
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
    sa = createSuffixArray(T);  
    
    global allolaplists
    allolaplists = [[]]   
           
    global rdbseq_offset;
    rT,rdbseq_offset = seq_db_text(1);
    rsa = createSuffixArray(rT);

    stPos = len(allolaplists);
    F,L,Occ,C = buildFMindex(sa,T);
    FMquery(F,L,Occ,C,sa,T,int(args.nOverlap),0,qseqarray);
    
   
    print('olapcount',stPos);  
    print('olapcount',len(allolaplists));    
    

    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs)); 
 
    stPos = len(allolaplists);
    FMquery(F,L,Occ,C,sa,T,int(args.nOverlap),0,newQseqs);
    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs));
    print(len(allqlabels));
    print(len(allolaplists));
    
    #Now work with reverse reads    
    rF,rL,rOcc,rC = buildFMindex(rsa,rT);

    stPos = len(allolaplists); 
    FMquery(rF,rL,rOcc,rC,rsa,rT,int(args.nOverlap),1,qseqarray);     
    
    newQseqs = addNewQseqs(stPos);
    print(len(newQseqs));
  
    stPos = len(allolaplists);
    FMquery(rF,rL,rOcc,rC,rsa,rT,int(args.nOverlap),1,newQseqs);
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

