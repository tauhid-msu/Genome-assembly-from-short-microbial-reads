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

def updateBackward(F,L,Occ,C,c_index,l,u):
    
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

def outputFirst10overlap(i,j,olaplen):
    if firstolaps<10:
        file_top10olap.write(">%s %s %d\n"%(qseqlabels[i],seqlabels[j],olaplen ) );
        file_top10olap.write("%s\n"%(qseqarray[i]))
        olapstring = "";
        #blanksp = len(qseqarray[i]) - j;
        blanksp = len(qseqarray[i])-olaplen+1;
        for k in range(blanksp):
            olapstring+="_";
        olapstring+=seqarray[j];
        #print(olapstring);
        file_top10olap.write("%s\n"%(olapstring)); 
        firstolaps += 1;


def outputOverlap(i,l,j,isqprefix):
    #print(i,l,j) 
    db_i = findDbseqid(l,isqprefix);
    curOlap = [];
    curOlap.append(i);
    curOlap.append(db_i);
    olaplen = len(qseqarray[i])-j-1;
    curOlap.append(olaplen);
    curOlap.append(isqprefix);  
    allolaplists.append(curOlap);
    #file_allolap.write("%s %s %d\n"%(qseqlabels[i],seqlabels[db_i],(len(qseqarray[i])-j-1) ) )
    #outputFirst10overlap(i,db_i,(len(qseqarray[i])-j));
           
    
def FMquery(F,L,Occ,C,ssa,T,th,isqprefix):
    #print("FMquery");
    dictionary = {"$":0,"A":1,"C":2,"G":3,"T":4}; 
    nqseq = qseqarray.shape[0];
    for i in range(nqseq):
        qlen = len(qseqarray[i])-1;
        j = qlen;
        curQseq = qseqarray[i];
        if isqprefix==1:     
            #print('>',curQseq);
            curQseq = ''.join(reversed(curQseq));    
            #print('<',curQseq);

        min_j = qlen; #min_j refers maximum overlap
        mxolap_ssa_index = 0;
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
                l1,u1 = updateBackward(F,L,Occ,C,0,l,u);
                if l1<=u1:
                    #print('lu$',l1,u1,'lu',l,u);
                    for rr in range(l,u+1,1):
                        if ssa[rr]>0:
                            #print(rr,ssa[rr],T[ssa[rr]-1]);
                            if T[ssa[rr]-1]=='$':
                                min_j = j;
                                mxolap_ssa_index = ssa[rr];
                                #outputOverlap(i,ssa[rr],j,isqprefix);
                        else:
                            #print(rr,ssa[rr],'$');
                            min_j = j;
                            mxolap_ssa_index = 0;
                            #outputOverlap(i,0,j,isqprefix);
                    

            c = curQseq[j]; 
            c_index = dictionary[c];
            l,u = updateBackward(F,L,Occ,C,c_index,l,u);
            #print(c,c_index,'lu',l,'->',u);
            j -= 1;
        #print the maximum overlap
        if min_j<qlen:          
            outputOverlap(i,mxolap_ssa_index,min_j,isqprefix);

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
    #seqarray = seqarray[0:args.nSeq];
    #seqlabels = seqlabels[0:args.nSeq];
    #print(seqarray.shape);

    global qseqarray
    global qseqlabels
    qseqarray,qseqlabels = read_seqfile(args.qseqfile);

    #print('\n'.join(bwm('abaaba$')))
    #text = str(seqarray[0]);
    text = "";
    global dbseq_offset;
    dbseq_offset = [];
    for i in range(seqarray.shape[0]):
        dbseq_offset.append(len(text));
        text +=seqarray[i];
        if i<(seqarray.shape[0]-1):
            text += "$";	
    dbseq_offset.append(len(text)-1);
    #print(dbseq_offset);
    T = text.strip()


    myT = []
    for chr in T:
        myT.append(ord(chr))
    sa = ksa(myT)
    #print T
    #print ["{0:^2}".format(chr) for chr in T]
    #print ["{0:^2}".format(i) for i, chr in enumerate(T)]
    #print sa
   
    
    global firstolaps 
    global allolaplists
    allolaplists = [[]]   
    firstolaps = 0;

    F,L,Occ,C = buildFMindex(sa,T);
    #print(C);
    #for i in range(len(sa)):
    #    print(i,sa[i],F[i],L[i],Occ[i]);
    FMquery(F,L,Occ,C,sa,T,int(args.nOverlap),0);


    rtext = "";
    global rdbseq_offset;
    rdbseq_offset = [];
    for i in range(seqarray.shape[0]):
        rdbseq_offset.append(len(rtext));
        rseq = ''.join(reversed(seqarray[i]));
        rtext +=rseq;
        if i<(seqarray.shape[0]-1):
            rtext += "$";
    rdbseq_offset.append(len(rtext)-1);
    #print(rdbseq_offset);
    rT = rtext.strip()

    myrT = []
    for chr in rT:
        myrT.append(ord(chr))
    rsa = ksa(myrT)

    #print rT
    #print ["{0:^2}".format(chr) for chr in rT]
    #print ["{0:^2}".format(i) for i, chr in enumerate(rT)]
    #print rsa

    rF,rL,rOcc,rC = buildFMindex(rsa,rT);
    #print(rC);
    #for i in range(len(rsa)):
    #    print(i,rsa[i],rF[i],rL[i],rOcc[i]);
    FMquery(rF,rL,rOcc,rC,rsa,rT,int(args.nOverlap),1);     


    file_allolap = open(args.outpath+"/overlap_list","w");
    for i in range(1,len(allolaplists),1):
        if allolaplists[i][3]==0:
            file_allolap.write("%s %s %d\n"%(qseqlabels[allolaplists[i][0]],seqlabels[allolaplists[i][1]], allolaplists[i][2] ) );
        elif allolaplists[i][3]==1:
            file_allolap.write("%s %s %d\n"%(seqlabels[allolaplists[i][1]],qseqlabels[allolaplists[i][0]], allolaplists[i][2] ) );
    file_allolap.close();
    
    file_top10olap = open(args.outpath+"/first_ten_overlaps","w")
    for i in range(1,len(allolaplists),1):
        if i<=10:
            olapstring = "";
            #query overlap as suffix
            if allolaplists[i][3]==0: 
                file_top10olap.write(">%s %s %d\n"%(qseqlabels[allolaplists[i][0]],seqlabels[allolaplists[i][1]],allolaplists[i][2]));
                file_top10olap.write("%s\n"%(qseqarray[allolaplists[i][0]]));
                blanksp = len(qseqarray[allolaplists[i][0]])-allolaplists[i][2];
                for k in range(blanksp):
                    olapstring+=" ";
                olapstring+=seqarray[allolaplists[i][1]]; 
            #query overlap as prefix
            elif allolaplists[i][3]==1:
                 file_top10olap.write(">%s %s %d\n"%(seqlabels[allolaplists[i][1]],qseqlabels[allolaplists[i][0]],allolaplists[i][2]));
                 file_top10olap.write("%s\n"%(seqarray[allolaplists[i][1]]));
                 blanksp = len(seqarray[allolaplists[i][1]])-allolaplists[i][2];
                 for k in range(blanksp):
                     olapstring+=" ";
                 olapstring+=qseqarray[allolaplists[i][0]]; 

            #blanksp = len(qseqarray[allolaplists[i][0]])-allolaplists[i][2];
            #for k in range(blanksp):
            #    olapstring+=" ";
            #olapstring+=seqarray[allolaplists[i][1]];
            file_top10olap.write("%s\n\n"%(olapstring));
            firstolaps += 1; 
    file_top10olap.close(); 	

if __name__ == '__main__':
    main()

