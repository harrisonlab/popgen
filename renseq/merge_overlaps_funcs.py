import numpy as np

def generate_kmers(seq,size):
    '''
    generate all the kmers from the sequencfe
    '''
    
    for i in xrange(len(seq)-size+1): yield seq[i:i+size]
        

def generate_baits(seq,size):
    '''
    chop seq up into tiled baits of the required size
    exclude the final bait if it's not full length
    '''
    
    nseq = len(seq)
    for i in xrange(0,nseq,size):
        bait = seq[i:i+size]
        if len(bait) != size: continue #exclude partial baits
        yield bait

def cleanup_seq(seq):
    '''
    convert to upper case
    check for unexpected symbols
    split at Ns
    '''
    
    #uppercase
    seq = seq.upper()
    
    #check for unexpected characters
    allow = 'ATCGN' #allowed symbols
    chk = [x for x in seq if x in allow]
    assert len(chk) == len(seq)
    
    #split at Ns
    seq_list = [x for x in seq.split('N') if x != '']
    
    return seq_list

def merge_overlaps(l):
    '''
    input is a list of [start,end] positions
    all assumed to be on the same reference sequence
    assuming python slice notation:
    start is the first position, end is one-past-the-end
    merging is triggered when any two intervals overlap
    or are immediately adjacent without any gap
    ie when start_a <= end_b or start_b <= end_a
    checks that start <= end for all pairs
    
    returns new merged list but does alter the original
    '''

    n = len(l)
    
    for i in xrange(n-1):
        hit_i = l[i]
        assert hit_i[0] <= hit_i[1]
        for j in xrange(i+1,n):
            hit_j = l[j]
            assert hit_j[0] <= hit_j[1]

            #  endi      startj     endj      starti
            if hit_i[1] < hit_j[0] or hit_j[1] < hit_i[0]: continue #not touching

            #merge i into j
            if hit_i[0] < hit_j[0]: hit_j[0] = hit_i[0]  #new start
            if hit_i[1] > hit_j[1]: hit_j[1] = hit_i[1]  #new end
            l[i] = None #flag for removal

            #skip any further comparisons with hit_i
            break

    #remove blank records
    return [x for x in l if x != None]

def ungapped_smithwaterman(seq1,seq2):
    '''
    ungapped smith waterman alignment
    score +1 if bases identical else 0
    '''
    
    n1,n2 = len(seq1),len(seq2)
    rows,cols = n1+1,n2+1
    
    #score matrix: H[row][col]
    H = np.zeros([rows,cols])

    #store global max H value
    Hmax = -1.0
    imax,jmax = None,None

    #calculate max score matrix
    for i in xrange(1,n1+1):#row/seq1
        for j in xrange(1,n2+1):#col/seq2
            if seq1[i-1] == seq2[j-1]:
                H[i][j] = H[i-1][j-1] + 1
            else:
                H[i][j] = H[i-1][j-1]

            if H[i][j] > Hmax:
                Hmax = H[i][j]
                imax,jmax = i,j
            
    #construct an optimal alignment from the traceback
    i,j = imax,jmax
    align1,align2 = '',''
    midline = ''
    while True:
        if H[i][j] == 0.0:
            #end of alignment
            break

        #diagonal
        align1 = align1 + seq1[i-1]
        align2 = align2 + seq2[j-1]
        if seq1[i-1] == seq2[j-1]:
            midline = midline + '|'
        else:
            midline = midline + ' '
        i -=1
        j -=1
            

    #switch alignments into correct orientation
    align1  = align1[::-1]
    midline = midline[::-1]
    align2  = align2[::-1]

    #print align1
    #print midline
    #print align2
    #print Hmax

    return align1,midline,align2,Hmax

def ungapped_smithwaterman_score(seq1,seq2):
    '''
    ungapped smith waterman alignment
    score +1 if bases identical else 0
    return score only
    '''
    
    n1,n2 = len(seq1),len(seq2)
    rows,cols = n1+1,n2+1
    
    #score matrix: H[row][col]
    H = np.zeros([rows,cols])

    #store global max H value
    Hmax = -1.0

    #calculate max score matrix
    for i in xrange(1,n1+1):#row/seq1
        for j in xrange(1,n2+1):#col/seq2
            if seq1[i-1] == seq2[j-1]:
                H[i][j] = H[i-1][j-1] + 1
            else:
                H[i][j] = H[i-1][j-1]

            if H[i][j] > Hmax:
                Hmax = H[i][j]
            
    return Hmax

def ungapped_smithwaterman_score_prealloc(seq1,seq2,H):
    '''
    ungapped smith waterman alignment
    score +1 if bases identical else 0
    return score only
    H must be preallocated
    '''
    
    n1,n2 = len(seq1),len(seq2)
    rows,cols = n1+1,n2+1
    
    #score matrix: H[row][col]
    H.fill(0.0)

    #store global max H value
    Hmax = -1.0

    #calculate max score matrix
    for i in xrange(1,n1+1):#row/seq1
        for j in xrange(1,n2+1):#col/seq2
            if seq1[i-1] == seq2[j-1]:
                H[i][j] = H[i-1][j-1] + 1.0
            else:
                H[i][j] = H[i-1][j-1]

            if H[i][j] > Hmax:
                Hmax = H[i][j]
            
    return Hmax

def exhaustive_search(seq1,seq2):
    '''
    generate all valid local sequence alignments
    retain the one with the maximum score
    intended for use testing smith waterman
    '''
    
    n1 = len(seq1)
    n2 = len(seq2)
    best = -1
    best_shift = None
    
    for shift in xrange(-(n2+n1+2),n2+n1+2):
        score = 0
        for i in xrange(n1):
            j = i + shift
            if j < 0 or j >= n2: continue
            if seq1[i] == seq2[j]: score += 1
        if score > best:
            best = score
            best_shift = shift
        
    return best,best_shift
