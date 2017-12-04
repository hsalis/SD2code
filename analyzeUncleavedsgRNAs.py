from Bio import SeqIO
import re, itertools, sys
import difflib

MININUM_MATCH_THRESHOLD = 20 #nucleotides
PATH_TO_FILES = "../"
global PATH_TO_FILES

def loadFastqFile(filename):

    handle = open(PATH_TO_FILES + filename)
    records = SeqIO.parse(handle, 'fastq')
    
    for record in records:
        yield record
    
    handle.close()

def generateRegPattern(patternList):
    pattern = re.compile('(' + '|'.join([x.upper() for x in patternList]) + ')')
    return pattern

def findStringMatches(readSeq, searchString, minimum_threshold):

    matchList = difflib.SequenceMatcher(None, readSeq, searchString).get_matching_blocks()
    if len(matchList) > 1:
        for match in matchList:
            if match[2] > minimum_threshold: return True
    
    return False

def doesMatch( (record_5p, record_3p, pattern, pos) ):
    readSeq = str(record_5p.seq) + str(record_3p.seq)
    match = re.search(pattern, readSeq)
    if match is None:
        return False
    else:
        return pos

def doesMatch2( (record_5p, record_3p, patternList, pos) ):
    readSeq = str(record_5p.seq) + str(record_3p.seq)
    matchList = []
    for pattern in patternList:
        if findStringMatches(readSeq, pattern, MININUM_MATCH_THRESHOLD):
            matchList.append( True )
        else:
            matchList.append( False )

    return matchList
                
def run(filenamePrefix, patternList):
    
    handle1 = open(PATH_TO_FILES + 'circuitFiltered_' + filenamePrefix + '_R1' + '.fastq','w')
    handle2 = open(PATH_TO_FILES + 'circuitFiltered_' + filenamePrefix + '_R2' + '.fastq','w')
    
    pattern = generateRegPattern(patternList)
    
    for i, (record_5p, record_3p) in enumerate(itertools.izip( loadFastqFile(filenamePrefix + '_R1' + '.fastq'), loadFastqFile(filenamePrefix + '_R2' + '.fastq') )):
        #print "Reading read #%s with sequences %s and %s" % (i, str(record_5p.seq), str(record_3p.seq))    
        print "Reading read #%s" % i
        result = doesMatch( (record_5p, record_3p, pattern, i) )
        if result is not False:
            SeqIO.write(record_5p, handle1, 'fastq')
            SeqIO.write(record_3p, handle2, 'fastq')

    handle1.close()
    handle2.close()

def analyze(filenamePrefix, patternList):
    
    handle1 = open(PATH_TO_FILES + 'circuitFiltered_' + filenamePrefix + '_R1' + '.fastq')
    handle2 = open(PATH_TO_FILES + 'circuitFiltered_' + filenamePrefix + '_R2' + '.fastq')
    
    patternCounts = {}
        
    for pattern in patternList:
        patternCounts[pattern] = 0
        pattern = generateRegPattern([pattern])
    
        for i, (record_5p, record_3p) in enumerate(itertools.izip( loadFastqFile(filenamePrefix + '_R1' + '.fastq'), loadFastqFile(filenamePrefix + '_R2' + '.fastq') )):
            #print "Reading read #%s with sequences %s and %s" % (i, str(record_5p.seq), str(record_3p.seq))    
            print "Reading read #%s" % i
            result = doesMatch( (record_5p, record_3p, pattern, i) )
            if result is not False:
                patternCounts[pattern] += 1

    handle1.close()
    handle2.close()
    return patternCounts
    
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        PATH_TO_FILES = sys.argv[1]    
    
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        if comm.Get_size() > 1:
            use_MPI = True
            from MPI_pool import Pool
            pool = Pool(MPI.COMM_WORLD)
        else:
            use_MPI = False
            pool = None
    except:
        use_MPI = False
        pool = None
        
    print "Using MPI? ", use_MPI
    
    #sgRNA guide sequences in "Digital logic circuits in yeast with CRISPR-dCas9 NOR gates" by Gander et. al.
    r1 =  'GGAACGTGATTGAATAACTT'
    r2 =  'ACCAACGCAAAAAGATTTAG'
    r3 =  'CATTGCCATACACCTTGAGG'
    r4 =  'GAAAATCACAACTCTACTGA'
    r5 =  'GAAGTCAGTTGACAGAGTCG'
    r6 =  'GTGGTAACTTGCTCCATGTC'
    r7 =  'CTTTACGTATAGGTTTAGAG'
    r8 =  'CGCATTTCCTATTCAAACTT'
    r9 =  'GCAACCCACAAATATCCAGT'
    r10 = 'GTGACATAAACATTCGACTC'
    r11 = 'GGGCAAAGAGACGCTTGTCG'
    r12 = 'GAAGTCATCGCTTCTTGTCG'
    r13 = 'GAGTTGACAAAGTATAACTT'
    r14 = 'GAAGTTTCAGAATCTCGACG'
    r15 = 'GGCTAGGATCCATCTGACTT'
    r16 = 'GCAACCATAGACTCTCCAGG'
    r17 = 'ACCACAACTGAGTCGAACCT'
    r18 = 'GGGTAGCAACACTCGTACTT'
    r19 = 'GTAAAAGATAACTCTGTTGC'
    r20 = 'TCTACCCGAGACTCAAACGG'
    
    RGRi = 'ggattctagaactagtggatctacaaaNNNNNNctgatgagtccgtgaggacgaaacgagtaagctcgtcNNNNNNNNNNNNNNNNNNNNgttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttggccggcatggtcccagcctcctcgctggcgccggctgggcaacatgcttcggcatggcgaatgggactgataccgtcgacctcgagtc'
    
    HHRz5p = 'ggattctagaactagtggatctacaaa'
    HHRz3p = 'ctgatgagtccgtgaggacgaaacgagtaagctcgtc'
    sgRNA_handle = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtgg'
    HDVRz = 'caccgagtcggtgcttttggccggcatggtcccagcctcctcgctggcgccggctgggcaacatgcttcggcatggcgaatgggactgataccgtcgacctcgagtc'
    
    
    
    
    basePatternList = [r1, r7, r9, HHRz5p, HHRz3p, sgRNA_handle, HDVRz]
    
    runDict = {}
    runDict['XOR_00_b1_t1'] = {'filename' : '4342742_rrna_free_reads_unmerged',
                               'patternList' : basePatternList}
    runDict['XOR_00_b1_t2'] = {'filename' : '4342743_rrna_free_reads_unmerged',
                                'patternList' : basePatternList}
    
    # runDict['XOR_00_b2_t1'] = {'filename' : '4342750_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList}
    # runDict['XOR_00_b2_t2'] = {'filename' : '4342751_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList}
    # runDict['XOR_00_b3_t1'] = {'filename' : '4342758_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList}
    # runDict['XOR_00_b3_t2'] = {'filename' : '4342759_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList}
    
    ====
    
    # runDict['XOR_01_b1_t1'] = {'filename' : '4342744_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    # runDict['XOR_01_b1_t2'] = {'filename' : '4342745_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    # runDict['XOR_01_b2_t1'] = {'filename' : '4342752_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    # runDict['XOR_01_b2_t2'] = {'filename' : '4342753_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    # runDict['XOR_01_b3_t1'] = {'filename' : '4342760_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    # runDict['XOR_01_b3_t2'] = {'filename' : '4342761_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r6}
    
    ====
    
    
    # runDict['XOR_10_b1_t1'] = {'filename' : '4342746_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}
    # runDict['XOR_10_b1_t2'] = {'filename' : '4342747_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}
    # runDict['XOR_10_b2_t1'] = {'filename' : '4342754_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}
    # runDict['XOR_10_b2_t2'] = {'filename' : '4342755_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}
    # runDict['XOR_10_b3_t1'] = {'filename' : '4342762_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}
    # runDict['XOR_10_b3_t2'] = {'filename' : '4342763_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3}

    ====

                                
    # runDict['XOR_11_b1_t1'] = {'filename' : '4342748_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    # runDict['XOR_11_b1_t2'] = {'filename' : '4342749_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    # runDict['XOR_11_b2_t1'] = {'filename' : '4342756_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    # runDict['XOR_11_b2_t2'] = {'filename' : '4342757_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    # runDict['XOR_11_b3_t1'] = {'filename' : '4342764_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    # runDict['XOR_11_b3_t2'] = {'filename' : '4342765_rrna_free_reads_unmerged',
                                # 'patternList' : basePatternList + r3 + r6}
    
    #====
    
    if use_MPI: pool.start()
    
    allPatternCounts = {}
    
    for (data, info) in runDict.items():
	    print info
        allPatternCounts[info['filename']] = analyze(info['filename'], info['patternList'])
    
    print allPatternCounts
    
    if use_MPI: pool.close()
