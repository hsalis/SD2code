from Bio import SeqIO
import re


def loadFastqFile(filename):

    handle = open(filename)
    records = SeqIO.parse(handle, 'fastq')
    
    for record in records:
        yield record

def generateRegPattern(patternList):
    pattern = re.compile('(' + '|'.join(patternList) + ')')
    return pattern

def doesMatch(readSeq, pattern, pos):
    
    match = re.search(pattern, readSeq)
    if match is None:
        return False
    else:
        return pos
        
def run(filename, patternList, output):

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        if comm.Get_size() > 1:
            use_MPI = True
            from MPI_pool import Pool
            pool = Pool(MPI.COMM_WORLD)
        else:
            use_MPI = False
    except:
        use_MPI = False
        
    print "Using MPI? ", use_MPI
    
    if use_MPI: pool.start()
    
    inputList = []
    
    pattern = generateRegPattern(patternList)
    
    for (i,record) in enumerate(loadFastqFile(filename)):
        inputList.append( (str(record.seq), pattern, i) )
    
    if use_MPI:
        resultList = pool.map(doesMatch, inputList)
    else:
        resultList = map(doesMatch, inputList)
    
    filteredReadIndices = [pos for pos in resultList if pos is not False]
    
    print filteredReadIndices
    #just for now
    if use_MPI: pool.close()
    

if __name__ == "__main__":
    
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
    
    
    patternList = [r1, r3, r6, r7, r9, HHrz5p, HHRz3p, sgRNA_handle, HDVRz]
    run('4342742_rrna_free_reads_unmerged_R1.fastq', patternList, 'testfile.data')
    