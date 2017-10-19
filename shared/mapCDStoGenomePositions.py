from Bio import SeqIO
from Bio.Blast import NCBIWWW as blast
from Bio.Seq import *
import pandas as pd

MAXCOUNTER = 10000

def returnSeqRecordsfromFasta(filename):

    handle = open(filename,'r')
    records = SeqIO.parse(handle,"fasta")
    for record in records:
        yield record

def returnOffTargetSequences(CSVfilename, genomeGenbankFilename):

    PRE_CUTOFF = 50
    POST_CUTOFF = 50

    genomeDict = {}
    handle = open(genomeGenbankFilename)
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        genomeDict[record.id] = record.seq
    handle.close()

    offTargetSites = []
    df = pd.read_csv(CSVfilename, sep = '\t')
    numRows = df.shape[0]

    print "LOG: Reading %s rows in DataFrame" % numRows
    
    for row in range(numRows):
        x = str(df.loc[row,:][0])
        dG = df.loc[row,:][3]
        
        x = x.replace('(','').replace(')','').replace("'","")
        words = x.split(', ')
        
        id = words[0]
        pos = int(words[1])
        strand = 1 #[2]
        
        seq = genomeDict[id]
        offTargetBindingSiteSequence = str(seq[max(0,pos-PRE_CUTOFF) : min(pos+POST_CUTOFF, len(seq) )])
        
        if row > MAXCOUNTER: break
        
    offTargetSites.append( {'sequence' : offTargetBindingSiteSequence, 'locus' : id, 'position' : pos, 'strand' : strand, 'dG_target' : dG} )
    return offTargetSites
    
def matchOffTargetSitesWithCDSs(offTargetSites, CDSFastaFilename):

    for CDSrecord in returnSeqRecordsfromFasta(CDSFastaFilename):
        
        gene = CDSrecord.id.split(" ")[0]
        CDSseq = str(CDSrecord.seq)
        
        print "LOG: Looking for off-target dCas9 binding sites in CDS %s" % gene
        
        for offTargetSite in offTargetSites:
            if offTargetSite['sequence'] in CDSseq or CDSseq in offTargetSite['sequence']:
                #This is a match
                if 'gene' in offTargetSite:
                    offTargetSite['geneList'].append(gene)
                else:
                    offTargetSite['geneList'] = [gene]
                
                print "LOG: Found an off-target site @ %s within CDS %s" % (offTargetSite['position'], gene)
        
    return offTargetSites
   
#BLAST sequence to selected Genome, return positions
# W303 is listed in NCBI's genome database, but there are *NO HITS* to W303 (?!)
def doBlast(seqRecord, genomeName):

    id = seqRecord.id
    gene = id.split(" ")[0]
    
    print "LOG: Running BLAST on gene %s" % gene
    print "LOG: Sequence is %s" % str(record.seq)
    
    result = blast.qblast(
        "blastx",
        "nr",
        str(seqRecord.seq),
        megablast=False,
        expect=1000,
        word_size=7,
        entrez_query=genomeName + '[Organism]',
        format_type = "Text"
        )

#        nucl_reward=1,
#        nucl_penalty=-3,
#        gapcosts="5 2",
    
    handle = open(gene + ".blast",'w')
    handle.write(result.read())
    handle.close()

if __name__ == "__main__":

    CSVfilename = '../dCas9_Circuit_Calculator/YeastGate_dCas9_Binding_Sites.csv'
    genomeGenbankFilename = '../resources/GCA_000292815.1_ASM29281v1_genomic.gbff'
    CDSFastaFilename = '../resources/W303_JRIU00000000_SGD_cds.fsa'
    outputFilename = '../dCas9_Circuit_Calculator/YeastGate_dCas9_Binding_Sites_In_CDSs'
    
    offTargetSites = returnOffTargetSequences(CSVfilename, genomeGenbankFilename)
    offTargetSites = matchOffTargetSitesWithCDSs(offTargetSites, CDSFastaFilename)
    
    df = pd.DataFrame(offTargetSites)
    df.to_csv(outputFilename + '.csv', sep = '\t')
    
    # filename = '../resources/W303_JRIU00000000_SGD_cds.fsa'
#     for record in returnSeqRecordsfromFasta(filename):
#         doBlast(record, 'Saccharomyces cerevisiae W303')
#         break
        
    
    
    