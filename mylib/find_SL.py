import pysam, argparse
from Bio.Seq import Seq
from subprocess import Popen, PIPE
import tqdm
import sys


 
#https://www.biostars.org/p/1890/
def bam_read_count(bamfile):
    """ Return a tuple of the number of mapped and unmapped reads in a bam file
        Bam file reference sorted and indexed
    """
    p = Popen(['samtools', 'idxstats', bamfile], stdout=PIPE)
    mapped = 0
    unmapped = 0
    for line in p.stdout:
        rname, rlen, nm, nu = line.rstrip().split()
        mapped += int(nm)
        unmapped += int(nu)
    #p.close()
    return (mapped, unmapped)


def test_none(barcode_region):
    '''
    function to test how many bases are
    soft-clipped in the barcode sequence. 
    the input is the output of .get_aligned_pairs()
    if the base is soft-clipped the second item 
    of the tuple is None
    '''
    #print(barcode_region)
    tot = [1 for n in barcode_region if str(n[1])=='None'] 
    #print(tot)
    #print(len(barcode_region))
    coverage = float(sum(tot)) / len(barcode_region)
    #print(coverage)
    #hard trheshold, we want at least 70% of 
    #barcode sequence soft clipped
    if coverage > 0.8:
        return True
    else:
        return False

def get_start(alignment, orient, barcode):
    #this return a tuple list for each base in the read
    #the position in the read and the corresponding genome position
    #genome position is None is softclipped
    pairs = alignment.get_aligned_pairs()
    #print(alignment.query_name)
    #print(pairs)
    start_barcode = alignment.seq.find(barcode)
    end_barcode = start_barcode+len(barcode)
    #extract only the aligned_pairs for the barcode
    barcode_region = pairs[start_barcode:end_barcode]
    #print(barcode_region)
    #print(
    #alignment.is_secondary, alignment.is_supplementary,
    #alignment.is_paired, alignment.is_reverse, alignment.is_qcfail, 
    #alignment.mapping_quality)
    if orient == 'F':
        #if the barcode is in forward orientation
        #the genome sequence start at the end of the
        #barcode sequence
        start = start_barcode+len(barcode)

    elif orient == 'R':
        #if the barcode is in reverse complement orientation
        #the genome sequence start at start of the
        #barcode sequence       
        start = start_barcode

    else:
        raise('error')
    try:
        #get finally the starting position
        #respect to the genome sequence
        #we remove one as python index from 0
        start = pairs[start-1][1]
    except:
        print(start, pairs)
        raise('error')
    return start, test_none(barcode_region)

def parse(bamfile, all_reads, barcode, barcode_rc):
    count=0
    for alignment in tqdm.tqdm(bamfile, total=all_reads, miniters=1000000):
        count+=1
        if alignment.mapping_quality >= 10:
            ref_chr = bamfile.get_reference_name(alignment.reference_id)
            #now we do something different if the barcode is found
            #forward or reverse complement
            #print(alignment.name)
            if barcode in str(alignment.seq):
                start, test = get_start(alignment, 'F', barcode)
                yield ','.join([ref_chr, str(start), str(test),'F'])
                
            elif barcode_rc in str(alignment.seq):
                start, test = get_start(alignment, 'R', barcode_rc)
                yield ','.join([ref_chr, str(start), str(test),'R']) 
            else:
                continue

            #print(start, ref_chr, test)
            #if count>=1:
                #break

def main(infile, barcode):

    outfile = infile.split('.')[0]+'_count.csv'

    with open(outfile, 'w') as outfile, pysam.Samfile(infile,"rb") as bamfile:

        outfile.write('chr,pos,code,orient\n')

        mapped, unmapped = bam_read_count(infile)
        print('parsing {mapped} mapped and {unmapped} unmapped reads'.format(
            mapped=mapped, unmapped=unmapped))
        
        all_reads = mapped + unmapped

        
        barcode_rc = str(Seq(barcode).reverse_complement())

        for n in parse(bamfile, all_reads, barcode, barcode_rc):
            outfile.write(n+'\n')

        outfile.close()
        bamfile.close()   





if __name__ == '__main__':
    
    #splice leader sequences
    #https://doi.org/10.1093/nar/gkq237
    #The last 14 nt (TCTGTACTATATTG) of the spliced leader (SL) sequence
    #  (AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG) are unique
    #_SLF = 'TCTGTACTATATTG'
    #_SLR = 'CAATATAGTACAGA'
    
    #iRNA barcode
    #_SLF = GCCTCGCGA
    #_SLR = TCGCGAGGC    

    infile = sys.argv[1]
    barcode = sys.argv[2]
    main(infile, barcode)
    print('done')