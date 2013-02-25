# gene_extractor -g RASSF4:assembly_files/hs_ref_GRCh37.p5_chr10.gbs -o RASSF4.fa ../../genomes/hs/hs_ref_GRCh37.p5_chr10.fa
CHRY_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/gbs/hs_ref_GRCh37.p5_chrY.gbs.gz"
CHRY_FA_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/seq/hs_ref_GRCh37.p5_chrY.fa.gz"
CHRY_GZ_FILENAME = "chrY.gbs.gz"
CHRY_FILENAME = "chrY.gbs"
CHRY_FA_GZ_FILENAME = "chrY.fa.gz"
CHRY_FA_FILENAME = "chrY.fa"
TEST_GENE = "IL9R"
GENE_EXTRACT = "IL9R.fa"
CANONICAL_RESULTS = "IL9R.hits"

def setup():
    from os import chdir
    import os.path as OP
    
    def report_hook(*args):
        from sys import stdout
        stdout.write(".")
    
    chdir("test")
    for (url,name) in [(CHRY_URL,CHRY_GZ_FILENAME),(CHRY_FA_URL,CHRY_FA_GZ_FILENAME)]: 
        if not OP.isfile(name):
            from urllib import urlretrieve
            print("Downloading %s" % name)
            urlretrieve(url,name,report_hook)
            print("")

def gunzip_chromosome():
    import os.path as OP
    import gzip
    
    if not OP.isfile(CHRY_GZ_FILENAME):
        raise Exception("Could not download %s" % CHRY_GZ_FILENAME)
    for (gzfile, bigfile) in [(CHRY_GZ_FILENAME,CHRY_FILENAME), (CHRY_FA_GZ_FILENAME,CHRY_FA_FILENAME)]:
        zipped = gzip.open(gzfile, "r")
        unzipped = open(bigfile,"w")
        unzipped.write(zipped.read())
    
def extract_gene():
    from gene_crawler import Assembly
    import gene_crawler.gene_extractor as GE
    
    asmb = Assembly(CHRY_FILENAME)
    gene = asmb.get_gene(TEST_GENE)
    (start,end) = gene.get_coords()
    print("Found %s at %s" % (gene.name,(start+1,end+1))) # human readable output = one-indexed
    chry_file = open(CHRY_FA_FILENAME,"r")
    for line in chry_file:
        if line.strip()[0] == ">":
            break
    test_gene = open(GENE_EXTRACT,"w")
    test_gene.write(">TEST GENE %s\n" % TEST_GENE)
    GE.extract_gene(chry_file, test_gene, start, end-start)
    
    return start
    
def find_sites(genomic_offset):
    from gene_crawler import gene_searcher
    
    print("Search Results:")
    gene_searcher.search(open(GENE_EXTRACT,"r"),[("AP1a","TGACTCA"), ("AP1b","TGAGTCA")],genomic_offset = genomic_offset)
    
    print("Expected Results:")
    print(open(CANONICAL_RESULTS,"r").read())
    
def run():
    setup()
    gunzip_chromosome()
    genomic_offset = extract_gene()
    find_sites(genomic_offset + 1) # the function assume the *user* supplied the offset, thus it is one-indexed
        
if __name__ == "__main__":
    run()