
def crawl_sequence(search_seq_name,search_seq,ref_seq,is_positive,filename="", genomic_offset = 0):
    """
    Does a regular expression search of the reference sequence, looking for instances of search_seq.

    Numbers printed are 1-index, as that is the expected index type that someone would use if they are
    reading genomic data.
    
    genomic_offset is a number that may be added to the coordinates to translate the position within
    the sequence to the position in the genome, if the sequence is a portion of the whole chromosome.
    """
    import re
    
    offset = len(ref_seq)
    
    fmt = search_seq_name + "\t%s\t%d\t%d"
    if is_positive:
        fmt += "\tF"
    else:
        fmt += "\tR"
    fmt  += "\t" + filename
                
    for hit in re.finditer(search_seq,ref_seq):
        (start, end) = (hit.start(), hit.end())
        string_coords = (start, end)
        if not is_positive:
            (end,start) = (offset - start, offset - end)
        (start,end) = (start + genomic_offset + 1, end + genomic_offset + 1)
        print(fmt % (hit.string[string_coords[0]:string_coords[1]],start,end))
    
    
def seq2regex(seq):
    """
    Converts a sequence string to a regular expression query string. 
    Lowercase characters are converted to wildcards.
    R is converted to AG
    Y is converted to CT
    S is converted to GC
    W is converted to AT
    K is converted to GT
    M is converted to AC
    B is converted to CGT
    D is converted to AGT
    H is converted to ACT
    V is converted to ACG
    N is converted to wildcard
    """
    
    search_seq = ""
    for i in seq:
        if i.islower():
            search_seq += "."
        if i in "R":
            search_seq += "[AG]"
        if i in "Y":
            search_seq += "[CT]"
        if i in "S":
            search_seq += "[GC]"
        if i in "W":
            search_seq += "[AT]"
        if i in "K":
            search_seq += "[GT]"
        if i in "M":
            search_seq += "[AC]"
        if i in "B":
            search_seq += "[CGT]"
        if i in "D":
            search_seq += "[AGT]"
        if i in "H":
            search_seq += "[ACT]"
        if i in "V":
            search_seq += "[ACG]"
        if i in "N":
            search_seq += "."
        else:
            search_seq += i
    return search_seq
    


def search(infile, search_seqs,should_forward_search = True, should_revcomp = True, genomic_offset = 0):
    from Bio import SeqIO
    for rec in SeqIO.parse(infile,"fasta"):
        rev_comp_str = None
        for (search_seq_name,search_seq) in search_seqs:
            if should_forward_search:
                crawl_sequence(search_seq_name,search_seq,str(rec.seq),True,infile.name,genomic_offset)
            if should_revcomp:
                if not rev_comp_str:
                    rev_comp_str = str(rec.seq.reverse_complement())
                crawl_sequence(search_seq_name,search_seq,rev_comp_str,False,infile.name,genomic_offset)
                