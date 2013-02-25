
def crawl_sequence(search_seq_name,search_seq,ref_seq,is_positive,filename="", genomic_offset = 0, output = None):
    """
    Does a regular expression search of the reference sequence, looking for instances of search_seq.

    Numbers printed are 1-index, as that is the expected index type that someone would use if they are
    reading genomic data.
    
    genomic_offset is a number that may be added to the coordinates to translate the position within
    the sequence to the position in the genome, if the sequence is a portion of the whole chromosome.
    """
    import re
    
    if not output:
        from sys import stdout
        output = stdout
    
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
            
        # Add offset. If the user supplied one, use it. If not, add 1 to be one-indexed.
        if genomic_offset != 0:
            start += genomic_offset
            end += genomic_offset
        else:
            start += 1
            end += 1
        
        output.write(fmt % (hit.string[string_coords[0]:string_coords[1]],start,end))
        output.write("\n")
    
    
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
        elif i in "R":
            search_seq += "[AG]"
        elif i in "Y":
            search_seq += "[CT]"
        elif i in "S":
            search_seq += "[GC]"
        elif i in "W":
            search_seq += "[AT]"
        elif i in "K":
            search_seq += "[GT]"
        elif i in "M":
            search_seq += "[AC]"
        elif i in "B":
            search_seq += "[CGT]"
        elif i in "D":
            search_seq += "[AGT]"
        elif i in "H":
            search_seq += "[ACT]"
        elif i in "V":
            search_seq += "[ACG]"
        elif i in "N":
            search_seq += "."
        else:
            search_seq += i
    
    return search_seq
    
def get_genome_offset(_string):
    from gene_crawler import FormatError
    padding_suffix = "bp_padding"
    string = _string[1:]
    tokens = string.split("|")
    if len(tokens) != 3:    
        raise FormatError("Invalid header: %s" % _string)
    
    padding_str = tokens[1]
    if not padding_suffix in tokens[1]:
        raise FormatError("Invalid header: %s" % _string)
    padding_idx = padding_str.replace(padding_suffix,"")
    padding = int(padding_idx)
    
    if not ".." in tokens[2]:
        raise FormatError("Invalid header: %s" % _string)
    coords = tokens[2].split("..")
    
    offset =  int(coords[0]) - padding
    if offset < 0:
        offset = 0
    return offset
    
        
    
    
def search(infile, search_seqs,should_forward_search = True, should_revcomp = True, genomic_offset = 0, output = None):
    """
    Searches Gene sequences for instances of each target sequences in the search_seqs list.
    
    search_seqs is a lsit of tuples, where the first element is the sequence name and the second element is the actual sequence
    
    genomic_offset is an integer that may be added to the output coordinates if the gene sequence is a segment taken from the genome.
    If genomic_offset is equal to -1, the offset will be extracted from the sequence header, assuming the header has the following format:
    ">GENENAME|Nbp padding|START..END"
    where N is the base pair padding size and START and END are the coordinates of the gene.
    """
    from Bio import SeqIO
    for rec in SeqIO.parse(infile,"fasta"):
        sequence_offset = genomic_offset
        if genomic_offset == -1:
            sequence_offset = get_genome_offset(rec.id)
                
        rev_comp_str = None
        for (search_seq_name,search_seq) in search_seqs:
            if should_forward_search:
                crawl_sequence(search_seq_name,search_seq,str(rec.seq),True,infile.name,sequence_offset,output)
            if should_revcomp:
                if not rev_comp_str:
                    rev_comp_str = str(rec.seq.reverse_complement())
                crawl_sequence(search_seq_name,search_seq,rev_comp_str,False,infile.name,sequence_offset,output)
                