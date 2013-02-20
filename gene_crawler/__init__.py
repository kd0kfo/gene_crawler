__all__ = ["gene_crawler"]


DEFAULT_ASSEMBLY_KEYS = ["locus","definition","source","organism",]
GBS_KEY_LEN = 12

class FormatError(Exception):
    pass

def load_keys(_file,key_list):
    """
    Takes a list of keys, searches the given file and returns a dict containing each 
    key's value.
    """
    
    retval = {}
    
    infile = _file
    if isinstance(_file,str):
        infile = open(_file,"r")
    
    last_key = None
    for line in infile:
        if line[-1] == '\n':
            line = line[0:-1]
        key = line[0:GBS_KEY_LEN].strip().lower()
        if not key:
            if last_key:
                retval[last_key] += "\n" + line[GBS_KEY_LEN:]
            continue
        if key in key_list:  
            retval[key] = line[GBS_KEY_LEN:] 
            
    return retval

class Gene():
    def __init__(self,loc = None, name = None):
        self.loc = loc
        self.name = name
        
    def get_coords(self):
        import re
        if not self.loc:
            return None
        locstr = self.loc
        for i in [">","<"]:
            locstr.replace(i,"")
        coords = re.search(r"(\d*)\.\.(\d*)",locstr)
        if not coords:
            raise FormatError("Invalid coordinate: '%s'" % self.loc)
        coords = coords.groups()
        if len(coords) != 2:
            raise FormatError("Invalid coordinate: '%s'" % self.loc)
        return (int(coords[0]),int(coords[1]))

def str2gene(string):
    string = string.strip()
    if not string:
        return None
    tokens = string.split("/")
    retval = Gene()
    if len(tokens[0]):
        retval.loc = tokens[0]
    return retval

class Assembly():
    def __init__(self,filename = None):
        self.data = {}
        self.filename = filename
        
        if filename:
            self.data = load_keys(filename,DEFAULT_ASSEMBLY_KEYS)
            
    def get_gene(self,gene_name):
        def extract_gene(string):
            gene = str2gene(string)
            if gene:
                gene.name = gene_name
            return gene 
        
        if not self.filename:
            return None
        
        have_gene = False
        gene_info = ""
        with open(self.filename,"r") as infile:
            for line in infile:
                if line[-1] == '\n':
                    line = line[0:-1]
                key = line[0:GBS_KEY_LEN].strip().lower()
                if not key: # Data within a key block
                    if have_gene: # We are in a gene block, concat the data
                        gene_info += line[GBS_KEY_LEN:]
                    continue
                if key == "gene": # Start of a gene block
                    if have_gene: # If this is true, two gene blocks are next to each other. Is this the one we are looking for?
                        if '/gene="%s"' % gene_name in gene_info:
                            return extract_gene(gene_info)
                    have_gene = True
                    gene_info = line[GBS_KEY_LEN:]
                    continue
                if have_gene:# In this case, a gene block just ended, because another block is starting.
                    if '/gene="%s"' % gene_name in gene_info:
                        return extract_gene(gene_info)
                    have_gene = False
                    gene_info = ""
                     
        return None
                    
        
        
            