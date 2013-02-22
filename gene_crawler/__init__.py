__all__ = ["Gene","Assembly"]


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
    """
    Gene Abstraction
    
    Stores Gene data.
    """
    def __init__(self,loc = "", name = "", synonym = "",note=""):
        self.loc = loc
        self.name = name
        self.synonym = synonym
        self.note = note
        self.direction = "+"
        if "complement" in loc:
            self.direction = "-"
        
    def get_coords(self):
        """
        Parses the location data and returns a tuple containing the beginning and end points
        of the gene
        """
        import re
        if not self.loc:
            return None
        locstr = self.loc
        for i in [">","<"]:
            locstr = locstr.replace(i,"")
        coords = re.search(r"(\d*)\.\.(\d*)",locstr)
        if not coords:
            raise FormatError("Invalid coordinate: '%s'" % self.loc)
        coords = coords.groups()
        if len(coords) != 2:
            raise FormatError("Invalid coordinate: '%s'" % self.loc)
        try:
            return (int(coords[0]) - 1,int(coords[1]) - 1) # GBS files, and therefore self.loc values, are one-indexed. Internally, we want to use zero-indexed values
        except ValueError as ve:
            print("Invalid coordinates: %s (%s)" % (locstr,coords))
            raise ve
    
def extract_genedata(datatype,geneinfo):
    """
    Gets the value of a specific type of data using regular expressions.
    """
    import re
    
    match = re.search(r"%s=\"([^\"]*)\"" % datatype,geneinfo)
    if match and match.groups():
        return match.groups()
    return None 

def str2gene(string):
    """
    Converts a block of gene data from a gbs file and converts it to a Gene object
    """
    string = string.strip()
    if not string:
        return None
    tokens = string.split("/")
    retval = Gene()
    if len(tokens[0]):
        retval.loc = tokens[0]
    if "complement" in retval.loc:
        retval.direction = "-"
    else:
        retval.direction = "+"
    
    name = extract_genedata("gene",string)
    if name:
        retval.name = name[0]
        
    note = extract_genedata("note",string)
    if note:
        retval.note = note[0]
        
    synonym = extract_genedata("gene_synonym",string)
    if synonym:
        retval.synonym = synonym[0]
        
    return retval

class Assembly():
    """
    Iterable abstraction of a GBS file for a chromosome assembly.
    
    When iterating, the next() function returns the next line of the gbs file
    as a tuple containing the key, which is the first 12 characters, and a value,
    which is the remaining characters.
    """
    def __init__(self,filename = None):
        self.data = {}
        self.filename = filename
        self.file = None
        
        if filename:
            self.data = load_keys(filename,DEFAULT_ASSEMBLY_KEYS)
            
    def get_gene(self,gene_name):
        """
        Finds a gene in the assembly using the gene's name.
        """
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
    
    def __iter__(self):
        return self
                    
    def open(self):
        """
        Opens the file associated with the assembly.
        """
        self.file = open(self.filename,"r")
       
    def opened(self):
        """
        Determines if the assembly file is open.
        """
        return self.file and not self.file.closed
    
    def closed(self):
        """
        Determines if the assembly file is closed.
        """
        return not self.opened()
    
        
    def next(self):
        """
        Reads the next line in the file. Returns they key-value pair for the LINE.
        """
        if self.closed():
            raise StopIteration
        line = self.file.readline()
        if not line:
            raise StopIteration
        if line[-1] == '\n':
            line = line[0:-1]
        key = line[0:GBS_KEY_LEN].strip().lower()
        val = line[GBS_KEY_LEN:].strip()
        return (key, val)