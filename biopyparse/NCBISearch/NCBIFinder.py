from Bio import SeqIO, Entrez

def ncbiSearchTaxon(name:str) -> list[dict] | None:
    '''Function that submit queries with given TaxonomyIDs to NCBI platform, section Taxonomy.
    Return data read'''
    handle = Entrez.esearch(db='taxonomy', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml")
    read:list[dict] = Entrez.read(handle)
    return read


def ncbiSearchNucleo(name:str) -> list[dict] | None:
    '''Function that submit queries with given a ScientificName to NCBI platform, section Nucleotide.
    Return data read'''
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='nucleotide', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    return read

def ncbiSearchProtein(name:str) -> list[dict] | None:
    '''Function that submit queries with given a ScientificName to NCBI platform, section Nucleotide.
    Return data read'''
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='protein', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="protein", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    return read