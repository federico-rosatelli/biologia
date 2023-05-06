import json
from pymongo import MongoClient
from Bio import SeqIO, Entrez

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["taxonomy_data"]
key="cc030996838fc52dd1a2653fad76bf5fe408"
Entrez.email = 'russo.1864451@studenti.uniroma1.it'

def ncbiSearchTaxon(name:str) ->list:
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='taxonomy', term=name, rettype='gb', retmode='text', retmax=10000, api_key=key)
    record = Entrez.read(handle, validate=False)
    handle.close()
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    #print(record["IdList"])
    handle = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml", api_key=key)
    read = Entrez.read(handle)
    return read


def aglo(name,rank):
    filter = {}
    dataResult = collection_data.find(filter)
    return

taxon = ncbiSearchTaxon("chlorella[next level]")
with open("taxonMeta.json","w") as jsw:
    json.dump(taxon,jsw,indent=4)
#print(taxon)

def finder(fn):
    taxons = []
    names = []
    rd = open(fn).readlines()
    for i in range(len(rd)):
        data = rd[i].split(";")
        name = data[3].split(" ")[0]
        if "_" in name:
            name = name.split("_")[0]
        if name != "" and name not in names:
            print(name,'\t\t',i,len(rd))
            try:
                taxon = ncbiSearchTaxon(f"{name}[next level]")
            except Exception as e:
                print(e)
                pass
            for t in taxon:
                nameT = t["ScientificName"]
                try:
                    newTaxon = ncbiSearchTaxon(f"{nameT}[next level]")
                    taxons = taxons + newTaxon
                except Exception as e:
                    print(e)
                    continue
            taxons = taxons+taxon
            names.append(name)
    return taxons

taxons = finder('data/databaseCsv/microAlgaeDatabase.csv')
with open("taxonMetaFull.json","w") as jsw:
    json.dump(taxons,jsw,indent=4)

with open("taxonMetaFull.json") as jsrd:
    taxons = json.load(jsrd)
print(len(taxons))
collection_data.insert_many(taxons,ordered=True)