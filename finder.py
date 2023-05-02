import re
from pymongo import MongoClient



# Questa classe comprende al suo interno i metodi di interrogazione del database
# attraverso l'interfaccia accessibile via web all'indirizzo localhost:3000.
# Per permetterne il funzionamento, é necessario prima attivare sia il
# servizio di mongodb che istanziare il server.js come illustrato in genbank.py. 



CLUSTER = "localhost:27017"


def findTaxon(key):
    client = MongoClient('localhost', 27017)
    db = client["Biologia"]
    collection_taxon = db["taxonomy_data"]
    rgx = re.compile(f'{key}', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxon.find(info)
    for i in dataFind:
        print(i)


#PER ORA SONO SOLO COPIE DI findTaxon

def findGenome(key):
    client = MongoClient('localhost', 27017)
    db = client["Biologia"]
    collection_taxon = db["taxonomy_data"]
    rgx = re.compile(f'{key}', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxon.find(info)
    for i in dataFind:
        print(i)

def findProtein(key):
    client = MongoClient('localhost', 27017)
    db = client["Biologia"]
    collection_taxon = db["taxonomy_data"]
    rgx = re.compile(f'{key}', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxon.find(info)
    for i in dataFind:
        print(i)

def findProduct(key):
    client = MongoClient('localhost', 27017)
    db = client["Biologia"]
    collection_taxon = db["taxonomy_data"]
    rgx = re.compile(f'{key}', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxon.find(info)
    for i in dataFind:
        print(i)

findTaxon("Eukaryota")
