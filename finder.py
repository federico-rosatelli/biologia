from pymongo import MongoClient
import re

CLUSTER = "localhost:27017"
def findTaxon(key):
    client = MongoClient('localhost', 27017)
    db = client["Biologia"]
    collection_taxon = db["taxonomy_data"]
    rgx = re.compile(f'*{key}*', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxon.find(info)
    for i in dataFind:
        print(i)

findTaxon("Eukaryota")
