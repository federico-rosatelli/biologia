
from pymongo import MongoClient
import csv

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["nucleotide_data"]

filter = {"Features.product": {"$exists": True}}
projection= {"_id":0, "Features.organism": 1, "Features.product": 1}
dataResult = collection_data.find(filter, projection)

list_res = list(dataResult)
flat_data= []
for diz in list_res:
    organism=diz["Features"][0]["organism"]
    indexProduct=1
    while (diz["Features"][indexProduct]=={}):
        indexProduct+=1
    product=diz["Features"][indexProduct]["product"]
    tupla=(organism, product)
    flat_data.append(tupla)
    continue

with open('SpecieProduct.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    for row in flat_data:
        writer.writerow(row)