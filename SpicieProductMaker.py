
import json
import time
from pymongo import MongoClient
import csv

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["nucleotide_data"]

# filter = {"Features.product": {"$exists": True}}
# projection= {"_id":0, "Features.organism": 1, "Features.product": 1}
# dataResult = collection_data.find(filter, projection)

# list_res = list(dataResult)
# flat_data= []
# for diz in list_res:
#     organism=diz["Features"][0]["organism"]
#     indexProduct=1
#     while (diz["Features"][indexProduct]=={}):
#         indexProduct+=1
#     product=diz["Features"][indexProduct]["product"]
#     tupla=(organism, product)
#     flat_data.append(tupla)
#     continue

# with open('SpecieProduct.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)

#     for row in flat_data:
#         writer.writerow(row)

filter = {"GBSeq_feature-table.GBFeature_key":"CDS","GBSeq_feature-table.GBFeature_quals.GBQualifier_name":"product"}
proj = {"GBSeq_feature-table":1,"GBSeq_organism":1}

dataResult = collection_data.find(filter, proj)
def jsonList():
    tot_data = {}
    for dataN in dataResult:
        print(dataN["GBSeq_organism"])
        if dataN["GBSeq_organism"] not in tot_data:
            tot_data[dataN["GBSeq_organism"]] = []
        for data in dataN["GBSeq_feature-table"]:
            if data["GBFeature_key"] == "CDS" or data["GBFeature_key"] == "rRNA":
                for prod in data["GBFeature_quals"]:
                    if prod["GBQualifier_name"] == "product" and (data["GBFeature_key"],prod["GBQualifier_value"]) not in tot_data[dataN["GBSeq_organism"]]:
                        tot_data[dataN["GBSeq_organism"]].append((data["GBFeature_key"],prod["GBQualifier_value"]))
        print(len(tot_data[dataN["GBSeq_organism"]]))

    with open('Species.json', "w") as jswr:
        json.dump(tot_data,jswr,indent=4)



def csvWrite(dataResult):
    tot_data = []
    num = []
    for dataN in dataResult:
        print(dataN["GBSeq_organism"])
        for data in dataN["GBSeq_feature-table"]:
            if data["GBFeature_key"] == "CDS" or data["GBFeature_key"] == "rRNA":
                for prod in data["GBFeature_quals"]:
                    if prod["GBQualifier_name"] == "product":
                        number = 1
                        if [prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()] in tot_data:
                            iii = tot_data.index([prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()])
                            num[iii] = num[iii]+1
                        else:
                            tot_data.append([prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()])
                            num.append(number)
        print(len(tot_data))
    tot_tot = []
    for i in range(len(tot_data)):
        tot_tot.append(tot_data[i]+[num[i]])
    with open('SpecieProduct.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile,delimiter="|")
        writer.writerows(tot_tot)

#csvWrite(dataResult)

datasTaxon = [
    (2961,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2961",True,True,False),
    (145388,"https://www.ncbi.nlm.nih.gov/genome/?term=txid145388",True,True,True),
    (72520,"https://www.ncbi.nlm.nih.gov/genome/?term=txid72520",True,True,True),
    (70448,"https://www.ncbi.nlm.nih.gov/genome/?term=txid70448",True,True,True),
    (44056,"https://www.ncbi.nlm.nih.gov/genome/?term=txid44056",True,True,True),
    (41880,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41880",True,True,False),
    (41875,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41875",True,True,True),
    (36894,"https://www.ncbi.nlm.nih.gov/genome/?term=txid36894",True,True,False),
    (35677,"https://www.ncbi.nlm.nih.gov/genome/?term=txid35677",True,True,True),
    (3077,"https://www.ncbi.nlm.nih.gov/genome/?term=txid3077",True,True,True),
    (1764295,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1764295",True,True,True),
    (1650286,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1650286",True,True,True),
    (1486918,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1486918",True,True,False),
    (1093141,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1093141",True,True,True),
    (574566,"https://www.ncbi.nlm.nih.gov/genome/?term=txid574566",True,True,True),
    (554065,"https://www.ncbi.nlm.nih.gov/genome/?term=txid554065",True,True,True),
    (426638,"https://www.ncbi.nlm.nih.gov/genome/?term=txid426638",True,True,True),
    (265525,"https://www.ncbi.nlm.nih.gov/genome/?term=txid265525",True,True,False),
    (257627,"https://www.ncbi.nlm.nih.gov/genome/?term=txid257627",True,True,False),
    (2009235,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2009235",True,True,False),
    (2730355,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2730355",True,True,False),
]

# new_collection = db["table_basic"]

# for data in datasTaxon:
#     dataPush = {
#         "Link":data[1],
#         "GBFF":data[2],
#         "FNA":data[3],
#         "GFF":data[4]
#     }
#     new_collection.update_one({"TaxId":str(data[0])},{"$push":{"Genomes":dataPush}})

new_new_collection = db["markdown"]

mrk = open("README.md").read()

new_new_collection.update_one({"Title":"Markdown v.0.0.1"},{"$set":{"Text":mrk}})