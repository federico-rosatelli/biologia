
import json
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
        writer = csv.writer(csvfile)
        writer.writerows(tot_tot)

csvWrite(dataResult)