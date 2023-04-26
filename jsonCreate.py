import json
# fileJ = "data/sourceJson/microAlgaeSource.json"

# with open(fileJ) as jsr:
#     data:list[str] = json.load(jsr)


# dataFinal = {}

# for d in data:
#     if "sp." in d:
#         splt = d.split("sp.")
#         if not splt[0] in dataFinal:
#             dataFinal[splt[0]] = []
#         dataFinal[splt[0]].append(splt[1])
#     else:
#         dataFinal[d] = []


# with open("uniqueGeneMicroAlgae.json","w") as jsw:
#     json.dump(dataFinal,jsw,indent=4)

from pymongo import MongoClient

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["protein_data"]

filter = {}
dataResult = collection_data.find(filter,{"Id":1})

for i in dataResult:
    print(i)


