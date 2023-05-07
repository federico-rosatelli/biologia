import csv
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
    print(record)
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

# taxon = ncbiSearchTaxon("chlorella[next level]")
# with open("taxonMeta.json","w") as jsw:
#     json.dump(taxon,jsw,indent=4)
#print(taxon)

def finderTaxon(fn):
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
                for tax in taxon:
                    if not collection_data.find_one({"TaxId":tax["TaxId"]}):
                        collection_data.insert_one(tax)
            except Exception as e:
                print(e)
                pass
            for t in taxon:
                nameT = t["ScientificName"]
                print(nameT,'taxon\t\t')
                try:
                    newTaxon = ncbiSearchTaxon(f"{nameT}[next level]")
                    for ntax in newTaxon:
                        if not collection_data.find_one({"TaxId":ntax["TaxId"]}):
                            collection_data.insert_one(ntax)
                except Exception as e:
                    print(e)
                    continue
            names.append(name)
    return

#taxons = finderTaxon('data/databaseCsv/microAlgaeDatabase.csv')
#db.taxonomy_data.deleteMany({"Lineage":{"$regex":"environmental samples"}})


# allBomb = collection_data.find({},{"ScientificName":1,"_id":0})

# listTuple = [[i["ScientificName"]] for i in allBomb]

# with open('OrganismList.csv', 'w') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerows(listTuple)


new_collection = db["taxonomy_tree"]

def algo():
    results = collection_data.find({})
    iter = 0
    for res in results:
        print(iter)
        iter += 1
        # if new_collection.find_one({"TaxId":res["ParentTaxId"]}):
        #     dataPush = {
        #             "TaxId":LineageEx["TaxId"],
        #             "Rank":LineageEx["Rank"],
        #             "ScientificName":LineageEx["ScientificName"],
        #             "SubClasses":[]
        #         }
        #     new_collection.insert_one(dataPush)
        for i in range(1, len(res["LineageEx"])):
            #for LineageEx in res["LineageEx"]:
            LineageEx = res["LineageEx"][i]
            if not new_collection.find_one({"TaxId":LineageEx["TaxId"]}):
                dataPush = {
                    "TaxId":LineageEx["TaxId"],
                    "Rank":LineageEx["Rank"],
                    "ScientificName":LineageEx["ScientificName"],
                    "SubClasses":[]
                }
                new_collection.insert_one(dataPush)
            if (i == 1):
                continue
            LineageExPrev = res["LineageEx"][i-1]

            if new_collection.find_one({"TaxId":LineageExPrev["TaxId"], "SubClasses":{"$elemMatch":{"TaxId":LineageEx["TaxId"]}}}):
                continue
            newDataPush = {
                    "TaxId":LineageEx["TaxId"],
                    "Rank":LineageEx["Rank"],
                    "ScientificName":LineageEx["ScientificName"],
                    }
            new_collection.update_one({"TaxId":LineageExPrev["TaxId"]},{"$push":{"SubClasses":newDataPush}})

        if not new_collection.find_one({"TaxId":res["ParentTaxId"], "SubClasses":{"$elemMatch":{"TaxId":res["TaxId"]}}}):
            newDataPush = {
                    "TaxId":res["TaxId"],
                    "Rank":res["Rank"],
                    "ScientificName":res["ScientificName"],
                    }
            new_collection.update_one({"TaxId":res["ParentTaxId"]},{"$push":{"SubClasses":newDataPush}})


#datas = new_collection.find_one({"ScientificName":"unclassified Chlorella"},{'_id':0})


def ncbiSearchNucleo(name:str) ->list:
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='nucleotide', term=name, rettype='gb', retmode='text', retmax=10000, api_key=key)
    record = Entrez.read(handle, validate=False)
    handle.close()
    print(f"Len of IDLIST:{len(record['IdList'])}")
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype='gb',retmode="xml",complexity=3, api_key=key)
    read = Entrez.read(handle)
    print(f"Len of EFETCH:{len(read)}")
    return read



taxon_collection = db["taxonomy_data"]
nucleo_collection = db["nucleotide_organism"]
# records = ncbiSearchNucleo("Scenedesmus bijugus")
# print(records[0])
csvOrganism = open('data/databaseCsv/microAlgaeDatabase.csv').readlines()
i = 0
for organism in csvOrganism:
    organism = organism.split(";")
    organism = organism[3]
    organism.strip()
    dataRank = taxon_collection.find_one({"ScientificName":organism})
    if not dataRank or dataRank["Rank"] != "species" or dataRank["Division"] == "Bacteria":
        continue
    if organism.split(" ")[1] == "sp." and len(organism.split(" ")) == 2:
        continue
    print(organism)
    if "_" in organism:
        organism = organism.split("_")[0]
    print(organism,f"{i} su {len(csvOrganism)}")
    i += 1
    try:
        records = ncbiSearchNucleo(organism)
    except Exception as e:
        print(e)
        continue
    for record in records:
        if not nucleo_collection.find_one({"GBSeq_locus":record["GBSeq_locus"]}):
            nucleo_collection.insert_one(record)


