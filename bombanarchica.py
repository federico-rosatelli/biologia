import csv
import json
from pymongo import MongoClient
from Bio import SeqIO, Entrez

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["taxonomy_data"]
Entrez.email = 'russo.1864451@studenti.uniroma1.it'
Entrez.api_key = "cc030996838fc52dd1a2653fad76bf5fe408"

def ncbiSearchTaxon(name:str) ->list:
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='taxonomy', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    # print(record)
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    #print(record["IdList"])
    handle = Entrez.efetch(db="taxonomy", id=record["IdList"], retmode="xml")
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

def finderTaxonFromFile(fn):
    names = []
    rd = open(fn).readlines()
    for i in range(len(rd)):
        data = rd[i].split(";")
        name = data[3].split(" ")[0]
        if "_" in name:
            name = name.split("_")[0]
        if name != "" and name not in names:
            finderTaxon(name)
        names.append(name)

ignore_names = ["environmental samples"]
def finderTaxon(name):
    if name in ignore_names:
        return
    try:
        taxon = ncbiSearchTaxon(f"{name}[next level]")
        for tax in taxon:
            if (not collection_data.find_one({"TaxId":tax["TaxId"]}) and tax["Rank"] == "species"):
                    print(tax["ScientificName"])
                    collection_data.insert_one(tax)
            finderTaxon(tax["ScientificName"])
    except Exception as e:
        return
    
#finderTaxonFromFile('data/databaseCsv/microAlgaeDatabase.csv')

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

#algo()
#datas = new_collection.find_one({"ScientificName":"unclassified Chlorella"},{'_id':0})


def ncbiSearchNucleo(name:str) ->list:
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='nucleotide', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    print(f"Len of IDLIST:{len(record['IdList'])}")
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    print(f"Len of EFETCH:{len(read)}")
    return read


def nucleoImport():
    taxon_collection = db["taxonomy_data"]
    nucleo_collection = db["nucleotide_organism"]
    records = ncbiSearchNucleo("Scenedesmus bijugus")
    print(records[0])
    csvOrganism = open('data/databaseCsv/microAlgaeDatabase.csv').readlines()

    listCompl = [csvOrganism[data].split(";")[3].strip() for data in range(0,200)]
    for i in range(200,len(csvOrganism)):
        organism = csvOrganism[i]
        
        organism = organism.split(";")
        organism = organism[3]
        organism.strip()
        dataRank = taxon_collection.find_one({"ScientificName":organism})
        if not dataRank or dataRank["Rank"] != "species" or dataRank["Division"] == "Bacteria":
            continue
        if organism.split(" ")[1] == "sp." and len(organism.split(" ")) == 2:
            continue
        if organism in listCompl:
            continue
        print(organism)
        if "_" in organism:
            organism = organism.split("_")[0]
        print(organism,f"{i} su {len(csvOrganism)}")
        
        listCompl.append(organism)
        try:
            records = ncbiSearchNucleo(organism)
        except Exception as e:
            print(e)
            continue
        for record in records:
            record.pop("GBSeq_sequence",None)
            if not nucleo_collection.find_one({"GBSeq_locus":record["GBSeq_locus"]}):
                nucleo_collection.insert_one(record)

# import re
# genus = collection_data.find({"Rank":"genus","Division":{"$not":re.compile("Bacteria")}})
# """SELECT * FROM COLLECTION WHERE RANK=genus AND NOT =bacteria"""

def efetchTaxon(id):
    handle = Entrez.efetch(db="taxonomy", id=id, rettype='gb',retmode="xml")
    read = Entrez.read(handle)
    return read

def genusList():
    import re
    new_collection = db["taxonomy_tree"]
    genus = new_collection.find({"Rank":"genus"})
    datas = [["Scientific Name","Taxon Id","Division"]]
    for gene in genus:
        print(gene["ScientificName"])
        taxons = efetchTaxon(gene["TaxId"])
        for taxon in taxons:
            if taxon["Division"] != "Bacteria":
                datas.append([taxon["ScientificName"],taxon["TaxId"],taxon["Division"]])

    with open('GenusList.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(datas)


nucleo_collection = db["nucleotide_organism"]
query = {"GBSeq_feature-table":{"$elemMatch":{"GBFeature_key":"CDS","GBFeature_quals":{"$elemMatch":{"GBQualifier_name":"protein_id"}}}}}
filter = {"GBSeq_feature-table":{"$elemMatch":{"GBFeature_key": "CDS"}}}
res = nucleo_collection.aggregate([
    {"$match": {'GBSeq_feature-table.GBFeature_key': 'CDS','GBSeq_feature-table.GBFeature_quals.GBQualifier_name':"protein_id"}},
    {"$project": {
        "GBSeq_feature-table": {"$filter": {
            "input": '$GBSeq_feature-table',
            "as": 'cds',
            "cond": {"$eq": ['$$cds.GBFeature_key', 'CDS']}
        }},
        "_id": 0,
    }}
])
#nucleos = nucleo_collection.count_documents(query,filter)
#print(nucleos)
proteins = []
count = 0
for nucleo in res:
    print(count)
    cdss = nucleo["GBSeq_feature-table"]
    for cds in cdss:
        for qual in cds["GBFeature_quals"]:
            if qual["GBQualifier_name"] == "protein_id":
                proteins.append(qual["GBQualifier_value"] in proteins)
    count += 1

def efetchProtein(ids):
    handle = Entrez.efetch(db="protein", id=ids, rettype='gb',retmode="xml")
    read = Entrez.read(handle)
    return read

results = efetchProtein(proteins)
nucleo_collection.insert_many(results)
