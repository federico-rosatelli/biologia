########################################################################################
# Libraries & API
########################################################################################

from pymongo import MongoClient
from Bio import SeqIO, Entrez
from time import ctime, perf_counter
from Bio.Seq import Seq
from matplotlib import pyplot as plt
import urllib.request as download
import pandas as pd
import csv
import json
import requests
import os
import shutil
import gzip
import re
import random
import json
import hashlib
import argparse
import sqlite3
import sys



########################################################################################
# Init and Declaration of Global variables
########################################################################################

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["taxonomy_data"]
new_collection = db["taxonomy_tree"]
ignore_names = ["environmental samples"]
Entrez.api_key = "cc030996838fc52dd1a2653fad76bf5fe408"
Entrez.email = ""
pattern = r"\"?([-a-zA-Z0-9.`?{}]+@\w+\.\w+)\"?"            # pattern passed for checking email format purpose



########################################################################################
# Utility classes and methods
########################################################################################

def test_email(your_pattern):
  '''Function to check for valid email address'''
  pattern = re.compile(your_pattern)
  emails = ["john@example.com", 
            "python-list@python.org", 
            "wha.t.`1an?ug{}ly@email.com"]                  # here is an example list of email to check it at the end
  for email in emails:
    if not re.match(pattern, email):
        print ("You failed to match %s" % (email))
    elif not your_pattern:
        print("Forgot to enter a pattern!")
    else:
        print("Pass")


def read_kbd_input(inputQueue):
    '''Function to detect input text'''
    print('Ready for keyboard input: ')
    while (True):
        input_str = input()
        inputQueue.put(input_str)



########################################################################################

class bcolors:
    '''Error codes, intended to be used inside <class Error> and <PrintWarning>'''
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    WARN_BOX = WARNING + '[!] '
    OK_BOX = OKBLUE + '[*] '



########################################################################################

class  PrintWarning:
    '''Class in cooperation with bcolors class. Various code errors and warning with different colors
    can be printed to track events and logs'''


    def __init__(self, type:int, error:str="") -> None:
        '''Constructor that assign colors to number IDs'''
        self.error = error
        self.color = None
        if type == 0:
            self.color = bcolors.HEADER
        if type == 1:
            self.color = bcolors.OKBLUE
        if type == 2:
            self.color = bcolors.OKCYAN
        if type == 3:
            self.color = bcolors.OKGREEN
        if type == 4:
            self.color = bcolors.WARNING
        if type == 5:
            self.color = bcolors.FAIL
        if type == 6:
            self.color = bcolors.ENDC
        if type == 7:
            self.color = bcolors.BOLD
        if type == 8:
            self.color = bcolors.UNDERLINE
        if type == 9:
            self.color = bcolors.WARN_BOX
        if type == 10:
            self.color = bcolors.OK_BOX


    def __str__(self) -> str:
        '''Internal use for print errors'''
        return self.error

    
    def stdout(self, *strings:any) -> None:
        '''Custom print methods of error codes and timestamp to advise developers'''
        self.error = ' '.join(str(i) for i in strings)
        plus = ""
        if self.error[:1] == "\n":
            plus = "\n"
            self.error = self.error[1:]
        print(bcolors.BOLD + f"{plus}[{ctime()}] " + bcolors.ENDC + self.color + self.error +  bcolors.ENDC)



#################################################################################
# Output files methods (mainly .cs)
#################################################################################

def productsTable():
    '''Function to produce SpecieProduct.csv in lists/ folder'''
    products = open("SpecieProduct.csv").readlines()
    complete_collection = db["table_complete"]
    basic_collection = db["table_basic"]
    for i in range(len(products)):
        splt = products[i].split("|")
        product,name,qua = splt[0],splt[1].capitalize(),splt[2].strip()
        print(i,name)
        inserBasic = {
            "ProductName":product
        }
        insertComplete = {
            "ProductName":product,
            "QtyProduct":qua
        }
        basic_collection.update_one({"ScientificName":name},{"$push":{"Products":inserBasic}})
        complete_collection.update_one({"ScientificName":name},{"$push":{"Products":insertComplete}})



#################################################################################
# Query methods
#################################################################################

def ncbiSearchTaxon(name:str) ->list:
    '''TESTING'''
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
    '''TESTING'''
    filter = {}
    dataResult = collection_data.find(filter)
    return


def finderTaxonFromFile(fn):
    '''TESTING'''
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


def finderTaxon(name):
    '''TESTING'''
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


def unwrappingSpecies():
    '''TESTING'''
    finderTaxonFromFile('data/databaseCsv/microAlgaeDatabase.csv')
    taxons = finderTaxon('data/databaseCsv/microAlgaeDatabase.csv')
    db.taxonomy_data.deleteMany({"Lineage":{"$regex":"environmental samples"}})
    allBomb = collection_data.find({},{"ScientificName":1,"_id":0})
    listTuple = [[i["ScientificName"]] for i in allBomb]
    with open('OrganismList.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(listTuple)
    return


def algo():
    '''TESTING'''
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


def ncbiSearchNucleo(name:str) ->list:
    '''TESTING'''
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
    '''TESTING'''
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


def efetchTaxon(id):
    '''TESTING'''
    handle = Entrez.efetch(db="taxonomy", id=id, rettype='gb',retmode="xml")
    read = Entrez.read(handle)
    return read


def genusList():
    '''TESTING'''
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


def testNucleo():
    '''TESTING'''
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
    nucleos = nucleo_collection.count_documents(query,filter)
    print(nucleos)
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
        '''TESTING'''
        handle = Entrez.efetch(db="protein", id=ids, rettype='gb',retmode="xml")
        read = Entrez.read(handle)
        return read

    results = efetchProtein(proteins)
    nucleo_collection.insert_many(results)
    return


def efetchGenome(taxId):
    '''TESTING'''
    res = requests.get(f"https://ncbi.nlm.nih.gov/genome/?term={taxId}").text
    tot = []
    try:
        findIndex1 = res.index("_genomic.gff.gz")
        lastHttp1 = res[:findIndex1]
        firstHttpIndex1 = lastHttp1.rfind("https://ftp.ncbi.nlm.nih.gov/genomes")
        tot.append(lastHttp1[firstHttpIndex1:]+"_genomic.gff.gz")
    except Exception as e:
        print(f"GFF {e}")

    try:
        findIndex2 = res.index("_genomic.gbff.gz")
        lastHttp2 = res[:findIndex2]
        firstHttpIndex2 = lastHttp2.rfind("https://ftp.ncbi.nlm.nih.gov/genomes")
        tot.append(lastHttp2[firstHttpIndex2:]+"_genomic.gbff.gz")
    except Exception as e:
        print(f"GBFF {e}")
    try:
        findIndex3 = res.index("_genomic.fna.gz")
        lastHttp3 = res[:findIndex3]
        firstHttpIndex3 = lastHttp3.rfind("https://ftp.ncbi.nlm.nih.gov/genomes")
        tot.append(lastHttp3[firstHttpIndex3:]+"_genomic.fna.gz")
    except Exception as e:
        print(f"FNA {e}")

    return tot


def GFF():
    '''TESTING'''
    # creation directory for downloading process
    directory = "genomes/"
    temp = "temp/"
    parent = "./data"
    path = os.path.join(parent, directory)
    try:
        os.mkdir(path)
        print("Directory '% s' created" % path)
    except OSError as e:
        print(e)

    path2 = os.path.join(parent, temp)
    try:
        os.mkdir(path2)
        print("Directory '% s' created" % path2)
    except OSError as e:
        print(e)
    # creation list taxid
    listId = []
    taxonomyD = db["taxonomy_data"]
    myquery = taxonomyD.find({},{"TaxId":1,"_id":0})
    for data in myquery:
        tax = f"txid{data['TaxId']}"
        # print(tax)
        listId.append(tax)
    # print(listId)
    for i in listId:
        print(f"\nSearching for GFF, GBFF and FNA of {i} on [Genome] NCBI platform")
        # dataLink.append((i,efetchGenome(i)))
        dst = path + i

        if os.path.exists(dst):
            print(f"WARNING: Path already exists")
            if os.path.exists(dst + "/" + i + ".fna.gz"):
                print(f"WARNING: File FNA already downloaded\n")
            if os.path.exists(dst + "/" + i + ".gff.gz"):
                print(f"WARNING: File GFF already downloaded\n")
            if os.path.exists(dst + "/" + i + ".gbff.gz"):
                print(f"WARNING: File GBFF already downloaded\n")
        if efetchGenome(i) == []:
            continue
        dataLink = efetchGenome(i)
        for l in dataLink:
            if "gff" in l:
                d = download.urlretrieve(l, f'./data/temp/{i}' + '.gff.gz')
            elif "gbff" in l:
                d = download.urlretrieve(l, f'./data/temp/{i}' + '.gbff.gz')
            elif "fna" in l:
                d = download.urlretrieve(l, f'./data/temp/{i}' + '.fna.gz')
            print(d)
        if os.path.exists(path + i):
            print(f"Directory on path {path + i} already exists")
        else:
            try:
                os.mkdir(path + i)
                print("Directory '% s' created" % dst)
            except OSError as e:
                print(e)
        # print(dst)
        try:
            shutil.move("./data/temp/" + i + ".fna.gz", dst)
            shutil.move("./data/temp/" + i + ".gbff.gz", dst)
            shutil.move("./data/temp/" + i + ".gff.gz", dst)
            print(f"File .gz copied to path {dst}")
            os.remove("./data/temp/" + i + ".gff.gz")
            os.remove("./data/temp/" + i + ".gbff.gz")
            os.remove("./data/temp/" + i + ".fna.gz")
            print(f"Deleted .gz file in path {path}")
        except OSError as e:
            print(e)
        print(f"Transfer process for {i} complete\n")
        # decompression
        # print(dst)
        for j in os.walk(dst + "/"):
            # print(j)
            out = ''
            if "gff" in j:
                out = open(dst + "/" + i + ".gff", "w")
                format = ".gff"
            elif "gbff" in j:
                out = open(dst + "/" + i +  ".gbff", "w")
                format = ".gbff"
            elif "fna" in j:
                out = open(dst + "/" + i + ".fna", "w")
                format = ".fna"
            with gzip.open(dst + "/" + i + format) as gz:
                out.write(gz.read().decode("utf-8"))
                out.close()
        # for file in dataLink:
        #     print(f'{file}')


def genomeFind(name:str):
    '''TESTING'''
    handle = Entrez.esearch(db='genome', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    # print(f"Len of IDLIST:{len(record['IdList'])}")
    # if len(record["IdList"]) == 0:
    #     raise Exception("List Empty")
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    json_object = json.loads(f"{read}")

    json_formatted_str = json.dumps(json_object, indent=2)
    # print(f"Len of EFETCH:{len(read)}")
    # return read
    print(json_formatted_str)


def ncbiSearchNucleo1(name:str) ->list:
    '''TESTING'''
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
    '''TESTING'''
    daora = False
    nucleo_collection = db["nucleotide_data"]
    all_data = collection_data.find({},{"TaxId":1,"_id":0})
    for data in all_data:
        if f"txid{data['TaxId']}" == "txid1411642":
            daora = True
        else:
            print("No",f"txid{data['TaxId']}")
        if daora:
            print(f"txid{data['TaxId']}")
            try:
                insertDatas = ncbiSearchNucleo1(f"txid{data['TaxId']}[Organism:exp]")
                for ins in insertDatas:
                    ins.pop("GBSeq_sequence",None)
                    nucleo_collection.insert_one(ins)
                #nucleo_collection.insert_many(insertDatas)
            except Exception as e:
                print(e)


def nucleoResult():
    '''TESTING'''
    nucleoImport()
    # "txid257627"
    # "txid257627"
    nucleo_collection = db["nucleotide_data"]
    tot = ncbiSearchNucleo1("txid257627[Organism:exp]")
    i = 0
    for t in tot:
        print(i)
        t.pop("GBSeq_sequence",None)
        nucleo_collection.insert_one(t)
        i += 1
    nucleo_collection = db["nucleotide_data"]
    dataAll = nucleo_collection.find({},{"GBSeq_locus":1,"_id":1})
    for data in dataAll:
        query = {
            "_id": { "$ne": data["_id"] },
            "GBSeq_locus": data["GBSeq_locus"]
        }
        print(nucleo_collection.find_one(query))
    data = db.nucleo_collection.aggregate([
        {"$group" : { "_id": "$GBSeq_locus", "count": { "$sum": 1 } } },
        {"$match": {"_id" :{ "$ne" : None } , "count" : {"$gt": 1} } }, 
        {"$project": {"name" : "$_id", "_id" : 0} }
    ])
    return


def parseToBasic() -> list:
    '''TESTING'''
    taxon_collection = db["taxonomy_data"]
    nucleo_collection = db["nucleotide_data"]
    tableBasic_collection = db["table_basic"]
    protein_collection = db["protein_data"]
    dataRank = taxon_collection.find({},{"ScientificName":1,"TaxId":1,"_id":0})
    control = 0
    for data in dataRank:
        id = data["TaxId"]
        print(control,data["ScientificName"])
        finder = {"GBSeq_feature-table.GBFeature_quals.GBQualifier_value":f"taxon:{id}"}
        totNucl = nucleo_collection.find(finder,{"GBSeq_locus":1, "_id":0})
        totNucl = list(totNucl)
        print("Tot Nucleotides",len(totNucl))
        finder = {"GBSeq_feature-table.GBFeature_quals.GBQualifier_value":f"taxon:{id}"}
        proteins = protein_collection.find(finder,{"GBSeq_locus":1,"_id":0})
        proteins = list(proteins)
        print("Tot Proteins",len(proteins))
        inserter = {
            "ScientificName":data["ScientificName"],
            "TaxId":data["TaxId"],
            "Nucleotides":totNucl,
            "Proteins":proteins
        }
        tableBasic_collection.insert_one(inserter)
        control += 1   
    
    # per ogni ScientificName in Nucleotide:
        # popola "nucleotide_basic" con
            # GBSeq_locus
    return


def fattoBene():
    '''TESTING'''
    tt = db["nucleotide_basic"]
    tt1 = db["nucleotide_basic1"]
    tot = tt.find({},{"_id":0})
    for t in tot:
        for k in t:
            print(k)
            var = {
                "ScientificName":k,
                "GBSeq_locus":t[k]
            }
            tt1.insert_one(var)

        
def ncbiSearchNucleo123(name:str) ->bool:
    '''TESTING'''
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='nucleotide', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    print(f"Len of IDLIST:{len(record['IdList'])} di {name}")

    if len(record['IdList']) == 0:
        return True
    # if len(record["IdList"]) == 0:
    #     raise Exception("List Empty")
    # handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    # read = Entrez.read(handle)
    # print(f"Len of EFETCH:{len(read)}")
    # return read
    return False


def ncbiProtein(name:str) ->list:
    '''TESTING'''
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db='protein', term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    print(f"Len of IDLIST:{len(record['IdList'])}")
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db="protein", id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    print(f"Len of EFETCH:{len(read)}")
    return read


def proteinFind():
    '''TESTING'''
    protein_collection = db["protein_data"]
    taxon_collection = db["taxonomy_data"]
    findTax = taxon_collection.find({},{"TaxId":1,"_id":0})
    for data in findTax:
        tax = f"txid{data['TaxId']}"
        print(tax)
        try:
            allData = ncbiProtein(tax)
            for prot in allData:
                prot.pop("GBSeq_sequence",None)
                protein_collection.insert_one(prot)
        except Exception as e:
            print(e)


def newCollectionBene():
    '''TESTING'''
    new_collection = db["table_complete1"]
    old_collection = db["table_basic"]
    nucleotide_collection = db["nucleotide_data"]
    protein_collection = db["protein_data"]
    findAll = old_collection.find({})
    for data in findAll:
        nucleotides = [t["GBSeq_locus"] for t in data["Nucleotides"]]
        proteins = [t["GBSeq_locus"] for t in data["Proteins"]]
        print(f"{data['ScientificName']}; N: {len(nucleotides)}; P:{len(proteins)}")
        allDataN = nucleotide_collection.find({"GBSeq_locus":{"$in":nucleotides}})
        allDataP = protein_collection.find({"GBSeq_locus":{"$in":proteins}})
        allDataN = list(allDataN)
        allDataP = list(allDataP)
        dataInsert = {
            "ScientificName":data["ScientificName"],
            "TaxId":data["TaxId"],
            "Nucleotides":allDataN,
            "Proteins":allDataP
        }
        try:
            new_collection.insert_one(dataInsert)
        except Exception:
            print("TROPPO GROSSO")
            dataInsert = {
                "ScientificName":data["ScientificName"],
                "TaxId":data["TaxId"],
                "Nucleotides":[],
                "Proteins":[]
            }
            new_collection.insert_one(dataInsert)
            # for n in allDataN:
            #     new_collection.update_one({"TaxId":data["TaxId"]},{"$push":{"Nucleotides":n}})
            # for n in allDataP:
            #     new_collection.update_one({"TaxId":data["TaxId"]},{"$push":{"Proteins":n}})








#DEPRECATED: finder.py

# # # # # # Questa classe comprende al suo interno i metodi di interrogazione del database
# # # # # # attraverso l'interfaccia accessibile via web all'indirizzo localhost:3000.
# # # # # # Per permetterne il funzionamento, é necessario prima attivare sia il
# # # # # # servizio di mongodb che istanziare il server.js come illustrato in genbank.py. 


# # # # # def findTaxon(key):
# # # # #     client = MongoClient('localhost', 27017)
# # # # #     db = client["Biologia"]
# # # # #     collection_taxon = db["taxonomy_data"]
# # # # #     rgx = re.compile(f'{key}', re.IGNORECASE)
# # # # #     info = {"Lineage":rgx}
# # # # #     dataFind = collection_taxon.find(info)
# # # # #     for i in dataFind:
# # # # #         print(i)


# # # # # PER ORA SONO SOLO COPIE DI findTaxon

# # # # # def findGenome(key):
# # # # #     client = MongoClient('localhost', 27017)
# # # # #     db = client["Biologia"]
# # # # #     collection_taxon = db["taxonomy_data"]
# # # # #     rgx = re.compile(f'{key}', re.IGNORECASE)
# # # # #     info = {"Lineage":rgx}
# # # # #     dataFind = collection_taxon.find(info)
# # # # #     for i in dataFind:
# # # # #         print(i)

# # # # # def findProtein(key):
# # # # #     client = MongoClient('localhost', 27017)
# # # # #     db = client["Biologia"]
# # # # #     collection_taxon = db["taxonomy_data"]
# # # # #     rgx = re.compile(f'{key}', re.IGNORECASE)
# # # # #     info = {"Lineage":rgx}
# # # # #     dataFind = collection_taxon.find(info)
# # # # #     for i in dataFind:
# # # # #         print(i)

# # # # # def findProduct(key):
# # # # #     client = MongoClient('localhost', 27017)
# # # # #     db = client["Biologia"]
# # # # #     collection_taxon = db["taxonomy_data"]
# # # # #     rgx = re.compile(f'{key}', re.IGNORECASE)
# # # # #     info = {"Lineage":rgx}
# # # # #     dataFind = collection_taxon.find(info)
# # # # #     for i in dataFind:
# # # # #         print(i)

# # # # # findTaxon("Eukaryota")






















#DEPRECATED: genbank.py


# # # # # #######################


# # # # # class Global:
# # # # #     '''Per l'accesso globale a costanti e risorse. I link qui preseti prelevano i dati direttamente
# # # # #     dalle piattaforme web che forniscono i database contenenti i dati genomici.'''


# # # # #     CLUSTER = "localhost:27017"
# # # # #     WEBSOURCE = {
# # # # #         "Algae":{
# # # # #             "Web":"https://shigen.nig.ac.jp/algae/download/downloadFile/Strain_e.txt",
# # # # #             "File":"data/databaseCsv/algaeDatabase.csv"
# # # # #         },
# # # # #         "MicroAlgae":{
# # # # #             "Web":"",
# # # # #             "File":"data/databaseCsv/microAlgaeDatabase.csv"
# # # # #         }
# # # # #     }
# # # # #     SQL = {
# # # # #         "Init":"init.sql",
# # # # #         "Store":"Biologia.db"
# # # # #     }
# # # # #     JSON = {
# # # # #         "Path":"data/sourceJson/",
# # # # #         "Type":"JSONFile"
# # # # #     }
# # # # #     ALIGNMENT = {
# # # # #         "Path":"data/alignments/",
# # # # #         "Type":"TXTFile"
# # # # #     }



# # # # # class Parsing(object):
# # # # #     '''Classe che si occupa di trattare i dati forniti da GenBank.
# # # # #     Il file su cui si opera si da essere per scontato in formato .gbk'''


# # # # #     def __init__(self, file_name:str=None) -> None:
# # # # #         '''Costruttore'''
# # # # #         if not file_name:
# # # # #             self.record = []
# # # # #         else:
# # # # #             self.record = SeqIO.parse(file_name, "genbank")


# # # # #     def parsing_gene(self, verbose:bool=False) -> list:
# # # # #         '''Metodo che ritorna due liste:
# # # # #         - In "gene_array" vi sono i dati analizzati: per ogni file.gbk dato in input ci sono X geni,
# # # # #           e per ogni gene avremo un dizionario dei campi di interesse, come il Name, l'ID, e così via;
# # # # #         - In "hex_array" vi sono i dati trattati riferiti alla Sequenza. Salviamo infatti sia la sequenza
# # # # #           di basi azotate "as is", sia la versione codificata con SHA256.'''
# # # # #         gene_array = []
# # # # #         hex_array = []
# # # # #         count = 0
# # # # #         last = -1
# # # # #         errors = 0
# # # # #         PrintWarning(8).stdout("\nSTART PARSING")               # messaggio informativo su console
# # # # #         print("~"*56+"\n")
# # # # #         t_last = perf_counter()
        
# # # # #         for gene in self.record:
# # # # #             count += 1
# # # # #             if int((count/86100)*100) != last and verbose:              # Gestito in base al verbose, attenzione alla lunghezza del file (...)
# # # # #                 t_now = perf_counter()
# # # # #                 print("[%-50s] %d%% %d" % ('='*((last+1)//2),last+1,(t_now-t_last)*(86100//count)),end="\r") #??????
# # # # #                 last +=1
# # # # #             #print(count)
# # # # #             json_gene, json_convert, errors = self.parseGene(gene)
# # # # #             hex_array.append(json_convert)
# # # # #             gene_array.append(json_gene)
# # # # #         PrintWarning(3).stdout("\nParsing completed\n") 
# # # # #         PrintWarning(5).stdout("Process completed with", errors, 'errors\n')
# # # # #         return gene_array, hex_array


# # # # #     def parseGene(self, gene) -> tuple[dict, dict, int]:
# # # # #         '''Metodo che ritorna una tupla strutturata nel modo seguente:
# # # # #         - dict json_gene        che contiene Nome e ID del gene analizzato;
# # # # #         - dict json_convert     che contiene sia la sequenza "raw" che tradotta in SHA256;
# # # # #         - int  error            che segnala quanti errori si sono verificati durante la ricerca, ovvero quante volte il dato è risultato mancante       
# # # # #         Questa tupla viene utilizzata all'interno del metodo parsing_gene() per popolare i due array gee_array e hex_array.'''
# # # # #         json_gene = {}
# # # # #         json_convert = {}
# # # # #         json_gene["Name"] = gene.name
# # # # #         json_gene["Id"] = gene.id
# # # # #         check = False
# # # # #         errors = 0
# # # # #         try:
# # # # #             if gene.seq != "":
# # # # #                 #print(count)
# # # # #                 check = True
# # # # #                 sh = self.sha_index(str(gene.seq))
# # # # #                 json_gene["Seq_Hex"] = sh["hex"]        #str(gene.seq)-> SHA256 funzione
# # # # #                 json_convert["Seq_Hex"] = sh["hex"]     #str(gene.seq)-> SHA256 funzione
# # # # #                 json_convert["Seq_Raw"] = sh["seq"]     #str(gene.seq)-> SHA256 funzione
# # # # #         except Exception as e:
# # # # #             if check:
# # # # #                 print(e,"Entra qua nella sequenza:",gene.id)
# # # # #             errors += 1
# # # # #             pass
# # # # #         json_gene["Description"] = gene.description
# # # # #         json_gene["Features"] = []
# # # # #         for feature in gene.features:
# # # # #             feature_type = {"Type":feature.type}
# # # # #             #json_gene["Features"][feature.type] = {}
# # # # #             for qual in feature.qualifiers:
# # # # #                 feature_type[qual] = feature.qualifiers.get(qual)[0]
# # # # #             feature_type["Location"] = (int(feature.location.start),int(feature.location.end))
# # # # #             json_gene["Features"].append(feature_type)
# # # # #         return json_gene, json_convert, errors
    

# # # # #     def sha_index(self, seq:str) -> dict:
# # # # #         '''Metodo che, dato in input una squenza, ritorna in un dizionario la sequenza e il suo codice hash.
# # # # #         ACGTU sono le possibili unità di una qualsiasi sequenza. '''
# # # # #         return {
# # # # #             'hex':hashlib.sha256(seq.encode('utf-8')).hexdigest(),
# # # # #             'seq':seq
# # # # #         }
    

# # # # #     def save_data(self, json_array:tuple) -> None:
# # # # #         '''Crea un file json contentene i dati parsati.'''
# # # # #         with open(f"{Global.JSON['Path']}datastruct.json","w") as js:
# # # # #             json.dump({
# # # # #                         "struct":json_array[0],
# # # # #                         "hex":json_array[1]
# # # # #                         },js,indent=4)


# # # # # #######################


# # # # # class Database:
# # # # #     '''In questa classe vengono gestiti i diversi tipi di DB su cui i dati
# # # # #     possono essere salvati. Compito di questa classe é l'inizializzazione
# # # # #     e la manipolazione dei dati, insieme alla possibilità di inserimento
# # # # #     delle credenziali per dialogare con i server delle maggiori piattaforme web.'''


# # # # #     def __init__(self, verbose=False, type:str="", email:str=None) -> None:
# # # # #         '''Costruttore'''
# # # # #         self.file = Global.SQL["Store"]
# # # # #         ip = Global.CLUSTER.split(":")[0]
# # # # #         port = int(Global.CLUSTER.split(":")[1])
# # # # #         self.client = MongoClient(f'{ip}', port)
# # # # #         self.connection = None
# # # # #         self.mongo_collections = (type+"_data",type+"_hex")
# # # # #         self.printWarning = PrintWarning(10)
# # # # #         self.dataSource = {}
# # # # #         self.totLength = 0
# # # # #         self.diffAlg = {}
# # # # #         self.verbose = verbose
# # # # #         self.parse = Parsing()
# # # # #         self.email = email
# # # # #         print(email)
# # # # #         if email:
# # # # #             Entrez.email = email
    
    
# # # # #     def listSource(self)->list:
# # # # #         '''Restituisce una lista delle fonti da cui vengono prelevati i dati.'''
# # # # #         return [key for key in self.dataSource]
    

# # # # #     def addEmail(self, email) -> None:
# # # # #         '''Aggiorna il campo email per identificazione.'''
# # # # #         self.email = email
# # # # #         Entrez.email = email


# # # # #     def save_on_mongo(self, tuple_of_array:tuple) -> None:
# # # # #         '''Questo metodo genera il database in MongoDB (v. Requisiti all'inizio del file o il README.md).
# # # # #         I dati vengono salvati in un db locale di nome Biologia, dove poi verranno eseguite tutte le operazioni.'''
# # # # #         PrintWarning(8).stdout("\nStart saving on database")
# # # # #         db = self.client["Biologia"]  # attenzione a quando si richiama il DB da riga di comando: case sensitive
# # # # #         print(self.mongo_collections)
# # # # #         collection_data = db[self.mongo_collections[0]]
# # # # #         collection_convert = db[self.mongo_collections[1]]
# # # # #         count = 0
# # # # #         last = -1
# # # # #         for i in range(len(tuple_of_array[0])):
# # # # #             count+=1
# # # # #             if int((count/len(tuple_of_array[0]))*100) != last:
# # # # #                 print("[%-50s] %d%%" % ('='*((last+1)//2),last+1),end="\r")
# # # # #                 last +=1
# # # # #             filter = {"Name":tuple_of_array[0][i]["Name"]}
# # # # #             check_name = collection_data.find_one(filter)
# # # # #             if not check_name:
# # # # #                 collection_data.insert_one(tuple_of_array[0][i])
# # # # #             if tuple_of_array[1][i] != {}:
# # # # #                 filter = {"Seq_Hex":tuple_of_array[1][i]["Seq_Hex"]}
# # # # #                 check_hex = collection_convert.find_one(filter)
# # # # #                 if not check_hex:
# # # # #                     collection_convert.insert_one(tuple_of_array[1][i])
# # # # #         self.printWarning.stdout("Saved correctly on mongo")
    

# # # # #     def save_one_in_mongo(self,struct:dict,collection:str)->None:
# # # # #         '''Come il metodo save_on_mongo(), ma con sensibilità per singolo record.
# # # # #         In questo modo, possiamo aggiungere manualmente i dati di interesse.'''
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data_name = collection + "_data"
# # # # #         collection_hex_name = collection + "_hex"
# # # # #         collection_data = db[collection_data_name]
# # # # #         collection_hex = db[collection_hex_name]
# # # # #         data = struct["data"]
# # # # #         hex = struct["hex"]
# # # # #         collection_data.insert_one(data)
# # # # #         collection_hex.insert_one(hex)
# # # # #         self.printWarning.stdout("Saved correctly on mongo")

        
# # # # #     def save_on_sql(self,tuple_of_array:tuple, file=None):
# # # # #         '''Questo metodo genera il database in SQL (v. Requisiti all'inizio del file o il README.md).
# # # # #         Sebbene funzionante, al momento questo ramo del progetto è in STANDBY'''
# # # # #         if file != None:
# # # # #             self.file = file
# # # # #         try:
# # # # #             self.connection = sqlite3.connect(self.file)
# # # # #             cursor = self.connection.cursor()
# # # # #             if not self.file_exists():
# # # # #                 scr = open("init.sql").read()
# # # # #                 cursor.executescript(scr)
# # # # #             for i in range(len(tuple_of_array[0])):
# # # # #                 data_tuple = (tuple_of_array[0][i]["Name"],tuple_of_array[0][i]["Id"],tuple_of_array[0][i]["Seq_Hex"],tuple_of_array[0][i]["Description"],i)
# # # # #                 find = f"SELECT Seq_Hex FROM Biologia WHERE Seq_Hex='{tuple_of_array[0][i]['Seq_Hex']}'"
# # # # #                 res = cursor.execute(find)
# # # # #                 if not res:
# # # # #                     table ="""INSERT INTO Biologia VALUES(?,?,?,?,?)"""
# # # # #                     cursor.execute(table,data_tuple)
# # # # #         except Exception as e:
# # # # #             self.connection.close()
# # # # #             return PrintWarning(4).stdout(f"Error Database Biologia.db: {e}")
# # # # #         self.connection.commit()
# # # # #         if not self.connection:
# # # # #             self.connection.close()
# # # # #             return PrintWarning(5).stdout(f"Error Writing Database Biologia.db")
# # # # #         self.connection.close()
# # # # #         return None


# # # # #     def get_data_from_mongo(self, info:dict, saveOnJson:bool=False, algae:bool=False, micro_algae:bool=False) -> dict:
# # # # #         '''Metodo di tipo "Getter", permette la consultazione del DB basato su MongoDB'''
# # # # #         if self.client == None:
# # # # #             ip = CLUSTER.split(":")[0]
# # # # #             port = int(CLUSTER.split(":")[1])
# # # # #             self.client = MongoClient(f'{ip}', port)
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data = db[self.mongo_collections[0]]
# # # # #         collection_convert = db[self.mongo_collections[1]]
# # # # #         finder = collection_data.find(info)
# # # # #         dataSource = {}
# # # # #         for x in finder:
# # # # #             #print(x["Features"])
# # # # #             source = None
# # # # #             for f in x["Features"]:
# # # # #                 if f["Type"] == "source":
# # # # #                     source = f["organism"]
# # # # #                     if source not in dataSource:
# # # # #                         dataSource[source] = []
# # # # #             dataSource[source].append(x["Id"])
# # # # #         self.totLength = len(dataSource)
# # # # #         self.dataSource = dataSource
# # # # #         if saveOnJson:
# # # # #             self.save_on_json()
# # # # #         if algae:
# # # # #             totAlgae,same = self.isAlgae(dataSource)
# # # # #             self.diffAlg["Algae"] = same
# # # # #             algaeSource = {}
# # # # #             for algae in totAlgae:
# # # # #                 algaeSource[algae] = dataSource[algae]
# # # # #             self.dataSource = algaeSource
# # # # #             PrintWarning(2).stdout(f"Total algae {len(self.dataSource)} out of {self.totLength} with {len(same)} different gene")
# # # # #             if saveOnJson:
# # # # #                 self.save_on_json(fileName="algaeSource.json")
# # # # #         if micro_algae:
# # # # #             # ? self.datasource vs. datasource
# # # # #             totAlgae,same = self.isMicroAlgae(dataSource)
# # # # #             self.diffAlg["MicroAlgae"] = same
# # # # #             algaeSource = {}
# # # # #             for algae in totAlgae:
# # # # #                 algaeSource[algae] = dataSource[algae]
# # # # #             self.dataSource = algaeSource
# # # # #             PrintWarning(2).stdout(f"Total micro algae {len(self.dataSource)} out of {self.totLength} with {len(same)} different gene")
# # # # #             if saveOnJson:
# # # # #                 self.save_on_json(fileName="microAlgaeSource.json")
# # # # #         #
# # # # #         # self.printDifferenceAlgae()
# # # # #         # dataDiff = self.checkDifferenceAlgae()
# # # # #         #return self.dataSource
# # # # #         self.confronto()
# # # # #         return self.dataSource
    

# # # # #     def alignmentSeq(self, f1_key:str, f2_key:str) -> str:
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data = db["genetic_data"]
# # # # #         collection_convert = db["hex_to_seq"]
# # # # #         info = {"Features":{"$elemMatch":{"Type":"source","organism":f1_key}}}
# # # # #         # FIND PROTEIN
# # # # #         finder = collection_data.find_one(info)
# # # # #         if not finder:
# # # # #             PrintWarning(5).stdout(f"Error searching {f1_key}: Organism Not Found")
# # # # #             return None, None
# # # # #         if not 'Seq_Hex' in finder:
# # # # #             PrintWarning(5).stdout(f"Error searching {f1_key}: Hex Sequence Not Found")
# # # # #             return None, None
# # # # #         hex1 = finder['Seq_Hex']
# # # # #         info2 = {"Features":{"$elemMatch":{"Type":"source","organism":f2_key}}}
# # # # #         finder2 = collection_data.find_one(info2)
# # # # #         if not finder2:
# # # # #             PrintWarning(5).stdout(f"Error searching {f2_key}: Organism Not Found")
# # # # #             return None, None
# # # # #         if not 'Seq_Hex' in finder2:
# # # # #             PrintWarning(5).stdout(f"Error searching {f2_key}: Hex Sequence Not Found")
# # # # #             return None, None
# # # # #         hex2 = finder2['Seq_Hex']
# # # # #         info_hex = {"Seq_Hex":hex1}
# # # # #         finder_hex = collection_convert.find_one(info_hex)
# # # # #         if not finder_hex:
# # # # #             PrintWarning(5).stdout(f"Error searching {f1_key}: Hex Sequence Not Found in the Database")
# # # # #             return None, None
# # # # #         seq1 = finder_hex['Seq_Raw']
# # # # #         info_hex = {"Seq_Hex":hex2}
# # # # #         finder_hex = collection_convert.find_one(info_hex)
# # # # #         if not finder_hex:
# # # # #             PrintWarning(5).stdout(f"Error searching {f2_key}: Hex Sequence Not Found in the Database")
# # # # #             return None, None
# # # # #         seq2 = finder_hex['Seq_Raw']
# # # # #         smwt = smith_waterman.Alignment(seq1,seq2,show_table=self.verbose)
# # # # #         s1,s2 = smwt.localAlignment()
# # # # #         return s1,s2
    

# # # # #     def saveFileSeq(self, name:str, seq1:str, seq2:str) -> None:
# # # # #         seq1_fasta = ""
# # # # #         for i in range(len(seq1)):
# # # # #             if i%100 == 0 and i != 0:
# # # # #                 seq1_fasta += "\n"
# # # # #             elif i%10 == 0 and i != 0:
# # # # #                 seq1_fasta += " "
# # # # #             seq1_fasta += seq1[i]
# # # # #         seq2_fasta = ""
# # # # #         for i in range(len(seq2)):
# # # # #             if i%100 == 0 and i != 0:
# # # # #                 seq2_fasta += "\n"
# # # # #             elif i%10 == 0 and i != 0:
# # # # #                 seq2_fasta += " "
# # # # #             seq2_fasta += seq2[i]
# # # # #         with open(name,"w") as wr:
# # # # #             wr.write(seq1_fasta+"\n\n//\n\n"+seq2_fasta)
    

# # # # #     def alignment(self,fileIn,fileOut) -> None:
# # # # #         rd = open(fileIn).readlines()
# # # # #         for i in range(len(rd)):
# # # # #             for j in range(i+1,len(rd)):
# # # # #                 PrintWarning(3).stdout(f"Align {rd[i].strip()} with {rd[j].strip()}")
# # # # #                 seq1,seq2 = self.alignmentSeq(rd[i].strip(),rd[j].strip())
# # # # #                 if seq1 and seq2:   
# # # # #                     PrintWarning(3).stdout(f"Sequence length: {len(seq1)}\n")
# # # # #                     self.saveFileSeq(f"{Global.ALIGNMENT['Path']}{fileOut.split('.')[0]}_{rd[i].strip()}_{rd[j].strip()}.txt",seq1,seq2)


# # # # #     def printDifferenceAlgae(self) -> list:
# # # # #         '''Stampa le algae non contenute in microalgae'''
# # # # #         diff = [self.diffAlg[key] for key in self.diffAlg]
# # # # #         print(len(diff[0])-len(diff[1]))


# # # # #     def checkDifferenceAlgae(self) -> list:
# # # # #         '''Effettua il check delle algae non contenute in microalgae'''
# # # # #         diff = [self.diffAlg[key] for key in self.diffAlg]
# # # # #         tot = {}
# # # # #         for i in range(len(diff)):
# # # # #             for k in range(len(diff[i])):
# # # # #                 if diff[i][k] not in tot:
# # # # #                     tot[diff[i][k]] = 1
# # # # #                 else:
# # # # #                     tot[diff[i][k]] += 1  
# # # # #         algTot = []
# # # # #         for i in tot:
# # # # #             if tot[i] == 1:
# # # # #                 algTot.append(i)
# # # # #         return algTot

    
# # # # #     def save_on_json(self, fileName:str="dataSource.json") -> None:
# # # # #         '''Salva i dati su un file in formato .json'''
# # # # #         with open(f"{Global.JSON['Path']}{fileName}","w") as js:
# # # # #             json.dump(self.dataSource,js,indent=4)
    

# # # # #     def isAlgae(self,dataSource,rewrite:bool=False) -> tuple:
# # # # #         ''' Questo metodo fa il confronto con i dati presenti nel db
# # # # #         locale e tra i record dei csv in database/Csv, prelevati
# # # # #         dal Database europeo per motivi di accessibilità.
# # # # #         Datasource è inserito nel path data/sourcejson/struct.txt ed
# # # # #         è fornito come esempio.
# # # # #         Se trova una corrispondenza, ritorna le liste delle Algae
# # # # #         trovate e delle specie uniche (senza eventuali codici) '''
# # # # #         if not os.path.isfile(Global.WEBSOURCE["Algae"]["File"]) or rewrite:
# # # # #             r = requests.get(Global.WEBSOURCE["Algae"]["Web"])
# # # # #             with open(Global.WEBSOURCE["Algae"]["File"], "w") as wr:
# # # # #                 dataContent = r.text
# # # # #                 dataContent = dataContent.split("\n")
# # # # #                 writer = csv.writer(wr)
# # # # #                 for i in range(len(dataContent)):
# # # # #                     dataContent[i] = dataContent[i].split("\t")
# # # # #                     writer.writerow(dataContent[i])
# # # # #         file = open(Global.WEBSOURCE["Algae"]["File"],"r")
# # # # #         dataAlgae = csv.reader(file)
# # # # #         dataAlgae = list(dataAlgae)
# # # # #         totAlgae = []
# # # # #         unico = []
# # # # #         for key in dataSource:
# # # # #             keySplit = key.split(" ")
# # # # #             if keySplit[0] == "cf.":
# # # # #                 val_key = key.split(" ")[1]
# # # # #             else:
# # # # #                 val_key = key.split(" ")[0]

# # # # #             for algae in dataAlgae:
# # # # #                 if len(algae) > 1 and val_key.lower() in algae[2].lower():
# # # # #                     totAlgae.append(key)
# # # # #                     if val_key.lower() not in unico:
# # # # #                         unico.append(val_key.lower())
# # # # #         return totAlgae,unico


# # # # #     def isMicroAlgae(self,dataSource) -> tuple:
# # # # #         ''' Questo metodo fa il confronto con i dati presenti nel db
# # # # #         locale e tra i record dei csv in database/Csv, prelevati
# # # # #         dal Database europeo per motivi di accessibilità.
# # # # #         Datasource è inserito nel path data/sourcejson/struct.txt ed
# # # # #         è fornito come esempio.
# # # # #         Se trova una corrispondenza, ritorna le liste delle MicroAlgae
# # # # #         trovate e delle specie uniche (senza eventuali codici) '''
# # # # #         rl = open(Global.WEBSOURCE["MicroAlgae"]["File"]).readlines()
# # # # #         dataAlgae = []
# # # # #         for i in range(len(rl)):
# # # # #             data = rl[i].split(";")
# # # # #             dataAlgae.append(data[1:])
# # # # #         totAlgae = []
# # # # #         unico = []
# # # # #         for key in dataSource:
# # # # #             keySplit = key.split(" ")
# # # # #             if keySplit[0] == "cf.":
# # # # #                 val_key = key.split(" ")[1]
# # # # #             else:
# # # # #                 val_key = key.split(" ")[0]
# # # # #             for algae in dataAlgae:
# # # # #                 if len(algae) > 1 and val_key.lower() in algae[2].lower():
# # # # #                     totAlgae.append(key)
# # # # #                     if val_key.lower() not in unico:
# # # # #                         unico.append(val_key.lower())
# # # # #         return totAlgae,unico
    

# # # # #     def toFasta(self,name):
# # # # #         print(name)
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data = db["genetic_data"]
# # # # #         collection_convert = db["hex_to_seq"]
# # # # #         info = {"Features":{"$elemMatch":{"Type":"source","organism":name}}}
# # # # #         # FIND PROTEIN
# # # # #         finder = collection_data.find(info)
# # # # #         if not finder:
# # # # #             PrintWarning(5).stdout(f"Error searching {name}: Organism Not Found")
# # # # #             return None, None
# # # # #         cds = []


# # # # #     def file_exists(self)-> bool:
# # # # #         return os.path.isfile(self.file) and not os.stat(self.file).st_size == 0

    
# # # # #     def proteinFind(self, id) -> dict:
# # # # #         '''Matching con db esistente di Nucleotide'''
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data_nucleotide = db["nucleotide_data"]
# # # # #         info = {"Id":id,"Features":{"$elemMatch":{"Type":"CDS"}}}
# # # # #         finder_data = collection_data_nucleotide.find_one(info)
# # # # #         if not finder_data:
# # # # #             PrintWarning(5).stdout(f"Error searching {id}: CDS Not Found")
# # # # #             return None
# # # # #         collection_data = db["protein_data"]
# # # # #         collection_convert = db["protein_hex"]
# # # # #         protein_id = None
# # # # #         for f in finder_data["Features"]:
# # # # #             if f["Type"] == "CDS":
# # # # #                 if "protein_id" in f:
# # # # #                     protein_id = f["protein_id"]
# # # # #         if not protein_id:
# # # # #             PrintWarning(5).stdout(f"Error searching {id}: Protein Not Found")
# # # # #             return None
# # # # #         info = {"Id":protein_id}
# # # # #         finder_data = collection_data.find_one(info)
# # # # #         if not finder_data:
# # # # #             PrintWarning(5).stdout(f"Error searching {protein_id}: Protein Not Found In Database...","\n","\t\tSearching on NCBI...")
# # # # #             p1,p2 = self.ncbiSearch(protein_id,"protein")
# # # # #             dataFind = {
# # # # #                 'data':p1,
# # # # #                 'hex':p2
# # # # #             }
# # # # #             self.save_one_in_mongo(dataFind,"protein")
# # # # #             return dataFind
# # # # #         info = {"Seq_Hex":finder_data["Seq_Hex"]}
# # # # #         finder_hex = collection_convert.find_one(info)
# # # # #         if not finder_hex:
# # # # #             PrintWarning(5).stdout(f"Error searching {finder_data['Seq_Hex']}: Hex Not Found")
# # # # #             return None
# # # # #         dataFind = {
# # # # #             'data':finder_data,
# # # # #             'hex':finder_hex
# # # # #         }
# # # # #         return dataFind
    

# # # # #     def taxonFind(self,id) -> dict:
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data_nucleotide = db["nucleotide_data"]
# # # # #         info = {"Id":id,"Features":{"$elemMatch":{"Type":"source"}}}
# # # # #         finder_data = collection_data_nucleotide.find_one(info)
# # # # #         if not finder_data:
# # # # #             PrintWarning(5).stdout(f"Error searching {id}: Source Not Found")
# # # # #             return None
# # # # #         collection_data = db["taxonomy_data"]
# # # # #         collection_convert = db["taxonomy_hex"]
# # # # #         taxon_meta = None
# # # # #         for f in finder_data["Features"]:
# # # # #             if f["Type"] == "source":
# # # # #                 if "db_xref" in f:
# # # # #                     taxon_meta = f["db_xref"]
# # # # #         if not taxon_meta:
# # # # #             PrintWarning(5).stdout(f"Error searching {id}: db_xref Not Found")
# # # # #             return None
# # # # #         taxon_id = taxon_meta.split(":")[1]
# # # # #         info = {"TaxId":taxon_id}
# # # # #         finder_data = collection_data.find_one(info)
# # # # #         if not finder_data:
# # # # #             PrintWarning(5).stdout(f"Error searching {taxon_id}: Taxonomy Not Found In Database...","\n","\t\tSearching on NCBI...")
# # # # #             dataFind = self.ncbiSearchTaxon(taxon_id,"taxonomy")
# # # # #             for data in dataFind:
# # # # #                 collection_data.insert_one(data)
# # # # #             return dataFind
# # # # #         return finder_data
    

# # # # #     def genomeFind(self,id) -> dict:
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data_nucleotide = db["nucleotide_data"]
# # # # #         info = {"Features":{"$elemMatch":{"Type":"CDS","db_xref":{"$exists":True}}}}
# # # # #         finder_data = collection_data_nucleotide.find(info)
# # # # #         genes = []
# # # # #         for finder in finder_data:
# # # # #             gene_id = None
# # # # #             for f in finder["Features"]:
# # # # #                 if f["Type"] == "gene":
# # # # #                     if "gene" in f:
# # # # #                         gene_id = f["gene"]
# # # # #             if gene_id != None:
# # # # #                 print(gene_id)
# # # # #                 dataFind = self.ncbiSearchGenome(gene_id,"genome")
# # # # #                 print(dataFind)
# # # # #                 # info = {"TaxId":gene_id}
# # # # #                 # finder_data = collection_data.find_one(info)
# # # # #         return finder_data


# # # # #     def ncbiSearch(self,id:str,database:str) -> tuple:
# # # # #         # NB: in futuro la composizione del link potrebbe cambiare nella sua struttura, in base alla gestione interna di NCBI.
# # # # #         # TO-DO: se questo metodo non ritorna i dati correttamente andrà aggiornata la procedura di reperimento dei dati.
# # # # #         handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
# # # # #         record = SeqIO.read(handle, "genbank")
# # # # #         p1,p2,errors = self.parse.parseGene(record)
# # # # #         return p1,p2


# # # # #     def ncbiSearchTaxon(self,id:str,database:str) ->tuple:
# # # # #         try:
# # # # #             handle = Entrez.efetch(db=database, id=id, retmode="xml")
# # # # #             read = Entrez.read(handle)
# # # # #         except Exception as e:
# # # # #             print(e)
# # # # #             os.sleep(20)
# # # # #         return read
    

# # # # #     def ncbiSearchGenome(self,id:str,database:str) ->tuple:
# # # # #         handle = Entrez.efetch(db=database, id=id, retmode="xml")
# # # # #         read = Entrez.read(handle)
# # # # #         #print(read)
# # # # #         # handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
# # # # #         # record = SeqIO.read(handle, "genbank")
# # # # #         return read


# # # # #     def ritornodicose(self) -> None:
# # # # #         finder = {"Id":"","Features":{"$elemMatch":{"Type":"CDS"}}}
# # # # #         db = self.client["Biologia"]
# # # # #         collection_data_nucleotide = db["nucleotide_data"]
# # # # #         data = collection_data_nucleotide.count_documents(finder)
# # # # #         print(data)


# # # # #     def confronto(self): #CAMBIA NOME
# # # # #         for key in self.dataSource:
# # # # #             for id in self.dataSource[key]:
# # # # #                 data = self.proteinFind(id)
# # # # #                 if not data:
# # # # #                     PrintWarning(5).stdout(f"ID:{id}")
# # # # #                     #return None
# # # # #                 else:
# # # # #                     PrintWarning(3).stdout(f"Protein ID:{id}")
                
# # # # #                 data = self.taxonFind(id)
# # # # #                 if not data:
# # # # #                     PrintWarning(5).stdout(f"ID:{id}")
# # # # #                     #return None
# # # # #                 else:
# # # # #                     PrintWarning(3).stdout(f"Taxon ID:{id}")     
# # # # #         #data = self.genomeFind(4)
# # # # #                 # if not data:
# # # # #                 #     PrintWarning(5).stdout(f"ID:{id}")
# # # # #                 #     #return None
# # # # #                 # else:
# # # # #                 #     PrintWarning(3).stdout(f"Genome ID:{id}")



# # # # # #######################
# # # # # ######---MAIN---#######
# # # # # #######################


# # # # # def main(args:dict) -> None:
# # # # #     v = False
# # # # #     print("Welcome! This script is intended to be used for parsing data from a .gbk file.\nAfter parsing process is completed, you can check and query results in the DB selected as argument.")
# # # # #     if len(sys.argv) == 1:
# # # # #         print("Please add an argument at least\nDigit -h or --help after calling genbank.py to show possible input and try again\nRETURN:")
# # # # #     if args["verbose"]:
# # # # #         v = True
# # # # #     type = "nucleotide"
# # # # #     if args["protein"]:
# # # # #         type = "protein"
# # # # #     email = None
# # # # #     if args["email"]:
# # # # #         email = args["email"]
# # # # #     d = Database(verbose=v,type=type,email=email)
# # # # #     if args["nosql_mongo"] and args["sqlite3"]:
# # # # #         return PrintWarning(5).stdout("Can't select both mongo-db and sql for storing")
# # # # #     if args["file"]:
# # # # #         p = Parsing(args["file"])
# # # # #         data,conv = p.parsing_gene()
# # # # #         if args["json"]:
# # # # #             p.save_data((data,conv))
# # # # #         if args["nosql_mongo"]:
# # # # #             d.save_on_mongo((data,conv))
# # # # #         if args["sqlite3"]:
# # # # #             d.save_on_sql((data,conv))
# # # # #     elif args["find"]:
# # # # #         finder = {"Features":{"$elemMatch":{"Type":"source"}}}
# # # # #         save = False
# # # # #         algae = False
# # # # #         micro = False
# # # # #         if args["json"]:
# # # # #             save = True
# # # # #         if args["algae"]:
# # # # #             algae = True
# # # # #         if args["micro_algae"]:
# # # # #             micro = True
# # # # #         data = d.get_data_from_mongo(finder,saveOnJson=save,algae=algae,micro_algae=micro)
# # # # #         if args["list"]:
# # # # #             PrintWarning(4).stdout(d.listSource())
# # # # #         elif args["alignment"]:
# # # # #             fileIn,fileOut = args["alignment"]
# # # # #             d.alignment(fileIn=fileIn,fileOut=fileOut)
# # # # #         if args["fasta"]:
# # # # #             #print(args["fasta"])
# # # # #             d.toFasta(args["fasta"])
# # # # #     #d.isMicorAlgae()
# # # # #     return



# # # # # #Costrutto che permette di passare da riga di comando l'opzione desiderata per l'inserimento dei dati in un tipo di DB
# # # # # if __name__ == "__main__":
# # # # #     parser = argparse.ArgumentParser(
# # # # #                     prog='Biologia Database Parsing',
# # # # #                     description='What the program does',
# # # # #                     epilog='Text at the bottom of help')
# # # # #     parser.add_argument('-v', '--verbose',
# # # # #                         action='store_true')
# # # # #     parser.add_argument('-f', '--file',
# # # # #                         help='uses a file .gbk in input to retrieve data and store in local db.')
# # # # #     parser.add_argument('-n', '--nosql-mongo',
# # # # #                         action='store_true',
# # # # #                         help='uses MongoDB as database to store data. Do not select both mongo-db and sql for storing')
# # # # #     parser.add_argument('-s', '--sqlite3',
# # # # #                         action='store_true',
# # # # #                         help='uses SQLite3 as database to store data. Do not select both mongo-db and sql for storing')
# # # # #     parser.add_argument('-j', '--json',
# # # # #                         action='store_true',
# # # # #                         help='uses json as file to store data')
# # # # #     parser.add_argument('--find',
# # # # #                         action='store_true',
# # # # #                         help='to search some key field in data')
# # # # #     parser.add_argument('-a','--algae',
# # # # #                         action='store_true',
# # # # #                         help='if we are using algae branch on data')
# # # # #     parser.add_argument('-m','--micro-algae',
# # # # #                         action='store_true',
# # # # #                         help='if we are using micro-algae branch on data')
# # # # #     parser.add_argument('-l','--list',
# # # # #                         action='store_true',
# # # # #                         help='show data as list')
# # # # #     parser.add_argument('-p','--protein',
# # # # #                         action='store_true',
# # # # #                         help='change search archive from nucleotide to protein')
# # # # #     parser.add_argument('--alignment', nargs=2,
# # # # #                         metavar=('fromfile', 'tofile'),
# # # # #                         help='Align all organisms in __fromfile__ and save it in __tofile__',
# # # # #                         )
# # # # #     parser.add_argument('--email',
# # # # #                         help='for authentication purposes (NCBI)')
# # # # #     parser.add_argument('--fasta',
# # # # #                         help='ask for FASTA format')
# # # # #     args = vars(parser.parse_args())
# # # # #     main(args)
    

















# DEPRECATED: jsonCreate.py


# # # # # fileJ = "data/sourceJson/microAlgaeSource.json"

# # # # # with open(fileJ) as jsr:
# # # # #     data:list[str] = json.load(jsr)


# # # # # dataFinal = {}

# # # # # for d in data:
# # # # #     if "sp." in d:
# # # # #         splt = d.split("sp.")
# # # # #         if not splt[0] in dataFinal:
# # # # #             dataFinal[splt[0]] = []
# # # # #         dataFinal[splt[0]].append(splt[1])
# # # # #     else:
# # # # #         dataFinal[d] = []


# # # # # with open("uniqueGeneMicroAlgae.json","w") as jsw:
# # # # #     json.dump(dataFinal,jsw,indent=4)



# # # # CLUSTER = "localhost:27017"
# # # # client = MongoClient('localhost', 27017)
# # # # db = client["Biologia"]
# # # # collection_data = db["protein_data"]

# # # # filter = {}
# # # # dataResult = collection_data.find(filter,{"Id":1})
# # # # collection_data = db["genetic_data"]

# # # # list_res = []


# # # # for protein in dataResult:
# # # #     filter = {"Features":{"$elemMatch":{"Type":"CDS","protein_id":protein["Id"]}}}
# # # #     result = collection_data.find_one(filter)
# # # #     if result:
# # # #         list_res.append(protein["Id"])

# # # # print(list_res)



















# NOT-DEPRECATED: SpicieProductMaker.py


# # # # import json
# # # # from pymongo import MongoClient
# # # # import csv

# # # # CLUSTER = "localhost:27017"
# # # # client = MongoClient('localhost', 27017)
# # # # db = client["Biologia"]
# # # # collection_data = db["nucleotide_data"]

# # # # # filter = {"Features.product": {"$exists": True}}
# # # # # projection= {"_id":0, "Features.organism": 1, "Features.product": 1}
# # # # # dataResult = collection_data.find(filter, projection)

# # # # # list_res = list(dataResult)
# # # # # flat_data= []
# # # # # for diz in list_res:
# # # # #     organism=diz["Features"][0]["organism"]
# # # # #     indexProduct=1
# # # # #     while (diz["Features"][indexProduct]=={}):
# # # # #         indexProduct+=1
# # # # #     product=diz["Features"][indexProduct]["product"]
# # # # #     tupla=(organism, product)
# # # # #     flat_data.append(tupla)
# # # # #     continue

# # # # # with open('SpecieProduct.csv', 'w', newline='') as csvfile:
# # # # #     writer = csv.writer(csvfile)

# # # # #     for row in flat_data:
# # # # #         writer.writerow(row)

# # # # filter = {"GBSeq_feature-table.GBFeature_key":"CDS","GBSeq_feature-table.GBFeature_quals.GBQualifier_name":"product"}
# # # # proj = {"GBSeq_feature-table":1,"GBSeq_organism":1}

# # # # dataResult = collection_data.find(filter, proj)
# # # # def jsonList():
# # # #     tot_data = {}
# # # #     for dataN in dataResult:
# # # #         print(dataN["GBSeq_organism"])
# # # #         if dataN["GBSeq_organism"] not in tot_data:
# # # #             tot_data[dataN["GBSeq_organism"]] = []
# # # #         for data in dataN["GBSeq_feature-table"]:
# # # #             if data["GBFeature_key"] == "CDS" or data["GBFeature_key"] == "rRNA":
# # # #                 for prod in data["GBFeature_quals"]:
# # # #                     if prod["GBQualifier_name"] == "product" and (data["GBFeature_key"],prod["GBQualifier_value"]) not in tot_data[dataN["GBSeq_organism"]]:
# # # #                         tot_data[dataN["GBSeq_organism"]].append((data["GBFeature_key"],prod["GBQualifier_value"]))
# # # #         print(len(tot_data[dataN["GBSeq_organism"]]))

# # # #     with open('Species.json', "w") as jswr:
# # # #         json.dump(tot_data,jswr,indent=4)



# # # # def csvWrite(dataResult):
# # # #     tot_data = []
# # # #     num = []
# # # #     for dataN in dataResult:
# # # #         print(dataN["GBSeq_organism"])
# # # #         for data in dataN["GBSeq_feature-table"]:
# # # #             if data["GBFeature_key"] == "CDS" or data["GBFeature_key"] == "rRNA":
# # # #                 for prod in data["GBFeature_quals"]:
# # # #                     if prod["GBQualifier_name"] == "product":
# # # #                         number = 1
# # # #                         if [prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()] in tot_data:
# # # #                             iii = tot_data.index([prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()])
# # # #                             num[iii] = num[iii]+1
# # # #                         else:
# # # #                             tot_data.append([prod["GBQualifier_value"].lower(), dataN["GBSeq_organism"].lower()])
# # # #                             num.append(number)
# # # #         print(len(tot_data))
# # # #     tot_tot = []
# # # #     for i in range(len(tot_data)):
# # # #         tot_tot.append(tot_data[i]+[num[i]])
# # # #     with open('SpecieProduct.csv', 'w', newline='') as csvfile:
# # # #         writer = csv.writer(csvfile,delimiter="|")
# # # #         writer.writerows(tot_tot)

# # # # #csvWrite(dataResult)

# # # # datasTaxon = [
# # # #     (2961,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2961",True,True,False),
# # # #     (145388,"https://www.ncbi.nlm.nih.gov/genome/?term=txid145388",True,True,True),
# # # #     (72520,"https://www.ncbi.nlm.nih.gov/genome/?term=txid72520",True,True,True),
# # # #     (70448,"https://www.ncbi.nlm.nih.gov/genome/?term=txid70448",True,True,True),
# # # #     (44056,"https://www.ncbi.nlm.nih.gov/genome/?term=txid44056",True,True,True),
# # # #     (41880,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41880",True,True,False),
# # # #     (41875,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41875",True,True,True),
# # # #     (36894,"https://www.ncbi.nlm.nih.gov/genome/?term=txid36894",True,True,False),
# # # #     (35677,"https://www.ncbi.nlm.nih.gov/genome/?term=txid35677",True,True,True),
# # # #     (3077,"https://www.ncbi.nlm.nih.gov/genome/?term=txid3077",True,True,True),
# # # #     (1764295,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1764295",True,True,True),
# # # #     (1650286,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1650286",True,True,True),
# # # #     (1486918,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1486918",True,True,False),
# # # #     (1093141,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1093141",True,True,True),
# # # #     (574566,"https://www.ncbi.nlm.nih.gov/genome/?term=txid574566",True,True,True),
# # # #     (554065,"https://www.ncbi.nlm.nih.gov/genome/?term=txid554065",True,True,True),
# # # #     (426638,"https://www.ncbi.nlm.nih.gov/genome/?term=txid426638",True,True,True),
# # # #     (265525,"https://www.ncbi.nlm.nih.gov/genome/?term=txid265525",True,True,False),
# # # #     (257627,"https://www.ncbi.nlm.nih.gov/genome/?term=txid257627",True,True,False),
# # # #     (2009235,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2009235",True,True,False),
# # # #     (2730355,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2730355",True,True,False),
# # # # ]

# # # # new_collection = db["table_basic"]

# # # # for data in datasTaxon:
# # # #     dataPush = {
# # # #         "Link":data[1],
# # # #         "GBFF":data[2],
# # # #         "FNA":data[3],
# # # #         "GFF":data[4]
# # # #     }
# # # #     new_collection.update_one({"TaxId":str(data[0])},{"$push":{"Genomes":dataPush}})














# NOT-DEPRECATED: smith_waterman.py





# # # # def printTable(table,gene,trace=[]):
# # # #     #print([i for i in range(len(table[0]))])
# # # #     indx = "  "
# # # #     for i in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
# # # #         indx += gene[1][i]+" "*3
# # # #     print(indx)
# # # #     print('_'*10*4)
# # # #     for i in range(trace[len(trace)-1][0],trace[len(trace)-1][0]+10):
# # # #         print("| ",end="")
# # # #         for j in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
# # # #             #print(table[i][j])
# # # #             if (i,j) in trace:
# # # #                 print(bcolors.OKGREEN+str(table[i][j])+bcolors.ENDC+" "*(2-len(str(table[i][j]))),end="| ")
# # # #             else:
# # # #                 print(str(table[i][j])+" "*(2-len(str(table[i][j]))),end="| ")
# # # #         if i != 0:
# # # #             print(gene[0][i])
# # # #         else:
# # # #             print("")
# # # #     print('\n')


# # # # def saveTable(table, trace=[]):
# # # #     if len(table)<40 or len(table[0])<40:
# # # #         return
# # # #     tableCopy = []
# # # #     for i in range(10,40):
# # # #         tbCp = []
# # # #         for j in range(10,40):
# # # #             tbCp.append(table[i][j])
# # # #         tableCopy.append(tbCp)
# # # #     color = [["w" for j in range(len(tableCopy[0]))]for i in range(len(tableCopy))]
# # # #     for i in range(len(tableCopy)):
# # # #         for j in range(len(tableCopy[i])):
# # # #             if (i+10,j+10) in trace:
# # # #                 color[i][j] = "#56b5fd"
# # # #     # for i in range(10,20):
# # # #     #     color[trace[len(trace)-i-1][0]][trace[len(trace)-i-1][1]] = "#56b5fd"
    
# # # #     fig,ax = plt.subplots()
# # # #     fig.patch.set_visible(False)
# # # #     ax.axis('off')
# # # #     ax.axis('tight')
# # # #     df = pd.DataFrame(tableCopy)
# # # #     ax.table(cellText=df.values,cellColours=color,colLabels=df.columns,loc='center')
# # # #     #fig.tight_layout()
# # # #     plt.savefig('table.png',bbox_inches='tight')
# # # #     #plt.show()


# # # # class Alignment:


# # # #     def __init__(self, *seqs:str, gap:int=1, show_table:bool=False) -> None:
# # # #         if len(seqs)<2:
# # # #             print(bcolors.FAIL+f"MinMaxError: input len for sequences must be >= 2, got {len(seqs)}"+bcolors.ENDC)
# # # #             return
# # # #         self.seqs = [Seq(i) for i in seqs]
# # # #         self.seq1 = seqs[0]
# # # #         self.seq2 = seqs[1]
# # # #         self.gap = gap
# # # #         self.show_table = show_table
    

# # # #     def __str__(self) -> str:
# # # #         return '\n'.join(self.seqs)


# # # #     def __len__(self) -> int:
# # # #         x = 1
# # # #         for i in self.seqs:
# # # #             x *= len(i)
# # # #         return x
    

# # # #     def __call__(self, seq:str) -> list:
# # # #         self.seqs.append(seq)
# # # #         return self.seqs


# # # #     def createScoreMatrix(self, lnSeq1:int, lnSeq2:int) -> tuple:
# # # #         score_matrix = [[0 for _ in range(lnSeq2 + 1)] for _ in range(lnSeq1 + 1)]

# # # #         max_score = 0
# # # #         max_index = (0,0)
# # # #         for i in range(1, lnSeq1 + 1):
# # # #             for j in range(1, lnSeq2 + 1):
# # # #                 score = 1 if self.seq1[i-1] == self.seq2[j-1] else -1

# # # #                 score_matrix[i][j] = max(
# # # #                     0,
# # # #                     score_matrix[i-1][j-1] + score,
# # # #                     score_matrix[i-1][j] - self.gap,
# # # #                     score_matrix[i][j-1] - self.gap
# # # #                 )

# # # #                 if score_matrix[i][j] > max_score:
# # # #                     max_score = score_matrix[i][j]
# # # #                     max_index = (i, j)
# # # #         return score_matrix,max_index


# # # #     def localAlignment(self, save_table:bool=False) -> tuple:
# # # #         if len(self.seqs)>2:
# # # #             print(bcolors.WARN_BOX+f"Warning! Only 2 arguments were expected, but got {len(self.seqs)}.\n\t-The algorithm will use only the first 2 sequences..."+bcolors.ENDC)
# # # #         aligned_seq1 = ""
# # # #         aligned_seq2 = ""

# # # #         score_matrix, max_index_score_matrix = self.createScoreMatrix(len(self.seq1),len(self.seq2))

# # # #         i, j = max_index_score_matrix
# # # #         traceback = []
# # # #         print(f"Table length. First sequence:{len(self.seq1)}, Second sequence:{len(self.seq2)}")
# # # #         print(f"Max term: {score_matrix[i][j]}; in index: {max_index_score_matrix}")
# # # #         f1 = 0
# # # #         while score_matrix[i][j] != 0:
# # # #             traceback.append((i,j))
# # # #             if score_matrix[i-1][j-1] >= score_matrix[i-1][j] and score_matrix[i-1][j-1] >= score_matrix[i][j-1]:
# # # #                 aligned_seq1 = self.seq1[i-1] + aligned_seq1
# # # #                 aligned_seq2 = self.seq2[j-1] + aligned_seq2
# # # #                 i, j = i-1, j-1
# # # #                 f1 += 1
# # # #             elif score_matrix[i-1][j] >= score_matrix[i-1][j-1] and score_matrix[i-1][j] >= score_matrix[i][j-1]:
# # # #                 aligned_seq1 = self.seq1[i-1] + aligned_seq1
# # # #                 aligned_seq2 = '-' + aligned_seq2
# # # #                 i -= 1
# # # #             else:
# # # #                 aligned_seq1 = '-' + aligned_seq1
# # # #                 aligned_seq2 = self.seq2[j-1] + aligned_seq2
# # # #                 j -= 1
        
# # # #         if self.show_table:
# # # #             printTable(score_matrix,(self.seq1,self.seq2),trace=traceback)
# # # #             if save_table:
# # # #                 saveTable(score_matrix,trace=traceback)
# # # #         print(f"{(f1/len(aligned_seq1))*100}% of alignment")
# # # #         return aligned_seq1, aligned_seq2
    

# # # #     def globalAlignment(self):
# # # #         return 0


# # # # # a = Alignment("AGTCCCTGATTTAGTCCCTGATTTAGTATTTAGTCCCTGATTTAGTATTTAGTCCCTGATTTAGTCCCTGATTT","TTTAGTCCCTGATTTAGTTTTAGTCCCTGATTTAGTTTTAGTCCCTGATTTAGT",show_table=True)
# # # # # a.localAlignment(save_table=True)



















#######################
######---MAIN---#######
#######################

def main(args:dict) -> None:

    # genus = collection_data.find({"Rank":"genus","Division":{"$not":re.compile("Bacteria")}})
    # """SELECT * FROM COLLECTION WHERE RANK=genus AND NOT =bacteria"""
    #datas = new_collection.find_one({"ScientificName":"unclassified Chlorella"},{'_id':0})

    # taxon = ncbiSearchTaxon("chlorella[next level]")
    # with open("taxonMeta.json","w") as jsw:
    #     json.dump(taxon,jsw,indent=4)
    #print(taxon)

    nb = input('Please insert email')
    #unwrappingSpecies()
    #algo()
    #aglo()
    #testNucleo()
    #GFF()
    #nucleoResult()
    #genomeFind("Chlorella variabilis")
    #parseToBasic()
    #fattoBene() 
    #proteinFind()
    #newCollectionBene()
    #productsTable()
    return