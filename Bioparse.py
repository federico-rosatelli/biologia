from pymongo import MongoClient
from Bio import SeqIO, Entrez
import csv
import json
import requests
import urllib.request as download
import os
import shutil
import gzip


# Dichiarazione e inizializzazione delle variabili globali #

CLUSTER = "localhost:27017"
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_data = db["taxonomy_data"]
ignore_names = ["environmental samples"]
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


def unwrappingSpecies():
    finderTaxonFromFile('data/databaseCsv/microAlgaeDatabase.csv')

    taxons = finderTaxon('data/databaseCsv/microAlgaeDatabase.csv')
    db.taxonomy_data.deleteMany({"Lineage":{"$regex":"environmental samples"}})


    allBomb = collection_data.find({},{"ScientificName":1,"_id":0})

    listTuple = [[i["ScientificName"]] for i in allBomb]

    with open('OrganismList.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(listTuple)
    return

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


def testNucleo():
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
        handle = Entrez.efetch(db="protein", id=ids, rettype='gb',retmode="xml")
        read = Entrez.read(handle)
        return read

    results = efetchProtein(proteins)
    nucleo_collection.insert_many(results)
    return


def efetchGenome(taxId):
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

        



def productsTable():
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






#######################
######---MAIN---#######
#######################

def main(args:dict) -> None:
    #unwrappingSpecies()
    #algo()
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