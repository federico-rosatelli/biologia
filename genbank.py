import csv
import os
from Bio import SeqIO
import json
import hashlib
import argparse
from pymongo import MongoClient
import sqlite3
from time import ctime,perf_counter
import requests
import pandas as pd

#######################
# List Address source db
# ARGS    = "https://shigen.nig.ac.jp/algae/download/downloadFile/Strain_e.txt"
# ARGS2   =
#######################

# Credits
# Authors: federico-rosatelli (Federico Rosatelli), AxnNxs (Mattia Di Gesaro)


# Riferimenti ufficiali
# biopython:    "https://biopython.org/wiki/Documentation"
# hashlib:      "https://docs.python.org/3/library/hashlib.html"
# argparse:     "https://docs.python.org/3/library/argparse.html"
# pymongo       "https://pymongo.readthedocs.io/en/stable/"

class Global:
    CLUSTER = "localhost:27017"
    WEBSOURCE = {
        "Algae":{
            "Web":"https://shigen.nig.ac.jp/algae/download/downloadFile/Strain_e.txt",
            "File":"algaeDatabase.csv"
        },
        "MicroAlgae":{
            "Web":"",
            "File":"microAlgaeProva.csv"
        }
    }
    SQL = {
        "Init":"init.sql",
        "Store":"Biologia.db"
    }


class  PrintWarning:
    # applicazione della classe di bcolors
    def __init__(self,type:int,error:str="") -> None:
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

    
    def stdout(self,*strings:any) ->None:
        self.error = ' '.join(str(i) for i in strings)
        plus = ""
        if self.error[:1] == "\n":
            plus = "\n"
            self.error = self.error[1:]
        print(bcolors.BOLD + f"{plus}[{ctime()}] " + bcolors.ENDC + self.color + self.error +  bcolors.ENDC)
    
    def __str__(self) -> str:
        return self.error

CLUSTER = "localhost:27017"     # apertura porta di default sul localhost che esegue il codice

class Parsing(object):
    """Parsing class for GenBank file format"""
    def __init__(self,file_name:str) -> None:
        
        self.record = SeqIO.parse(file_name, "genbank")

    def parsing_gene(self) -> list:
        """Ritorna due liste: 
        - nella prima vi sono i dati analizzati: per ogni file.gbk dato in input ci sono  di geni,
          e per ogni gene avremo un dizionario dei campi di interesse, come il Name, l'ID, e così via;
        - nella seconda vi sono i dati trattati riferiti alla Sequenza. Salviamo infatti sia la sequenza
          di basi azotate "as is", sia la versione codificata con SHA256."""  
        gene_array = []
        hex_array = []
        count = 0
        last = -1
        errors = 0
        PrintWarning(8).stdout("\nSTART PARSING")
        print("~"*56+"\n")
        t_last = perf_counter()
        for gene in self.record:
            count += 1
            if int((count/86100)*100) != last:
                t_now = perf_counter()
                print("[%-50s] %d%% %d" % ('='*((last+1)//2),last+1,(t_now-t_last)*(86100//count)),end="\r") #??????
                last +=1
            #print(count)
            json_gene = {}
            json_convert = {}
            json_gene["Name"] = gene.name
            json_gene["Id"] = gene.id
            check = False
            try:
                if gene.seq != "":
                    #print(count)
                    check = True
                    sh = self.sha_index(str(gene.seq))
                    json_gene["Seq_Hex"] = sh["hex"]  #str(gene.seq)-> SHA256 funzione
                    json_convert["Seq_Hex"] = sh["hex"]  #str(gene.seq)-> SHA256 funzione
                    json_convert["Seq_Raw"] = sh["seq"]  #str(gene.seq)-> SHA256 funzione
                    
            except Exception as e:
                if check:
                    print(e,"Entra qua nella sequenza:",gene.id)
                errors += 1
                pass
            json_gene["Description"] = gene.description
            json_gene["Features"] = []
            for feature in gene.features:
                feature_type = {"Type":feature.type}
                #json_gene["Features"][feature.type] = {}
                for qual in feature.qualifiers:
                    feature_type[qual] = feature.qualifiers.get(qual)[0]
                feature_type["Location"] = (int(feature.location.start),int(feature.location.end))
                json_gene["Features"].append(feature_type)
            hex_array.append(json_convert)
            gene_array.append(json_gene)
        PrintWarning(3).stdout("\nParsing completed\n") 
        PrintWarning(5).stdout("Process completed with",errors,'errors\n')
        return gene_array,hex_array
    
    def sha_index(self,seq:str) -> dict:
        # ritorna, in un dizionario, il codice hash di una sequenza di basi azotate (ACGTU le possibili unità)
        return {
            'hex':hashlib.sha256(seq.encode('utf-8')).hexdigest(),
            'seq':seq
        }
    
    def save_data(self, json_array:tuple) -> None:
        # crea un file json contentene i dati parsati - salvare in DB?
        with open("datastruct.json","w") as js:
            json.dump({
                        "struct":json_array[0],
                        "hex":json_array[1]
                        },js,indent=4)




class Database:
    def __init__(self) -> None:
        self.file = Global.SQL["Store"]
        self.client = None
        self.connection = None
        self.printWarning = PrintWarning(10)
        self.dataSource = {}
        self.totLength = 0
    
    # def save_on_mongo(self)
    def save_on_mongo(self,tuple_of_array:tuple) -> None:
        PrintWarning(8).stdout("\nStart saving on database")
        ip = Global.CLUSTER.split(":")[0]
        port = int(Global.CLUSTER.split(":")[1])
        self.client = MongoClient(f'{ip}', port)
        db = self.client["Biologia"]  # attenzione a quando si richiama il DB da riga di comando: case sensitive
        collection_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        count = 0
        last = -1
        for i in range(len(tuple_of_array[0])):
            count+=1
            if int((count/len(tuple_of_array[0]))*100) != last:
                print("[%-50s] %d%%" % ('='*((last+1)//2),last+1),end="\r")
                last +=1
            filter = {"Name":tuple_of_array[0][i]["Name"]}
            check_name = collection_data.find_one(filter)
            if not check_name:
                collection_data.insert_one(tuple_of_array[0][i])
            if tuple_of_array[1][i] != {}:
                filter = {"Seq_Hex":tuple_of_array[1][i]["Seq_Hex"]}
                check_hex = collection_convert.find_one(filter)
                if not check_hex:
                    collection_convert.insert_one(tuple_of_array[1][i])
        self.printWarning.stdout("Saved correctly on mongo")
        
    def save_on_sql(self,tuple_of_array:tuple,file=None):
        # filename to form database
        if file != None:
            self.file = file
        try:
            self.connection = sqlite3.connect(self.file)
            cursor = self.connection.cursor()
            if not self.file_exists():
                scr = open("init.sql").read()
                cursor.executescript(scr)
            for i in range(len(tuple_of_array[0])):
                data_tuple = (tuple_of_array[0][i]["Name"],tuple_of_array[0][i]["Id"],tuple_of_array[0][i]["Seq_Hex"],tuple_of_array[0][i]["Description"],i)
                find = f"SELECT Seq_Hex FROM Biologia WHERE Seq_Hex='{tuple_of_array[0][i]['Seq_Hex']}'"
                res = cursor.execute(find)
                if not res:
                    table ="""INSERT INTO Biologia VALUES(?,?,?,?,?)"""
                    cursor.execute(table,data_tuple)
        except Exception as e:
            self.connection.close()
            return PrintWarning(4).stdout(f"Error Database Biologia.db: {e}")
        self.connection.commit()
        if not self.connection:
            self.connection.close()
            return PrintWarning(5).stdout(f"Error Writing Database Biologia.db")
        self.connection.close()
        return None

    def get_data_from_mongo(self,info:dict,saveOnJson:bool=False,algae:bool=False,micro_algae:bool=False) -> dict:
        if self.client == None:
            ip = CLUSTER.split(":")[0]
            port = int(CLUSTER.split(":")[1])
            self.client = MongoClient(f'{ip}', port)
        db = self.client["Biologia"]
        collection_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        finder = collection_data.find(info)
        dataSource = {}
        
        for x in finder:
            #print(x["Features"])
            source = None
            for f in x["Features"]:
                if f["Type"] == "source":
                    source = f["organism"]
                    if source not in dataSource:
                        dataSource[source] = []
            dataSource[source].append(x["Id"])

        self.totLength = len(dataSource)
        self.dataSource = dataSource
        self.diffAlg = {}
        if saveOnJson:
            self.save_on_json()

        if algae:
            totAlgae,same = self.isAlgae(dataSource)
            self.diffAlg["Algae"] = same
            algaeSource = {}
            for algae in totAlgae:
                algaeSource[algae] = dataSource[algae]
            self.dataSource = algaeSource
            PrintWarning(2).stdout(f"Total algae {len(self.dataSource)} out of {self.totLength} with {len(same)} different gene")
            if saveOnJson:
                self.save_on_json(fileName="algaeSource.json")

        if micro_algae:
            # ? self.datasource vs. datasource
            totAlgae,same = self.isMicroAlgae(dataSource)
            self.diffAlg["MicroAlgae"] = same
            algaeSource = {}
            for algae in totAlgae:
                algaeSource[algae] = dataSource[algae]
            self.dataSource = algaeSource
            PrintWarning(2).stdout(f"Total micro algae {len(self.dataSource)} out of {self.totLength} with {len(same)} different gene")
            if saveOnJson:
                self.save_on_json(fileName="microAlgaeSource.json")

        #
        self.printDifferenceAlgae()
        dataDiff = self.checkDifferenceAlgae()
        print(len(dataDiff))

        return self.dataSource

    def printDifferenceAlgae(self) -> list:

        # list1 = [1,2,3,5,7]
        # list2 = [1,2,4,6,8]
        # list3 = [1,4,6]

        # listatot = [3,5,6,7,8]
        # 


        #print algae non contenute in microalgae
        diff = [self.diffAlg[key] for key in self.diffAlg]
        print(len(diff[0])-len(diff[1]))


    def checkDifferenceAlgae(self) -> list:
        #check algae non contenute in microalgae
        diff = [self.diffAlg[key] for key in self.diffAlg]
        tot = {}
        for i in range(len(diff)):
            for k in range(len(diff[i])):
                if diff[i][k] not in tot:
                    tot[diff[i][k]] = 1
                else:
                    tot[diff[i][k]] += 1  
        algTot = []
        for i in tot:
            if tot[i] == 1:
                algTot.append(i)
        return algTot


        
        # for i in range(len(diff)):
        #     print("lista ", i, " in stampa... ")
        #     diffsort = sorted(diff[i])
        #     for k in range(len(diffsort)):
        #         print(diffsort[k])
        # print("adesso i dati presenti soltanto una volta")
        # algTotsort = sorted(algTot)
        # print(algTotsort)
        # print(len(algTotsort))
        # return algTotsort

        # [1,2,3,4,5,6],[2,5,8,2,6],[1,5,3]
    
    def save_on_json(self,fileName:str="dataSource.json") ->None:
        with open(fileName,"w") as js:
            json.dump(self.dataSource,js,indent=4)
    
    def isAlgae(self,dataSource,rewrite:bool=False) ->list:
        if not os.path.isfile(Global.WEBSOURCE["Algae"]["File"]) or rewrite:
            r = requests.get(Global.WEBSOURCE["Algae"]["Web"])
            with open(Global.WEBSOURCE["Algae"]["File"], "w") as wr:
                dataContent = r.text
                dataContent = dataContent.split("\n")
                writer = csv.writer(wr)
                for i in range(len(dataContent)):
                    dataContent[i] = dataContent[i].split("\t")
                    writer.writerow(dataContent[i])
        file = open(Global.WEBSOURCE["Algae"]["File"],"r")
        dataAlgae = csv.reader(file)
        dataAlgae = list(dataAlgae)
        totAlgae = []
        unico = []
        for key in dataSource:
            keySplit = key.split(" ")
            if keySplit[0] == "cf.":
                val_key = key.split(" ")[1]
            else:
                val_key = key.split(" ")[0]

            for algae in dataAlgae:
                if len(algae) > 1 and val_key.lower() in algae[2].lower():
                    totAlgae.append(key)
                    if val_key.lower() not in unico:
                        unico.append(val_key.lower())
        return totAlgae,unico

    def isMicroAlgae(self,dataSource):
        rl = open(Global.WEBSOURCE["MicroAlgae"]["File"]).readlines()
        dataAlgae = []
        for i in range(len(rl)):
            data = rl[i].split(";")
            dataAlgae.append(data[1:])
        totAlgae = []
        unico = []
        for key in dataSource:
            keySplit = key.split(" ")
            if keySplit[0] == "cf.":
                val_key = key.split(" ")[1]
            else:
                val_key = key.split(" ")[0]

            for algae in dataAlgae:
                if len(algae) > 1 and val_key.lower() in algae[2].lower():
                    totAlgae.append(key)
                    if val_key.lower() not in unico:
                        unico.append(val_key.lower())
        return totAlgae,unico
        





    def file_exists(self)-> bool:
        return os.path.isfile(self.file) and not os.stat(self.file).st_size == 0


class bcolors:
    # codici di errore, intesa da utilizzare entro class Error
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





def main(args:dict) -> None:
    d = Database()
    if args["nosql_mongo"] and args["sqlite3"]:
        return PrintWarning(5).stdout("Can't select both mongo-db and sql for storing")
    if args["file"]:
        p = Parsing(args["file"])
        data,conv = p.parsing_gene()
        if args["json"]:
            p.save_data((data,conv))
        
        if args["nosql_mongo"]:
            d.save_on_mongo((data,conv))
        if args["sql"]:
            d.save_on_sql((data,conv))
    if args["find"]:
        finder = {"Features":{"$elemMatch":{"Type":"source"}}}
        save = False
        algae = False
        micro = False
        if args["json"]:
            save = True
        if args["algae"]:
            algae = True
        if args["micro_algae"]:
            micro = True
        data = d.get_data_from_mongo(finder,saveOnJson=save,algae=algae,micro_algae=micro)
    #d.isMicorAlgae()
    
    return



    

if __name__ == "__main__":
    # costrutto che permette di passare da riga di comando l'opzione desiderata per l'inserimento dei dati 
    # in un tipo di DB
    parser = argparse.ArgumentParser(
                    prog='Biologia Database Parsing',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('-n', '--nosql-mongo',
                        action='store_true')
    parser.add_argument('-s', '--sqlite3',
                        action='store_true')
    parser.add_argument('-j', '--json',
                        action='store_true')
    parser.add_argument('--find',
                        action='store_true')
    parser.add_argument('-a','--algae',
                        action='store_true')
    parser.add_argument('-m','--micro-algae',
                        action='store_true')
    parser.add_argument('-f', '--file')
    args = vars(parser.parse_args())
    main(args)
    



