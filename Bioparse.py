########################################################################################
# Libraries & API
########################################################################################

from pymongo import MongoClient
from Bio import SeqIO, Entrez
from time import ctime, perf_counter
from Bio.Seq import Seq
from matplotlib import pyplot as plt
from validate_email import validate_email
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

Global = ""                                                 # Deprecated?
client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_taxonomy_data = db["taxonomy_data"]
collection_nucleotide_data = db["nucleotide_data"]
collection_taxonomy_tree = db["taxonomy_tree"]
collection_table_basic = db["table_basic"]
collection_protein_data = db["protein_data"]
ignore_names = ["environmental samples"]                    # taxonomy ID to be avoided in collection_taxonomy_data
Entrez.api_key = "cc030996838fc52dd1a2653fad76bf5fe408"
Entrez.email = ""
regex = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,7}\b'

CLUSTER = "localhost:27017"
WEBSOURCE = {
    "Algae":{
        "Web":"https://shigen.nig.ac.jp/algae/download/downloadFile/Strain_e.txt",
        "File":"data/databaseCsv/algaeDatabase.csv"
    },
    "MicroAlgae":{
        "Web":"",
        "File":"data/databaseCsv/microAlgaeDatabase.csv"
    }
}
SQL = {
    "Init":"data/DBs/init.sql",
    "Store":"Biologia.db"
}
JSON = {
    "Path":"data/src/",
    "Type":"JSONFile"
}
ALIGNMENT = {
    "Path":"data/src/alignments/",
    "Type":"TXTFile"
}



########################################################################################
# Utility classes and methods
########################################################################################

def check_email() -> str:
    nb = input('Please insert email: ')
    while(re.fullmatch(regex, nb) == None):
        nb = input('Invalid Email, please retry: ')
        if (re.fullmatch(regex, nb) != None):
            break
    print(f"Valid Email, Passing --> {nb} <-- to NCBI platform")
    return nb

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
        return
    

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
        return

########################################################################################

class Alignment:
    '''Class for SMITH-WATERMAN algorithm implementation'''


    def __init__(self, *seqs:str, gap:int=1, show_table:bool=False) -> None:
        '''Constructor'''
        if len(seqs)<2:
            print(bcolors.FAIL+f"MinMaxError: input len for sequences must be >= 2, got {len(seqs)}"+bcolors.ENDC)
            return
        self.seqs = [Seq(i) for i in seqs]
        self.seq1 = seqs[0]
        self.seq2 = seqs[1]
        self.gap = gap
        self.show_table = show_table
    

    def __str__(self) -> str:
        '''Internal purposes'''
        return '\n'.join(self.seqs)


    def __len__(self) -> int:
        '''Internal purposes'''
        x = 1
        for i in self.seqs:
            x *= len(i)
        return x
    

    def __call__(self, seq:str) -> list:
        '''Internal purposes'''
        self.seqs.append(seq)
        return self.seqs


    def createScoreMatrix(self, lnSeq1:int, lnSeq2:int) -> tuple:
        '''Function that return a Matrix with score per gene'''
        score_matrix = [[0 for _ in range(lnSeq2 + 1)] for _ in range(lnSeq1 + 1)]
        max_score = 0
        max_index = (0,0)
        for i in range(1, lnSeq1 + 1):
            for j in range(1, lnSeq2 + 1):
                score = 1 if self.seq1[i-1] == self.seq2[j-1] else -1

                score_matrix[i][j] = max(
                    0,
                    score_matrix[i-1][j-1] + score,
                    score_matrix[i-1][j] - self.gap,
                    score_matrix[i][j-1] - self.gap
                )

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_index = (i, j)
        return score_matrix,max_index


    def localAlignment(self, save_table:bool=False) -> tuple:
        '''Function that use createScoreMatrix() to operate on genes and make appropriate localAlignment according to Smith-Waterman'''
        if len(self.seqs)>2:
            print(bcolors.WARN_BOX+f"Warning! Only 2 arguments were expected, but got {len(self.seqs)}.\n\t-The algorithm will use only the first 2 sequences..."+bcolors.ENDC)
        aligned_seq1 = ""
        aligned_seq2 = ""
        score_matrix, max_index_score_matrix = self.createScoreMatrix(len(self.seq1),len(self.seq2))
        i, j = max_index_score_matrix
        traceback = []
        print(f"Table length. First sequence:{len(self.seq1)}, Second sequence:{len(self.seq2)}")
        print(f"Max term: {score_matrix[i][j]}; in index: {max_index_score_matrix}")
        f1 = 0
        while score_matrix[i][j] != 0:
            traceback.append((i,j))
            if score_matrix[i-1][j-1] >= score_matrix[i-1][j] and score_matrix[i-1][j-1] >= score_matrix[i][j-1]:
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                i, j = i-1, j-1
                f1 += 1
            elif score_matrix[i-1][j] >= score_matrix[i-1][j-1] and score_matrix[i-1][j] >= score_matrix[i][j-1]:
                aligned_seq1 = self.seq1[i-1] + aligned_seq1
                aligned_seq2 = '-' + aligned_seq2
                i -= 1
            else:
                aligned_seq1 = '-' + aligned_seq1
                aligned_seq2 = self.seq2[j-1] + aligned_seq2
                j -= 1
        
        if self.show_table:
            printTable(score_matrix,(self.seq1,self.seq2),trace=traceback)
            if save_table:
                saveTable(score_matrix,trace=traceback)
        print(f"{(f1/len(aligned_seq1))*100}% of alignment")
        return aligned_seq1, aligned_seq2
    

    def globalAlignment(self):
        '''TEST check'''
        return 0


#FOR TESTING PURPOSE:
# a = Alignment("AGTCCCTGATTTAGTCCCTGATTTAGTATTTAGTCCCTGATTTAGTATTTAGTCCCTGATTTAGTCCCTGATTT",
# "TTTAGTCCCTGATTTAGTTTTAGTCCCTGATTTAGTTTTAGTCCCTGATTTAGT",show_table=True)
# a.localAlignment(save_table=True)


########################################################################################

class Parsing(object):
    '''Class that takes in input a .gbk file for parsing purposes.
    It has to be considered Legacy and needs to be modified before use'''


    def __init__(self, file_name:str=None) -> None:
        '''Constructor'''
        if not file_name:
            self.record = []
        else:
            self.record = SeqIO.parse(file_name, "genbank")


    def parsing_gene(self, verbose:bool=False) -> list:
        '''Method that returns two lists:
         - In "gene_array" there are the analyzed data: for each file.gbk given in input there are X genes,
           and for each gene we will have a dictionary of fields of interest, such as the Name, the ID, and so on;
         - In "hex_array" there are the processed data referred to the Sequence. In fact, let's save both the sequence
           of nitrogenous bases "as is", both the version encoded with SHA256'''
        gene_array = []
        hex_array = []
        count = 0
        last = -1
        errors = 0
        PrintWarning(8).stdout("\nSTART PARSING")               # messaggio informativo su console
        print("~"*56+"\n")
        t_last = perf_counter()
        
        for gene in self.record:
            count += 1
            if int((count/86100)*100) != last and verbose:              # Gestito in base al verbose, attenzione alla lunghezza del file (...)
                t_now = perf_counter()
                print("[%-50s] %d%% %d" % ('='*((last+1)//2),last+1,(t_now-t_last)*(86100//count)),end="\r") #??????
                last +=1
            #print(count)
            json_gene, json_convert, errors = self.parseGene(gene)
            hex_array.append(json_convert)
            gene_array.append(json_gene)
        PrintWarning(3).stdout("\nParsing completed\n") 
        PrintWarning(5).stdout("Process completed with", errors, 'errors\n')
        return gene_array, hex_array


    def parseGene(self, gene) -> tuple[dict, dict, int]:
        '''Method that returns a tuple structured as follows:
         - dict json_gene which contains Name and ID of the analyzed gene;
         - dict json_convert which contains both the "raw" and the translated SHA256 sequence;
         - int error which reports how many errors occurred during the search, i.e. how many times the data was missing
         This tuple is used inside the parsing_gene() method to populate the two arrays gee_array and hex_array'''
        json_gene = {}
        json_convert = {}
        json_gene["Name"] = gene.name
        json_gene["Id"] = gene.id
        check = False
        errors = 0
        try:
            if gene.seq != "":
                #print(count)
                check = True
                sh = self.sha_index(str(gene.seq))
                json_gene["Seq_Hex"] = sh["hex"]        #str(gene.seq)-> SHA256 funzione
                json_convert["Seq_Hex"] = sh["hex"]     #str(gene.seq)-> SHA256 funzione
                json_convert["Seq_Raw"] = sh["seq"]     #str(gene.seq)-> SHA256 funzione
        except Exception as e:
            if check:
                print(e,"Enters here:",gene.id)
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
        return json_gene, json_convert, errors
    

    def sha_index(self, seq:str) -> dict:
        '''Method which, given a sequence as input, returns the sequence and its hash code in a dictionary.
         ACGTU are the possible units of any sequence'''
        return {
            'hex':hashlib.sha256(seq.encode('utf-8')).hexdigest(),
            'seq':seq
        }
    

    def save_data(self, json_array:tuple) -> None:
        '''Function that creates a JSON file with parsed data'''
        with open(f"{Global.JSON['Path']}datastruct.json","w") as js:
            json.dump({
                        "struct":json_array[0],
                        "hex":json_array[1]
                        },js,indent=4)


########################################################################################


class Database:
    '''In this class are managed the different types of DB on which the data
     they can be saved. The task of this class is initialization
     and data manipulation, along with the ability to input
     credentials to communicate with the servers of the major web platforms.'''


    def __init__(self, verbose=False, type:str="", email:str=None) -> None:
        '''Constructor'''
        self.file = Global.SQL["Store"]
        ip = Global.CLUSTER.split(":")[0]
        port = int(Global.CLUSTER.split(":")[1])
        self.client = MongoClient(f'{ip}', port)
        self.connection = None
        self.mongo_collections = (type+"_data",type+"_hex")
        self.printWarning = PrintWarning(10)
        self.dataSource = {}
        self.totLength = 0
        self.diffAlg = {}
        self.verbose = verbose
        self.parse = Parsing()
        self.email = email
        print(email)
        if email:
            Entrez.email = email
    
    
    def listSource(self)->list:
        '''Returns a list of sources from which data is fetched.'''
        return [key for key in self.dataSource]
    

    def addEmail(self, email) -> None:
        '''Update the email field for identification.'''
        self.email = email
        Entrez.email = email


    def save_on_mongo(self, tuple_of_array:tuple) -> None:
        '''This method generates the database in MongoDB (see Requirements at the top of the file or the README.md).
         The data is saved in a local db called Biology, where all the operations will then be performed.'''
        PrintWarning(8).stdout("\nStart saving on database")
        db = self.client["Biologia"]  # attenzione a quando si richiama il DB da riga di comando: case sensitive
        print(self.mongo_collections)
        collection_taxonomy_data = db[self.mongo_collections[0]]
        collection_convert = db[self.mongo_collections[1]]
        count = 0
        last = -1
        for i in range(len(tuple_of_array[0])):
            count+=1
            if int((count/len(tuple_of_array[0]))*100) != last:
                print("[%-50s] %d%%" % ('='*((last+1)//2),last+1),end="\r")
                last +=1
            filter = {"Name":tuple_of_array[0][i]["Name"]}
            check_name = collection_taxonomy_data.find_one(filter)
            if not check_name:
                collection_taxonomy_data.insert_one(tuple_of_array[0][i])
            if tuple_of_array[1][i] != {}:
                filter = {"Seq_Hex":tuple_of_array[1][i]["Seq_Hex"]}
                check_hex = collection_convert.find_one(filter)
                if not check_hex:
                    collection_convert.insert_one(tuple_of_array[1][i])
        self.printWarning.stdout("Saved correctly on mongo")
    

    def save_one_in_mongo(self,struct:dict,collection:str)->None:
        '''Like the save_on_mongo() method, but with per-record sensitivity.
         In this way, we can manually add the data of interest.'''
        db = self.client["Biologia"]
        collection_taxonomy_data = collection + "_data"
        collection_hex_name = collection + "_hex"
        collection_taxonomy_data = db[collection_taxonomy_data]
        collection_hex = db[collection_hex_name]
        data = struct["data"]
        hex = struct["hex"]
        collection_taxonomy_data.insert_one(data)
        collection_hex.insert_one(hex)
        self.printWarning.stdout("Saved correctly on mongo")

        
    def save_on_sql(self,tuple_of_array:tuple, file=None):
        '''This method generates the database in SQL (see Requirements at the beginning of the file or the README.md).
         Although functional, at the moment this branch of the project is in STANDBY'''
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


    def get_data_from_mongo(self, info:dict, saveOnJson:bool=False, algae:bool=False, micro_algae:bool=False) -> dict:
        '''"Getter" type method, allows consultation of the DB based on MongoDB'''
        if self.client == None:
            ip = CLUSTER.split(":")[0]
            port = int(CLUSTER.split(":")[1])
            self.client = MongoClient(f'{ip}', port)
        db = self.client["Biologia"]
        collection_taxonomy_data = db[self.mongo_collections[0]]
        collection_convert = db[self.mongo_collections[1]]
        finder = collection_taxonomy_data.find(info)
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
        # self.printDifferenceAlgae()
        # dataDiff = self.checkDifferenceAlgae()
        #return self.dataSource
        self.confronto()
        return self.dataSource
    

    def alignmentSeq(self, f1_key:str, f2_key:str) -> str:
        '''Function that use Smith-Waterman Alignment to make comparison between sequences'''
        db = self.client["Biologia"]
        collection_taxonomy_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        info = {"Features":{"$elemMatch":{"Type":"source","organism":f1_key}}}
        # FIND PROTEIN
        finder = collection_taxonomy_data.find_one(info)
        if not finder:
            PrintWarning(5).stdout(f"Error searching {f1_key}: Organism Not Found")
            return None, None
        if not 'Seq_Hex' in finder:
            PrintWarning(5).stdout(f"Error searching {f1_key}: Hex Sequence Not Found")
            return None, None
        hex1 = finder['Seq_Hex']
        info2 = {"Features":{"$elemMatch":{"Type":"source","organism":f2_key}}}
        finder2 = collection_taxonomy_data.find_one(info2)
        if not finder2:
            PrintWarning(5).stdout(f"Error searching {f2_key}: Organism Not Found")
            return None, None
        if not 'Seq_Hex' in finder2:
            PrintWarning(5).stdout(f"Error searching {f2_key}: Hex Sequence Not Found")
            return None, None
        hex2 = finder2['Seq_Hex']
        info_hex = {"Seq_Hex":hex1}
        finder_hex = collection_convert.find_one(info_hex)
        if not finder_hex:
            PrintWarning(5).stdout(f"Error searching {f1_key}: Hex Sequence Not Found in the Database")
            return None, None
        seq1 = finder_hex['Seq_Raw']
        info_hex = {"Seq_Hex":hex2}
        finder_hex = collection_convert.find_one(info_hex)
        if not finder_hex:
            PrintWarning(5).stdout(f"Error searching {f2_key}: Hex Sequence Not Found in the Database")
            return None, None
        seq2 = finder_hex['Seq_Raw']
        smwt = Alignment(seq1,seq2,show_table=self.verbose)
        s1,s2 = smwt.localAlignment()
        return s1,s2
    

    def saveFileSeq(self, name:str, seq1:str, seq2:str) -> None:
        '''Function used in alignment() that given in input two sequence FASTA return a file with those sequences'''
        seq1_fasta = ""
        for i in range(len(seq1)):
            if i%100 == 0 and i != 0:
                seq1_fasta += "\n"
            elif i%10 == 0 and i != 0:
                seq1_fasta += " "
            seq1_fasta += seq1[i]
        seq2_fasta = ""
        for i in range(len(seq2)):
            if i%100 == 0 and i != 0:
                seq2_fasta += "\n"
            elif i%10 == 0 and i != 0:
                seq2_fasta += " "
            seq2_fasta += seq2[i]
        with open(name,"w") as wr:
            wr.write(seq1_fasta+"\n\n//\n\n"+seq2_fasta)
    

    def alignment(self, fileIn, fileOut) -> None:
        '''[DEPRECATED] Function that given a file sequence make an alignement'''
        rd = open(fileIn).readlines()
        for i in range(len(rd)):
            for j in range(i+1,len(rd)):
                PrintWarning(3).stdout(f"Align {rd[i].strip()} with {rd[j].strip()}")
                seq1,seq2 = self.alignmentSeq(rd[i].strip(),rd[j].strip())
                if seq1 and seq2:   
                    PrintWarning(3).stdout(f"Sequence length: {len(seq1)}\n")
                    self.saveFileSeq(f"{Global.ALIGNMENT['Path']}{fileOut.split('.')[0]}_{rd[i].strip()}_{rd[j].strip()}.txt",seq1,seq2)


    def printDifferenceAlgae(self) -> list:
        '''Print the algae not contained in microalgae'''
        diff = [self.diffAlg[key] for key in self.diffAlg]
        print(len(diff[0])-len(diff[1]))


    def checkDifferenceAlgae(self) -> list:
        '''Check the algae not contained in microalgae'''
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

    
    def save_on_json(self, fileName:str="dataSource.json") -> None:
        '''Save data to a file in .json format'''
        with open(f"{Global.JSON['Path']}{fileName}","w") as js:
            json.dump(self.dataSource,js,indent=4)
    

    def isAlgae(self,dataSource,rewrite:bool=False) -> tuple:
        ''' [DEPRECATED] This method compares with the data present in the db
        local and among the csv records in database/Csv, fetched
        from the European Database for accessibility reasons.
        Datasource is inserted in the path data/sourcejson/struct.txt ed
        is provided as an example.
        If it finds a match, it returns the Algae lists
        found and unique species (without any codes) '''
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


    def isMicroAlgae(self,dataSource) -> tuple:
        ''' [DEPRECATED] This method compares with the data present in the db
        local and among the csv records in database/Csv, fetched
        from the European Database for accessibility reasons.
        Datasource is inserted in the path data/sourcejson/struct.txt ed
        is provided as an example.
        If it finds a match, it returns the MicroAlgae lists
        found and unique species (without any codes) '''
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
    

    def toFasta(self,name):
        '''Test'''
        print(name)
        db = self.client["Biologia"]
        collection_taxonomy_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        info = {"Features":{"$elemMatch":{"Type":"source","organism":name}}}
        # FIND PROTEIN
        finder = collection_taxonomy_data.find(info)
        if not finder:
            PrintWarning(5).stdout(f"Error searching {name}: Organism Not Found")
            return None, None
        cds = []

    
    def proteinFind(self, id) -> dict:
        '''Matching con db esistente di Nucleotide'''
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        info = {"Id":id,"Features":{"$elemMatch":{"Type":"CDS"}}}
        finder_data = collection_data_nucleotide.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {id}: CDS Not Found")
            return None
        collection_taxonomy_data = db["protein_data"]
        collection_convert = db["protein_hex"]
        protein_id = None
        for f in finder_data["Features"]:
            if f["Type"] == "CDS":
                if "protein_id" in f:
                    protein_id = f["protein_id"]
        if not protein_id:
            PrintWarning(5).stdout(f"Error searching {id}: Protein Not Found")
            return None
        info = {"Id":protein_id}
        finder_data = collection_taxonomy_data.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {protein_id}: Protein Not Found In Database...","\n","\t\tSearching on NCBI...")
            p1,p2 = self.ncbiSearch(protein_id,"protein")
            dataFind = {
                'data':p1,
                'hex':p2
            }
            self.save_one_in_mongo(dataFind,"protein")
            return dataFind
        info = {"Seq_Hex":finder_data["Seq_Hex"]}
        finder_hex = collection_convert.find_one(info)
        if not finder_hex:
            PrintWarning(5).stdout(f"Error searching {finder_data['Seq_Hex']}: Hex Not Found")
            return None
        dataFind = {
            'data':finder_data,
            'hex':finder_hex
        }
        return dataFind
    

    def taxonFind(self,id) -> dict:
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        info = {"Id":id,"Features":{"$elemMatch":{"Type":"source"}}}
        finder_data = collection_data_nucleotide.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {id}: Source Not Found")
            return None
        collection_convert = db["taxonomy_hex"]
        taxon_meta = None
        for f in finder_data["Features"]:
            if f["Type"] == "source":
                if "db_xref" in f:
                    taxon_meta = f["db_xref"]
        if not taxon_meta:
            PrintWarning(5).stdout(f"Error searching {id}: db_xref Not Found")
            return None
        taxon_id = taxon_meta.split(":")[1]
        info = {"TaxId":taxon_id}
        finder_data = collection_taxonomy_data.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {taxon_id}: Taxonomy Not Found In Database...","\n","\t\tSearching on NCBI...")
            dataFind = self.ncbiSearch(taxon_id,"taxonomy")
            for data in dataFind:
                collection_taxonomy_data.insert_one(data)
            return dataFind
        return finder_data
    

    def genomeFind(self,id) -> dict:
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        info = {"Features":{"$elemMatch":{"Type":"CDS","db_xref":{"$exists":True}}}}
        finder_data = collection_data_nucleotide.find(info)
        genes = []
        for finder in finder_data:
            gene_id = None
            for f in finder["Features"]:
                if f["Type"] == "gene":
                    if "gene" in f:
                        gene_id = f["gene"]
            if gene_id != None:
                print(gene_id)
                dataFind = self.ncbiSearchGenome(gene_id,"genome")
                print(dataFind)
                # info = {"TaxId":gene_id}
                # finder_data = collection_taxonomy_data.find_one(info)
        return finder_data


    def ritornodicose(self) -> None:
        finder = {"Id":"","Features":{"$elemMatch":{"Type":"CDS"}}}
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        data = collection_data_nucleotide.count_documents(finder)
        print(data)


    def confronto(self): #CAMBIA NOME
        for key in self.dataSource:
            for id in self.dataSource[key]:
                data = self.proteinFind(id)
                if not data:
                    PrintWarning(5).stdout(f"ID:{id}")
                    #return None
                else:
                    PrintWarning(3).stdout(f"Protein ID:{id}")
                
                data = self.taxonFind(id)
                if not data:
                    PrintWarning(5).stdout(f"ID:{id}")
                    #return None
                else:
                    PrintWarning(3).stdout(f"Taxon ID:{id}")     
        #data = self.genomeFind(4)
                # if not data:
                #     PrintWarning(5).stdout(f"ID:{id}")
                #     #return None
                # else:
                #     PrintWarning(3).stdout(f"Genome ID:{id}")

    def ncbiSearch(self,id:str,database:str) -> tuple:
        # NB: in futuro la composizione del link potrebbe cambiare nella sua struttura, in base alla gestione interna di NCBI.
        # TO-DO: se questo metodo non ritorna i dati correttamente andrà aggiornata la procedura di reperimento dei dati.
        handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        p1,p2,errors = self.parse.parseGene(record)
        return p1,p2


    def ncbiSearchTaxon(self,id:str,database:str) ->tuple:
        try:
            handle = Entrez.efetch(db=database, id=id, retmode="xml")
            read = Entrez.read(handle)
        except Exception as e:
            print(e)
            os.sleep(20)
        return read
    

    def ncbiSearchGenome(self,id:str,database:str) ->tuple:
        handle = Entrez.efetch(db=database, id=id, retmode="xml")
        read = Entrez.read(handle)
        #print(read)
        # handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
        # record = SeqIO.read(handle, "genbank")
        return read


########################################################################################

def file_exists(self)-> bool:
    '''check file existence in path'''
    return os.path.isfile(self.file) and not os.stat(self.file).st_size == 0


def csvWrite(dataResult):
    '''Custom function that, using csv library, return a file with field of interests
    from a given list of data, parsed as NCBI Json structure'''
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
    return


def printTable(table, gene, trace=[]):
    '''Function implemented in Smith-Waterman that help to display Score Matrix '''
    #print([i for i in range(len(table[0]))])
    indx = "  "
    for i in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
        indx += gene[1][i]+" "*3
    print(indx)
    print('_'*10*4)
    for i in range(trace[len(trace)-1][0],trace[len(trace)-1][0]+10):
        print("| ",end="")
        for j in range(trace[len(trace)-1][1],trace[len(trace)-1][1]+10):
            #print(table[i][j])
            if (i,j) in trace:
                print(bcolors.OKGREEN+str(table[i][j])+bcolors.ENDC+" "*(2-len(str(table[i][j]))),end="| ")
            else:
                print(str(table[i][j])+" "*(2-len(str(table[i][j]))),end="| ")
        if i != 0:
            print(gene[0][i])
        else:
            print("")
    print('\n')
    return


def saveTable(table, trace=[]):
    '''Aid function that is used in Smith-Waterman and serve the purpose to save in .png format
    the result of printTable()'''
    if len(table)<40 or len(table[0])<40:
        return
    tableCopy = []
    for i in range(10,40):
        tbCp = []
        for j in range(10,40):
            tbCp.append(table[i][j])
        tableCopy.append(tbCp)
    color = [["w" for j in range(len(tableCopy[0]))]for i in range(len(tableCopy))]
    for i in range(len(tableCopy)):
        for j in range(len(tableCopy[i])):
            if (i+10,j+10) in trace:
                color[i][j] = "#56b5fd"
    # for i in range(10,20):
    #     color[trace[len(trace)-i-1][0]][trace[len(trace)-i-1][1]] = "#56b5fd"
    fig,ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(tableCopy)
    ax.table(cellText=df.values,cellColours=color,colLabels=df.columns,loc='center')
    #fig.tight_layout()
    plt.savefig('table.png',bbox_inches='tight')
    #plt.show()
    return


def SpecieProductMaker():
    '''Function that scan local DB in collection taxonomy_data and generate a file
    lists/SpecieProduct.csv that contains all product (proteins, transcriptomes,
    genes etc.) listed per species.'''
    filter = {"Features.product": {"$exists": True}}
    projection= {"_id":0, "Features.organism": 1, "Features.product": 1}
    dataResult = collection_taxonomy_data.find(filter, projection)
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
    filter = {"GBSeq_feature-table.GBFeature_key":"CDS","GBSeq_feature-table.GBFeature_quals.GBQualifier_name":"product"}
    proj = {"GBSeq_feature-table":1,"GBSeq_organism":1}
    dataResult = collection_nucleotide_data.find(filter, proj)


    def jsonList():
        '''Inner function (not used) that create a json file from the precedent processed data'''
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
        return
    
    return



#################################################################################
# Online Query methods (NCBI), require login method to be saved before
#################################################################################

def ncbiSearch(name:str, db="nucleotide") ->list:
    '''Function that submit queries with given a ScientificName to NCBI platform, section Nucleotide.
    Return data read. Default DB selected is nucleotide, can be changed referring to NCBI.'''
    # handle = Entrez.efetch(db="taxonomy", Lineage=name, retmode="xml")
    # read = Entrez.read(handle)
    handle = Entrez.esearch(db, term=name, rettype='gb', retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    print(f"Len of IDLIST:{len(record['IdList'])}")
    if len(record["IdList"]) == 0:
        raise Exception("List Empty")
    handle = Entrez.efetch(db, id=record["IdList"], rettype='gb',retmode="xml",complexity=1)
    read = Entrez.read(handle)
    print(f"Len of EFETCH:{len(read)}")
    return read


def finderTaxon(name):
    '''Function that, given a Taxonomy ID, add the relative specie to collection_taxonomy_data branch in Biologia DB,
    if it is not already present in list'''
    if name in ignore_names:
        return
    try:
        taxon = ncbiSearch(f"{name}[next level]", "taxonomy")
        for tax in taxon:
            if (not collection_taxonomy_data.find_one({"TaxId":tax["TaxId"]}) and tax["Rank"] == "species"):
                    print(tax["ScientificName"])
                    collection_taxonomy_data.insert_one(tax)
            finderTaxon(tax["ScientificName"])
    except Exception as e:
        return


def genomeRetrieve():
    '''Function that search on NCBI, for all occurrencies in taxonomy_data, if
    there are links to GFF, GBFF and FNA files containing genomic sequences
    for eache specie; in that case, carry out the download in ./data/temp
    folder, unwrap the files and move the extracted content in ./data/genomes/<txid>
    where <txid> is the relative id found in taxonomy_data for that specie'''
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


def nucleoImport():
    '''TESTING'''
    taxon_collection = db["taxonomy_data"]
    nucleo_collection = db["nucleotide_organism"]
    records = ncbiSearch("Scenedesmus bijugus")
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
            records = ncbiSearch(organism)
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


def nucleoImport():
    '''Sub-function that add new entries in nucleotide_data collection searching for them
    on NCBI platform'''
    daora = False
    nucleo_collection = db["nucleotide_data"]
    all_data = collection_taxonomy_data.find({},{"TaxId":1,"_id":0})
    for data in all_data:
        if f"txid{data['TaxId']}" == "txid1411642":
            daora = True
        else:
            print("No",f"txid{data['TaxId']}")
        if daora:
            print(f"txid{data['TaxId']}")
            try:
                insertDatas = ncbiSearch(f"txid{data['TaxId']}[Organism:exp]")
                for ins in insertDatas:
                    ins.pop("GBSeq_sequence",None)
                    nucleo_collection.insert_one(ins)
                #nucleo_collection.insert_many(insertDatas)
            except Exception as e:
                print(e)


def nucleoResult():
    '''Function that lightens the records of nucleotide_data collection taking only
    datas of interest.'''
    nucleoImport()
    # "txid257627"
    # "txid257627"
    nucleo_collection = db["nucleotide_data"]
    tot = ncbiSearch("txid257627[Organism:exp]")
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


def finderTaxonFromFile(fn):
    '''Function that, given a file path in input, scan for species and datas.
    File has to be formatted like those present in path ./data./databaseCsv'''
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


def unwrappingSpecies():
    '''Function that, using finderTaxonFromFile scan, nake a new OrganismList.csv file
    that contain contents parsed and cleaned in a readable format for humans.'''
    finderTaxonFromFile('data/databaseCsv/microAlgaeDatabase.csv')
    taxons = finderTaxon('data/databaseCsv/microAlgaeDatabase.csv')
    collection_taxonomy_data.delete_many({"Lineage":{"$regex":"environmental samples"}})
    parsed_file = collection_taxonomy_data.find({},{"ScientificName":1,"_id":0})
    listTuple = [[i["ScientificName"]] for i in parsed_file]
    with open('./lists/OrganismList.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(listTuple)
    return


def genusList():
    '''Function that generate a .csv file that list all Ranks present in
    taxonomy_tree collection, taking off some of them by avoiding to
    insert non-microalgae'''
    genus = collection_taxonomy_tree.find({"Rank":"genus"})
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
    taxon_collection = collection_taxonomy_data
    findTax = taxon_collection.find({},{"TaxId":1,"_id":0})
    for data in findTax:
        tax = f"txid{data['TaxId']}"
        print(tax)
        try:
            allData = ncbiProtein(tax)
            for prot in allData:
                prot.pop("GBSeq_sequence",None)
                collection_protein_data.insert_one(prot)
        except Exception as e:
            print(e)


def findTaxon(key):
    '''TEST'''
    rgx = re.compile(f'{key}', re.IGNORECASE)
    info = {"Lineage":rgx}
    dataFind = collection_taxonomy_data.find(info)
    for i in dataFind:
        print(i)



#################################################################################
# Offline Query methods (MongoDB)
#################################################################################


def productsTable():
    '''Function that take data from lists/SpecieProduct.csv and parse
    new collections (table_complete and table_basic) with those datas
    [Legacy procedure]'''
    products = open("./lists/SpecieProduct.csv").readlines()
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


def searchForRank(name = '', rank = ''):
    '''Function that search for input with a mongo-query in local DB'''
    filter = {}
    if (name != '' and rank != ''):
        filter = {name, rank}
    elif (name != ''):
        filter = {name}
    elif (rank != ''):
        filter = {rank}
    dataResult = collection_taxonomy_data.find(filter)
    return


def taxTreeMaker():
    '''Function that retrieves datas from collection_taxonomy_data of MongoDB and
    create a new collection, taxonomy_tree, that contain the Lineage for
    each specie'''
    results = collection_taxonomy_data.find({})
    iter = 0
    for res in results:
        print(iter)
        iter += 1
        # if collection_taxonomy_tree.find_one({"TaxId":res["ParentTaxId"]}):
        #     dataPush = {
        #             "TaxId":LineageEx["TaxId"],
        #             "Rank":LineageEx["Rank"],
        #             "ScientificName":LineageEx["ScientificName"],
        #             "SubClasses":[]
        #         }
        #     collection_taxonomy_tree.insert_one(dataPush)
        for i in range(1, len(res["LineageEx"])):
            #for LineageEx in res["LineageEx"]:
            LineageEx = res["LineageEx"][i]
            if not collection_taxonomy_tree.find_one({"TaxId":LineageEx["TaxId"]}):
                dataPush = {
                    "TaxId":LineageEx["TaxId"],
                    "Rank":LineageEx["Rank"],
                    "ScientificName":LineageEx["ScientificName"],
                    "SubClasses":[]
                }
                collection_taxonomy_tree.insert_one(dataPush)
            if (i == 1):
                continue
            LineageExPrev = res["LineageEx"][i-1]

            if collection_taxonomy_tree.find_one({"TaxId":LineageExPrev["TaxId"], "SubClasses":{"$elemMatch":{"TaxId":LineageEx["TaxId"]}}}):
                continue
            newDataPush = {
                    "TaxId":LineageEx["TaxId"],
                    "Rank":LineageEx["Rank"],
                    "ScientificName":LineageEx["ScientificName"],
                    }
            collection_taxonomy_tree.update_one({"TaxId":LineageExPrev["TaxId"]},{"$push":{"SubClasses":newDataPush}})
        if not collection_taxonomy_tree.find_one({"TaxId":res["ParentTaxId"], "SubClasses":{"$elemMatch":{"TaxId":res["TaxId"]}}}):
            newDataPush = {
                    "TaxId":res["TaxId"],
                    "Rank":res["Rank"],
                    "ScientificName":res["ScientificName"],
                    }
            collection_taxonomy_tree.update_one({"TaxId":res["ParentTaxId"]},{"$push":{"SubClasses":newDataPush}})


def parseToBasic() -> list:
    '''TESTING'''
    taxon_collection = db["taxonomy_data"]
    nucleo_collection = db["nucleotide_data"]
    tableBasic_collection = db["table_basic"]
    protein_collection = collection_protein_data
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


def newCollectionBene():
    '''TESTING'''
    new_collection = db["table_complete1"]
    old_collection = collection_taxonomy_tree
    nucleotide_collection = db["nucleotide_data"]
    protein_collection = collection_protein_data
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


def updateByGenomes(datasTaxon=''):
    '''Function that take a list where FNA files were found on NCBI platform
    through genomeRetrieve() and insert those link on local DB. It can be
    speciefied a new list in input'''
    if datasTaxon == '':
        datasTaxon = [(2961,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2961",True,True,False),
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
            (2730355,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2730355",True,True,False)]
    for data in datasTaxon:
        dataPush = {
            "Link":data[1],
            "GBFF":data[2],
            "FNA":data[3],
            "GFF":data[4]
        }
        collection_table_basic.update_one({"TaxId":str(data[0])},{"$push":{"Genomes":dataPush}})
    return































#DEPRECATED: genbank.py






    

















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








    # genus = collection_taxonomy_data.find({"Rank":"genus","Division":{"$not":re.compile("Bacteria")}})
    # """SELECT * FROM COLLECTION WHERE RANK=genus AND NOT =bacteria"""
    #datas = collection_taxonomy_tree.find_one({"ScientificName":"unclassified Chlorella"},{'_id':0})

    # taxon = ncbiSearch("chlorella[next level]", "taxonomy")
    # with open("taxonMeta.json","w") as jsw:
    #     json.dump(taxon,jsw,indent=4)
    #print(taxon)






########################################################################################
########################################################################################

def main(args:dict) -> None:
    print("Welcome! This program is intended to be used for parsing data from NCBI or a file.\n"
          "After parsing process is completed, you can check and query results in the DB selected as argument.")
    count = 0
    if len(sys.argv) == 1:
        count += 1
        print("Please add an argument at least.\n")
        if (count >= 3):
            print("Need a help? Digit -h or --help after calling genbank.py.\n")
    
    name = input('Please insert something to search: ')
    
    if (True):   #qui andrà il check per i metodi che fanno ricerche online
        Entrez.email = check_email()
    #unwrappingSpecies()
    #taxTreeMaker()                 #previously named algo()
    #searchForRank()                #previously named aglo()
    #testNucleo()
    #genomeRetrieve()               #previously named GFF()
    #nucleoResult()
    #genomeFind(name)
    #parseToBasic()
    #fattoBene() 
    #proteinFind()
    #newCollectionBene()
    #productsTable()
    #csvWrite(dataResult)
    #updateByGenomes()
    #findTaxon("Eukaryota")
    return



# Driver Code
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Biologia Database Parsing',
        description='What the program does',
        epilog='Text at the bottom of help')
    args = vars(parser.parse_args())
    main(args)


# # # # # def main(args:dict) -> None:
# # # # #     v = False
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
# # # # #     #d.isMicroAlgae()
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