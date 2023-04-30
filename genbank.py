# Libraries & API

import csv
import os
import random
import json
import hashlib
import argparse
import sqlite3
import requests
import pandas as pd
import smith_waterman
from Bio import SeqIO, Entrez
from pymongo import MongoClient
from time import ctime, perf_counter


#######################
# List Address source db
# ARGS    = "https://shigen.nig.ac.jp/algae/download/downloadFile/Strain_e.txt"
# ARGS2   = 
#######################

# Credits
# Authors: federico-rosatelli (Federico Rosatelli), AxnNxs (Mattia Di Gesaro)

# Riferimenti ufficiali & Requisiti
# Anaconda:     "https://www.anaconda.com"
# Biopython:    "https://biopython.org/wiki/Documentation"
# Hashlib:      "https://docs.python.org/3/library/hashlib.html"
# Argparse:     "https://docs.python.org/3/library/argparse.html"
# Pymongo:      "https://pymongo.readthedocs.io/en/stable/"

# Fonti ufficiali:
# 
# Shigen:       "https://shigen.nig.ac.jp"         DB giapponese
# NCBI:         "https://www.ncbi.nlm.nih.gov/"    DB americano
#
# La maggior parte dei dati presenti localmente nel Database interno
# sono stati scaricati manualmente a causa delle restrizioni sui
# download.



# Altri link
# Drive:        "https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW"
# Github:       "https://github.com/federico-rosatelli/biologia"



# NOTES
#
# MONGODB conserva i dati implicitamente in una memoria virtuale. Per trasportarli da un sistema
# a un altro è necessario utilizzare il comando mongodump per generare una cartella contenente
# il db di interesse e mongorestore sulla nuova postazione, una volta importata la cartella generata.
#
#
#
# La rappresentazione FASTA e FASTQ
# 
# FASTA conserva soltanto la Sequenza di nucleotidi o amminoacidi, codificando ogni gene in singole lettere
# per indice di posizione. Nella rappresentazione in Genbank, troviamo tale dato nel file JSON che salviamo
# in locale, alla voce translation per ogni Coding Sequence sotto ogni Specie, secondo la seguente gerarchia:
# 
# SPECIE
#   FEATURES
#       CDS
#           /translation="LSLAVGTTITLASYHWLL[...]""
#
# FASTQ è un "quality score" che associa alla sequenza, per ogni indice di posizione, un valore 
# qualitativo codificato in ASCII. Un esempio a seguire:
# @SRR64[...]       Name Sequence
# CCTCGTCTA[...]    DNA Sequence
# +SRR64[...]       Quality address
# BBBBBFFFF[...]    Quality Score
#
#
#
# Sequenziamento genetico
# E' il processo di determinazione dell'ordine dei nucleotidi (Adenina, Citosina, Guanina e Timina) che 
# costituiscono il frammento di DNA in analisi. Le tecniche principali di Squenziamento sono
# Sanger e NGS.
#
#
#
# Allineamento genetico
# E' il processo di confronto di due o più sequenze di DNA o proteine per identificare regioni identiche
# o simili per individuare eventuali relazioni funzionali, strutturali o filogenetiche.
# Le tecniche di allineamento prevedono il confronto globale e locale.
# Un'applicazione algoritmica di allineamento locale è data da Smith Waterman.
# Un'applicazione algoritmica di allineamento globale é data da Needleman-Wunsch.


# Global constants
#apertura porta di default sul localhost che esegue MongoDB
CLUSTER = "localhost:27017"     


#######################


class Global:
    '''Per l'accesso globale a costanti e risorse. I link qui preseti prelevano i dati direttamente
    dalle piattaforme web che forniscono i database contenenti i dati genomici.'''


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
        "Init":"init.sql",
        "Store":"Biologia.db"
    }
    JSON = {
        "Path":"data/sourceJson/",
        "Type":"JSONFile"
    }
    ALIGNMENT = {
        "Path":"data/alignments/",
        "Type":"TXTFile"
    }


#######################


class bcolors:
    '''Codici di errore, intesa da utilizzare entro class Error e PrintWarning'''


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


#######################


class  PrintWarning:
    """ Applicazione della classe di bcolors. Con questa classe possiamo stampare a video
    vari codici di errore/warning con colori diversi e tenere traccia visivamente di vari eventi."""


    def __init__(self, type:int, error:str="") -> None:
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

    
    def stdout(self, *strings:any) -> None:
        '''Metodo custom di stampa a video dei nostri codici di errore, in cui riportiamo il timestamp
        e la descrizione dell'errore con un colore fissato per utilizzo informativo. '''
        self.error = ' '.join(str(i) for i in strings)
        plus = ""
        if self.error[:1] == "\n":
            plus = "\n"
            self.error = self.error[1:]
        print(bcolors.BOLD + f"{plus}[{ctime()}] " + bcolors.ENDC + self.color + self.error +  bcolors.ENDC)
    

    def __str__(self) -> str:
        '''Utilizzato internamente per richiamare gli errori dal Costruttore'''
        return self.error


#######################


class Parsing(object):
    '''Classe che si occupa di trattare i dati forniti da GenBank.
    Il file su cui si opera si da essere per scontato in formato .gbk'''


    def __init__(self, file_name:str=None) -> None:
        '''Costruttore'''
        if not file_name:
            self.record = []
        else:
            self.record = SeqIO.parse(file_name, "genbank")


    def parsing_gene(self, verbose:bool=False) -> list:
        '''Metodo che ritorna due liste:
        - In "gene_array" vi sono i dati analizzati: per ogni file.gbk dato in input ci sono X geni,
          e per ogni gene avremo un dizionario dei campi di interesse, come il Name, l'ID, e così via;
        - In "hex_array" vi sono i dati trattati riferiti alla Sequenza. Salviamo infatti sia la sequenza
          di basi azotate "as is", sia la versione codificata con SHA256.'''
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
        '''Metodo che ritorna una tupla strutturata nel modo seguente:
        - dict json_gene        che contiene Nome e ID del gene analizzato;
        - dict json_convert     che contiene sia la sequenza "raw" che tradotta in SHA256;
        - int  error            che segnala quanti errori si sono verificati durante la ricerca, ovvero quante volte il dato è risultato mancante       
        Questa tupla viene utilizzata all'interno del metodo parsing_gene() per popolare i due array gee_array e hex_array.'''
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
        return json_gene, json_convert, errors
    

    def sha_index(self, seq:str) -> dict:
        '''Metodo che, dato in input una squenza, ritorna in un dizionario la sequenza e il suo codice hash.
        ACGTU sono le possibili unità di una qualsiasi sequenza. '''
        return {
            'hex':hashlib.sha256(seq.encode('utf-8')).hexdigest(),
            'seq':seq
        }
    

    def save_data(self, json_array:tuple) -> None:
        '''Crea un file json contentene i dati parsati.'''
        with open(f"{Global.JSON['Path']}datastruct.json","w") as js:
            json.dump({
                        "struct":json_array[0],
                        "hex":json_array[1]
                        },js,indent=4)


#######################


class Database:
    '''In questa classe vengono gestiti i diversi tipi di DB su cui i dati
    possono essere salvati. Compito di questa classe é l'inizializzazione
    e la manipolazione dei dati, insieme alla possibilità di inserimento
    delle credenziali per dialogare con i server delle maggiori piattaforme web.'''


    def __init__(self, verbose=False, type:str="", email:str=None) -> None:
        '''Costruttore'''
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
        '''Restituisce una lista delle fonti da cui vengono prelevati i dati.'''
        return [key for key in self.dataSource]
    

    def addEmail(self, email) -> None:
        '''Aggiorna il campo email per identificazione.'''
        self.email = email
        Entrez.email = email


    def save_on_mongo(self, tuple_of_array:tuple) -> None:
        '''Questo metodo genera il database in MongoDB (v. Requisiti all'inizio del file o il README.md).
        I dati vengono salvati in un db locale di nome Biologia, dove poi verranno eseguite tutte le operazioni.'''
        PrintWarning(8).stdout("\nStart saving on database")
        db = self.client["Biologia"]  # attenzione a quando si richiama il DB da riga di comando: case sensitive
        print(self.mongo_collections)
        collection_data = db[self.mongo_collections[0]]
        collection_convert = db[self.mongo_collections[1]]
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
    

    def save_one_in_mongo(self,struct:dict,collection:str)->None:
        '''Come il metodo save_on_mongo(), ma con sensibilità per singolo record.
        In questo modo, possiamo aggiungere manualmente i dati di interesse.'''
        db = self.client["Biologia"]
        collection_data_name = collection + "_data"
        collection_hex_name = collection + "_hex"
        collection_data = db[collection_data_name]
        collection_hex = db[collection_hex_name]
        data = struct["data"]
        hex = struct["hex"]
        collection_data.insert_one(data)
        collection_hex.insert_one(hex)
        self.printWarning.stdout("Saved correctly on mongo")

        
    def save_on_sql(self,tuple_of_array:tuple, file=None):
        '''Questo metodo genera il database in SQL (v. Requisiti all'inizio del file o il README.md).
        Sebbene funzionante, al momento questo ramo del progetto è in STANDBY'''
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
        '''Metodo di tipo "Getter", permette la consultazione del DB basato su MongoDB'''
        if self.client == None:
            ip = CLUSTER.split(":")[0]
            port = int(CLUSTER.split(":")[1])
            self.client = MongoClient(f'{ip}', port)
        db = self.client["Biologia"]
        collection_data = db[self.mongo_collections[0]]
        collection_convert = db[self.mongo_collections[1]]
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
        #self.confronto()
        return self.dataSource
    

    def alignmentSeq(self, f1_key:str, f2_key:str) -> str:
        db = self.client["Biologia"]
        collection_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        info = {"Features":{"$elemMatch":{"Type":"source","organism":f1_key}}}
        # FIND PROTEIN
        finder = collection_data.find_one(info)
        if not finder:
            PrintWarning(5).stdout(f"Error searching {f1_key}: Organism Not Found")
            return None, None
        if not 'Seq_Hex' in finder:
            PrintWarning(5).stdout(f"Error searching {f1_key}: Hex Sequence Not Found")
            return None, None
        hex1 = finder['Seq_Hex']
        info2 = {"Features":{"$elemMatch":{"Type":"source","organism":f2_key}}}
        finder2 = collection_data.find_one(info2)
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
        smwt = smith_waterman.Alignment(seq1,seq2,show_table=self.verbose)
        s1,s2 = smwt.localAlignment()
        return s1,s2
    

    def saveFileSeq(self, name:str, seq1:str, seq2:str) -> None:
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
    

    def alignment(self,fileIn,fileOut) -> None:
        rd = open(fileIn).readlines()
        for i in range(len(rd)):
            for j in range(i+1,len(rd)):
                PrintWarning(3).stdout(f"Align {rd[i].strip()} with {rd[j].strip()}")
                seq1,seq2 = self.alignmentSeq(rd[i].strip(),rd[j].strip())
                if seq1 and seq2:   
                    PrintWarning(3).stdout(f"Sequence length: {len(seq1)}\n")
                    self.saveFileSeq(f"{Global.ALIGNMENT['Path']}{fileOut.split('.')[0]}_{rd[i].strip()}_{rd[j].strip()}.txt",seq1,seq2)


    def printDifferenceAlgae(self) -> list:
        '''Stampa le algae non contenute in microalgae'''
        diff = [self.diffAlg[key] for key in self.diffAlg]
        print(len(diff[0])-len(diff[1]))


    def checkDifferenceAlgae(self) -> list:
        '''Effettua il check delle algae non contenute in microalgae'''
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
        with open(f"{Global.JSON['Path']}{fileName}","w") as js:
            json.dump(self.dataSource,js,indent=4)
    

    def isAlgae(self,dataSource,rewrite:bool=False) -> tuple:
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
        print(name)
        db = self.client["Biologia"]
        collection_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        info = {"Features":{"$elemMatch":{"Type":"source","organism":name}}}
        # FIND PROTEIN
        finder = collection_data.find(info)
        if not finder:
            PrintWarning(5).stdout(f"Error searching {name}: Organism Not Found")
            return None, None
        cds = []


    def file_exists(self)-> bool:
        return os.path.isfile(self.file) and not os.stat(self.file).st_size == 0

    
    def proteinFind(self, id) -> dict:
        '''Matching con db esistente di Nucleotide'''
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        info = {"Id":id,"Features":{"$elemMatch":{"Type":"CDS"}}}
        finder_data = collection_data_nucleotide.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {id}: CDS Not Found")
            return None
        collection_data = db["protein_data"]
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
        finder_data = collection_data.find_one(info)
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
        collection_data = db["taxonomy_data"]
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
        finder_data = collection_data.find_one(info)
        if not finder_data:
            PrintWarning(5).stdout(f"Error searching {taxon_id}: Taxonomy Not Found In Database...","\n","\t\tSearching on NCBI...")
            dataFind = self.ncbiSearchTaxon(taxon_id,"taxonomy")
            for data in dataFind:
                collection_data.insert_one(data)
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
                # finder_data = collection_data.find_one(info)
        return finder_data


    def ncbiSearch(self,id:str,database:str) -> tuple:
        # NB: in futuro la composizione del link potrebbe cambiare nella sua struttura, in base alla gestione interna di NCBI.
        # TO-DO: se questo metodo non ritorna i dati correttamente andrà aggiornata la procedura di reperimento dei dati.
        handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        p1,p2 = self.parse.parseGene(record)
        return p1,p2


    def ncbiSearchTaxon(self,id:str,database:str) ->tuple:
        handle = Entrez.efetch(db=database, id=id, retmode="xml")
        read = Entrez.read(handle)
        return read
    

    def ncbiSearchGenome(self,id:str,database:str) ->tuple:
        handle = Entrez.efetch(db=database, id=id, retmode="xml")
        read = Entrez.read(handle)
        #print(read)
        # handle = Entrez.efetch(db=database, id=id,rettype="gb", retmode="text")
        # record = SeqIO.read(handle, "genbank")
        return read


    def ritornodicose(self) -> None:
        finder = {"Id":"","Features":{"$elemMatch":{"Type":"CDS"}}}
        db = self.client["Biologia"]
        collection_data_nucleotide = db["nucleotide_data"]
        data = collection_data_nucleotide.count_documents(finder)
        print(data)


    def confronto(self): #CAMBIA NOME
        # for key in self.dataSource:
        #     for id in self.dataSource[key]:
                # data = self.proteinFind(id)
                # if not data:
                #     PrintWarning(5).stdout(f"ID:{id}")
                #     #return None
                # else:
                #     PrintWarning(3).stdout(f"Protein ID:{id}")
                
                # data = self.taxonFind(id)
                # if not data:
                #     PrintWarning(5).stdout(f"ID:{id}")
                #     #return None
                # else:
                #     PrintWarning(3).stdout(f"Taxon ID:{id}")     
        data = self.genomeFind(4)
                # if not data:
                #     PrintWarning(5).stdout(f"ID:{id}")
                #     #return None
                # else:
                #     PrintWarning(3).stdout(f"Genome ID:{id}")



#######################
######---MAIN---#######
#######################


def main(args:dict) -> None:
    v = False
    if args["verbose"]:
        v = True
    type = "nucleotide"
    if args["protein"]:
        type = "protein"
    email = None
    if args["email"]:
        email = args["email"]
    d = Database(verbose=v,type=type,email=email)
    if args["nosql_mongo"] and args["sqlite3"]:
        return PrintWarning(5).stdout("Can't select both mongo-db and sql for storing")
    if args["file"]:
        p = Parsing(args["file"])
        data,conv = p.parsing_gene()
        if args["json"]:
            p.save_data((data,conv))
        if args["nosql_mongo"]:
            d.save_on_mongo((data,conv))
        if args["sqlite3"]:
            d.save_on_sql((data,conv))
    elif args["find"]:
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
        if args["list"]:
            PrintWarning(4).stdout(d.listSource())
        elif args["alignment"]:
            fileIn,fileOut = args["alignment"]
            d.alignment(fileIn=fileIn,fileOut=fileOut)
        if args["fasta"]:
            #print(args["fasta"])
            d.toFasta(args["fasta"])
    #d.isMicorAlgae()
    return



#Costrutto che permette di passare da riga di comando l'opzione desiderata per l'inserimento dei dati in un tipo di DB
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='Biologia Database Parsing',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('-v', '--verbose',
                        action='store_true')
    parser.add_argument('-f', '--file')
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
    parser.add_argument('-l','--list',
                        action='store_true')
    parser.add_argument('-p','--protein',
                        action='store_true')
    parser.add_argument('--alignment',nargs=2,
                        metavar=('fromfile', 'tofile'),
                        help='Align all organisms in __fromfile__ and save it in __tofile__',
                        )
    parser.add_argument('--email')
    parser.add_argument('--fasta')
    args = vars(parser.parse_args())
    main(args)
    



