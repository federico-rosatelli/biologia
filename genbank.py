from Bio import SeqIO
import json
import hashlib
import argparse
from pymongo import MongoClient
import sqlite3


# Credits
# Authors: federico-rosatelli (Federico Rosatelli), AxnNxs (Mattia Di Gesaro)


# Riferimenti ufficiali
# biopython:    "https://biopython.org/wiki/Documentation"
# hashlib:      "https://docs.python.org/3/library/hashlib.html"
# argparse:     "https://docs.python.org/3/library/argparse.html"
# pymongo       "https://pymongo.readthedocs.io/en/stable/"


class Error:
    # applicazione della classe di bcolors
    def __init__(self,type:int) -> None:
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
        self.error = None

    
    def newError(self,error) ->None:
        self.error = error
        print(self.color + self.error + bcolors.ENDC)
    
    def __str__(self) -> str:
        return self.error

CLUSTER = "localhost:27017"     # apertura porta di default sul localhost che esegue il codice

class Parsing(object):
    def __init__(self,file_name:str) -> None:
        # funzione di inizializzazione del dato in input
        self.record = SeqIO.parse(file_name, "genbank")

    def parsing_gene(self) -> list:
        # Ritorna due liste: 
        # - nella prima vi sono i dati analizzati: per ogni file.gbk dato in input ci sono  di geni,
        #   e per ogni gene avremo un dizionario dei campi di interesse, come il Name, l'ID, e così via;
        # - nella seconda vi sono i dati trattati riferiti alla Sequenza. Salviamo infatti sia la sequenza
        #   di basi azotate "as is", sia la versione codificata con SHA256.  
        gene_array = []
        hex_array = []
        for gene in self.record:
            json_gene = {}
            json_convert = {}
            json_gene["Name"] = gene.name
            json_gene["Id"] = gene.id
            sh = self.sha_index(str(gene.seq))
            json_gene["Seq_Hex"] = sh["hex"]  #str(gene.seq)-> SHA256 funzione
            #json_gene["Seq_Raw"] = sh["seq"]  #str(gene.seq)-> SHA256 funzione
            json_convert["Seq_Hex"] = sh["hex"]  #str(gene.seq)-> SHA256 funzione
            json_convert["Seq_Raw"] = sh["seq"]  #str(gene.seq)-> SHA256 funzione
            json_gene["Description"] = gene.description
            json_gene["Features"] = {}
            for feature in gene.features:
                json_gene["Features"][feature.type] = {}
                for qual in feature.qualifiers:
                    json_gene["Features"][feature.type][qual] = feature.qualifiers.get(qual)[0]
                json_gene["Features"][feature.type]["Location"] = (int(feature.location.start),int(feature.location.end))
            gene_array.append(json_gene)
            hex_array.append(json_convert)
        return gene_array,hex_array
    
    def sha_index(self,seq:str) -> dict:
        # ritorna, in un dizionario, il codice hash di una sequenza di basi azotate (ACGTU le possibili unità)
        return {
            'hex':hashlib.sha256(seq.encode()).hexdigest(),
            'seq':seq
        }
    
    def save_data(self, json_array:tuple) -> None:
        # crea un file json contentene i dati parsati - salvare in DB?
        print("AIIII")
        with open("datastruct.json","w") as js:
            json.dump({
                        "struct":json_array[0],
                        "hex":json_array[1]
                        },js,indent=4)




class Database:
    def __init__(self) -> None:
        pass
    
    # def save_on_mongo(self)
    def save_on_mongo(self,tuple_of_array:tuple) -> None:
        ip = CLUSTER.split(":")[0]
        port = int(CLUSTER.split(":")[1])
        client = MongoClient(f'{ip}', port)
        db = client["Biologia"]  # attenzione a quando si richiama il DB da riga di comando: case sensitive
        collection_data = db["genetic_data"]
        collection_convert = db["hex_to_seq"]
        for i in range(len(tuple_of_array[0])):
            filter = {"Name":tuple_of_array[0][i]["Name"]}
            check_name = collection_data.find_one(filter)
            if not check_name:
                collection_data.insert_one(tuple_of_array[0][i])
            filter = {"Seq_Hex":tuple_of_array[1][i]["Seq_Hex"]}
            check_hex = collection_convert.find_one(filter)
            if not check_hex:
                collection_convert.insert_one(tuple_of_array[1][i])
        client.close()
        
    # def save_on_sql(self)
    def save_on_sql(self,tuple_of_array:tuple) -> Error:
        # filename to form database
        
        file = "Biologia.db"
        try:
            conn = sqlite3.connect(file)
            print("Database Biologia.db formed.")
            cursor = conn.cursor()
            scr = open("init.sql").read()
            cursor.executescript(scr)
            list_of_tuples = [(tuple_of_array[0][i]["Name"],tuple_of_array[0][i]["Id"],tuple_of_array[0][i]["Seq_Hex"],tuple_of_array[0][i]["Description"],i) for i in range(len(tuple_of_array))]
            #cursor.executemany("INSERT INTO movie VALUES(?, ?, ?)", data)
            table ="""INSERT INTO Biologia VALUES(?,?,?,?,?)"""
            # cursor.execute(table)
            # cursor.execute()
            cursor.executemany(table,list_of_tuples)
            print("ASASAS")
        except Exception as e:
            conn.close()
            return Error(4).newError(f"Error Database Biologia.db: {e}")
        conn.commit()
        if conn:
            conn.close()
            return Error(5).newError(f"Error Writing Database Biologia.db: {e}")
        return None


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





def main(args:dict) -> Error:
    if args["mongodb"] and args["sql"]:
        return Error(5).newError("Can't select both mongo-db and sql for storing")
    if args["file"] == None:
        return Error(5).newError("File must not be empty")
    p = Parsing(args["file"])
    data,conv = p.parsing_gene()
    if args["json"]:
        p.save_data((data,conv))
    d = Database()
    if args["mongodb"]:
        d.save_on_mongo((data,conv))
    if args["sql"]:
        d.save_on_sql((data,conv))
    
    return



    

if __name__ == "__main__":
    # costrutto che permette di passare da riga di comando l'opzione desiderata per l'inserimento dei dati 
    # in un tipo di DB
    parser = argparse.ArgumentParser(
                    prog='Biologia Database Parsing',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument('-m', '--mongodb',
                        action='store_true')
    parser.add_argument('-s', '--sql',
                        action='store_true')
    parser.add_argument('-j', '--json',
                        action='store_true')
    parser.add_argument('-f', '--file')
    args = vars(parser.parse_args())
    main(args)
    



