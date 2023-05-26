from Database import Database
from NCBISearch import NCBIFinder

class BioPyParse:
    def __init__(self,verbose:bool|None=False) -> None:
        self.database = None
        self.ncbi = NCBIFinder
        self.verbose = verbose
        pass

    def __str__(self) -> str:
        return "Biologia Computazionale Team Sapienza 2023"
    
    def newDatabase(self,databaseName:str,clientIp:str="localhost",clientPort:int=27017) -> None:
        self.database = Database.New(clientIp,clientPort)
        self.database.addDatabaseName(databaseName)
        return
    
    def showCollections(self):
        return self.database.showCollections()
    
    def importTaxonFromList(self,listTaxonNames:list[str],collectionName:str|None = "taxonomy_data",next_level:bool|None=True) -> None:
        '''Function that, given a list of Taxonomy ID, add the relative specie to collection_data branch in Biologia DB,
            if it is not already present in list'''
        str_next_leve = "[next level]" if next_level else ""
        for name in listTaxonNames:
            try:
                taxon = self.ncbi.ncbiSearchTaxon(f"{name}{str_next_leve}")
                if self.verbose:
                    print(f"Length of {name}: {len(taxon)}")
                for tax in taxon:
                    if (not self.database.findOneCollection(collectionName,{"TaxId":tax["TaxId"]}) and tax["Rank"] == "species"):
                        if self.verbose:
                            print(f"Scientific Name: {tax['ScientificName']}")
                        self.database.insertOneCollection(collectionName,tax)
            except Exception as e:
                print(f"Raise Exception: {e} searching: {name}")
    
    def generateTaxonomyTree(self,collectionName:str|None = "taxonomy_tree",taxonomyCollection:str|None = "taxonomy_data") -> None:
        '''Function that retrieves datas from collection_data of MongoDB and
            create a new collection, taxonomy_tree, that contain the Lineage for
            each specie'''
        results = self.database.findManyCollection(taxonomyCollection,{})
        iter = 0
        for res in results:
            iter += 1
            if self.verbose:
                print(f"Organism: {res['ScientificName']}, index: {iter}")
            for i in range(1, len(res["LineageEx"])):
                LineageEx = res["LineageEx"][i]
                if not self.database.findOneCollection(collectionName,{"TaxId":LineageEx["TaxId"]}):
                    dataPush = {
                        "TaxId":LineageEx["TaxId"],
                        "Rank":LineageEx["Rank"],
                        "ScientificName":LineageEx["ScientificName"],
                        "SubClasses":[]
                    }
                    self.database.insertOneCollection(collectionName,dataPush)
                if (i == 1):
                    continue
                LineageExPrev = res["LineageEx"][i-1]

                if self.database.findOneCollection(collectionName,{"TaxId":LineageExPrev["TaxId"], "SubClasses":{"$elemMatch":{"TaxId":LineageEx["TaxId"]}}}):
                    continue
                newDataPush = {
                        "TaxId":LineageEx["TaxId"],
                        "Rank":LineageEx["Rank"],
                        "ScientificName":LineageEx["ScientificName"],
                        }
                self.database.updateOneCollection(collectionName,{"TaxId":LineageExPrev["TaxId"]},{"$push":{"SubClasses":newDataPush}})
            if not self.database.findOneCollection(collectionName,{"TaxId":res["ParentTaxId"], "SubClasses":{"$elemMatch":{"TaxId":res["TaxId"]}}}):
                newDataPush = {
                        "TaxId":res["TaxId"],
                        "Rank":res["Rank"],
                        "ScientificName":res["ScientificName"],
                        }
                self.database.updateOneCollection(collectionName,{"TaxId":res["ParentTaxId"]},{"$push":{"SubClasses":newDataPush}})
    
    def importNucleotide(self,collectionName:str|None = "nucleotide_data",taxonomyCollection:str|None = "taxonomy_data"):
        all_taxons_id = self.database.findManyCollection(taxonomyCollection,{},{"TaxId":1})
        for taxTerm in all_taxons_id:
            taxId = taxTerm["TaxId"]
            try:
                nucleos = self.ncbi.ncbiSearchNucleo(f"txid{taxId}")
                if self.verbose:
                    print(f"Length of {taxId}: {len(nucleos)}")
                for nucleo in nucleos:
                    nucleo.pop("GBSeq_sequence",None)
                    self.database.insertOneCollection(collectionName,nucleo)
            except Exception as e:
                print(f"Raise Exception: {e} searching: {taxId}")
    
    def importProtein(self,collectionName:str|None = "protein_data",taxonomyCollection:str|None = "taxonomy_data"):
        all_taxons_id = self.database.findManyCollection(taxonomyCollection,{},{"TaxId":1})
        for taxTerm in all_taxons_id:
            taxId = taxTerm["TaxId"]
            try:
                nucleos = self.ncbi.ncbiSearchNucleo(f"txid{taxId}")
                if self.verbose:
                    print(f"Length of {taxId}: {len(nucleos)}")
                for nucleo in nucleos:
                    nucleo.pop("GBSeq_sequence",None)
                    self.database.insertOneCollection(collectionName,nucleo)
            except Exception as e:
                print(f"Raise Exception: {e} searching: {taxId}")
    
    def createTable(self,collectionName:str|None = "table",taxonomyCollection:str|None = "taxonomy_data",nucleotideCollection:str|None = "nucleotide_data",proteinCollection:str|None = "protein_data"):
        filter = {"GBSeq_feature-table.GBFeature_key":"CDS","GBSeq_feature-table.GBFeature_quals.GBQualifier_name":"product"}
        proj = {"GBSeq_feature-table":1,"GBSeq_organism":1}
        dataResult = self.database.findManyCollection(nucleotideCollection,filter, proj)
        dataProduct = {}
        for dataN in dataResult:
            if dataN["GBSeq_organism"] not in dataProduct:
                dataProduct[dataN["GBSeq_organism"]] = []
            for data in dataN["GBSeq_feature-table"]:
                if data["GBFeature_key"] == "CDS" or data["GBFeature_key"] == "rRNA":
                    for prod in data["GBFeature_quals"]:
                        if prod["GBQualifier_name"] == "product":
                            isFind = False
                            for k in range(len(dataProduct[dataN["GBSeq_organism"]])):
                                if dataProduct[dataN["GBSeq_organism"]][k]["ProductName"] == prod["GBQualifier_value"]:
                                    dataProduct[dataN["GBSeq_organism"]][k]["QtyProduct"] = dataProduct[dataN["GBSeq_organism"]][k]["QtyProduct"] + 1
                                    isFind = True
                            if not isFind:
                                dataProduct[dataN["GBSeq_organism"]].append({
                                    "ProductName" : prod["GBQualifier_value"],
                                    "QtyProduct" : 1
                                    })
        #Chi ha levato Country?
        dataRank = self.database.findManyCollection(taxonomyCollection,{},{"ScientificName":1,"TaxId":1,"_id":0})
        control = 0
        for data in dataRank:
            id = data["TaxId"]
            if self.verbose:
                print(f"Organism: {data['ScientificName']}, Index: {control}")
            finder = {"GBSeq_feature-table.GBFeature_quals.GBQualifier_value":f"taxon:{id}"}
            
            totNuclBase = self.database.findManyCollection(nucleotideCollection,finder,{"GBSeq_locus":1, "_id":0})
            totNuclBase = list(totNuclBase)
            totNuclCompl = self.database.findManyCollection(nucleotideCollection,finder)
            totNuclCompl = list(totNuclCompl)

            totProtBase = self.database.findManyCollection(proteinCollection,finder,{"GBSeq_locus":1, "_id":0})
            totProtBase = list(totProtBase)
            totProtCompl = self.database.findManyCollection(proteinCollection,finder)
            totProtCompl = list(totProtCompl)
            
            inserterBase = {
                "ScientificName":data["ScientificName"],
                "TaxId":data["TaxId"],
                "Nucleotides":totNuclBase,
                "Proteins":totProtBase,
                "Products": dataProduct[data["ScientificName"]],
                "Country": []
            }
            inserterCompl = {
                "ScientificName":data["ScientificName"],
                "TaxId":data["TaxId"],
                "Nucleotides":totNuclCompl,
                "Proteins":totProtCompl,
                "Products": dataProduct[data["ScientificName"]],
                "Country": []
            }
            self.database.insertOneCollection(f"{collectionName}_basic",inserterBase)
            self.database.insertOneCollection(f"{collectionName}_basic",inserterCompl)
            control += 1   





def findSpeciesFromFile(filepath:str,splitColumns:str|None=";",indexSearch:int|None=0,fromLine:int|None=0) -> list[str]:
    '''Function that, given a file path in input, scan for species and datas.\n
    Parameters:
        - filepath: string with absolute or relative path of the file
        - splitColumns: string for splitting the columns
        - indexSearch: integer which column find the species
        - fromLine: integer which line start the search
    '''
    names:list[str] = []
    rd = open(filepath).readlines()
    for i in range(fromLine,len(rd)):
        rd[i] = rd[i].strip()
        data = rd[i].split(splitColumns)
        name = data[indexSearch].split(" ")[0]
        if "_" in name:
            name = name.split("_")[0]
        names.append(name)
    return names
