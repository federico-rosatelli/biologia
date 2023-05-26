from pymongo import MongoClient





class Database:
    def __init__(self,clientIp:str,clientPort:int) -> None:
        self.clientIp = clientIp
        self.clientPort = clientPort
        self.client = MongoClient(self.clientIp,self.clientPort)
        self.databaseName = ""
        self.db = None
        pass
    def __str__(self) -> str:
        return f'{self.clientIp}:{self.clientPort}'
    
    def addDatabaseName(self,name:str) -> None:
        self.databaseName = name
        self.db = self.client[self.databaseName]
    
    def showCollections(self) -> list[str]:
        return self.client[self.databaseName].list_collection_names()
    
    def findOneCollection(self,collection:str,filter:dict,proj:dict|None={}) -> dict | None:
        return self.db[collection].find_one(filter,proj)
    
    def findManyCollection(self,collection:str,filter:dict,proj:dict|None={}) -> dict | None:
        return list(self.db[collection].find(filter,proj))
    
    def insertOneCollection(self,collection:str,inserter:dict) -> None:
        self.db[collection].insert_one(inserter)
    
    def updateOneCollection(self,collection:str,filter:dict,update:dict) -> None:
        self.db[collection].update_one(filter,update)

def New(clientIp:str,clientPort:int) -> Database:
    return Database(clientIp,clientPort)