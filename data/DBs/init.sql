CREATE TABLE IF NOT EXISTS Biologia(
    Name VARCHAR(30),
    Id VARCHAR(30),
    Seq_Hex VARCHAR(60),
    Description VARCHAR(128),
    Features VARCHAR(30)
);

CREATE TABLE IF NOT EXISTS Source(
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    Organism VARCHAR(60),
    Organelle VARCHAR(30), 
    Mol_Type VARCHAR(30),
    Isolate VARCHAR(30),
    Db_Xref VARCHAR(30),
    LocationStart INT,
    LocationEnd INT
);

CREATE TABLE IF NOT EXISTS HexSeq(
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    Hex VARCHAR(60),
    Raw TEXT
);
