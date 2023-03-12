CREATE TABLE IF NOT EXISTS Biologia(
    id INT PRIMARY KEY AUTO_INCREMENT,
    Name VARCHAR(30),
    Id VARCHAR(30),
    Seq_Index VARCHAR(60),
    Description VARCHAR(128),
    Features VARCHAR(30)
);

CREATE TABLE IF NOT EXISTS Source(
    id INT PRIMARY KEY AUTO_INCREMENT,
    Organism VARCHAR(60),
    Organelle VARCHAR(30), 
    Mol_Type VARCHAR(30),
    Isolate VARCHAR(30),
    Db_Xref VARCHAR(30),
    LocationStart INT,
    LocationEnd INT
);

CREATE TABLE IF NOT EXISTS HexSeq(
    id INT PRIMARY KEY AUTO_INCREMENT,
    Hex VARCHAR(60),
    Raw TEXT
);
