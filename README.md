# MicroAlgae Parsing & Searching

Credits
[`federico-rosatelli`](https://github.com/federico-rosatelli) [`Mat`](https://github.com/AxnNxs) [`Loriv3`](https://github.com/Loriv3) [`Samsey`](https://github.com/Samseys)



# Requisiti:

Anaconda:             "https://www.anaconda.com"

Biopython:            "https://biopython.org/wiki/Documentation"

Hashlib:              "https://docs.python.org/3/library/hashlib.html"

Argparse:             "https://docs.python.org/3/library/argparse.html"

Pymongo:              "https://pymongo.readthedocs.io/en/stable/"

MongoDB CE "Jammy":   "https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/"

Nodejs:               "https://github.com/nodejs/help/wiki/Installation"

NPM:                  "https://docs.npmjs.com/downloading-and-installing-node-js-and-npm"

Nodejs Express:       "http://expressjs.com/"

GoLang:               "https://go.dev/doc/install"

VUEjs:                "https://cli.vuejs.org/guide/installation.html"

Docker:               "https://docs.docker.com/engine/install/ubuntu/"

Docker compose:       "https://docs.docker.com/desktop/install/linux-install/" (Integrato in Desktop)



# Fonti ufficiali:

Shigen:               "https://shigen.nig.ac.jp"         DB giapponese

NCBI:                 "https://www.ncbi.nlm.nih.gov/"    DB americano

La maggior parte dei dati presenti localmente nel Database interno
sono stati scaricati manualmente a causa delle restrizioni sui
download.



# Altri link:

Drive:                "https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW"

Github:               "https://github.com/federico-rosatelli/biologia"




# Guida al primo utilizzo
Per la prima installazione su un qualsiasi PC, seguire i seguenti passaggi (si raccomanda Ubuntu):
 
   - Scaricare il progetto tramite git, o copiare il progetto all'indirizzo inserito tra "Altri link" sopra questo paragrafo.
   - Installare tutti i requisiti, settando le variabili di ambiente a livello globale
   - Inizializzare MongoDB:
        Dare i seguenti comandi da terminale:
  
        ``sudo chown -R mongodb:mongodb /var/lib/mongodb``

        ``sudo chown mongodb:mongodb /tmp/mongodb-27017.sock``

        dopodiché:
        
        ``sudo systemctl start mongod``
        
        per rendere automatico l'avvio di mongod in fase di accensione dell'SO:
        
        ``sudo systemctl enable mongod.service``
        
   - runnare lo script genbank.py, specificando un file .gbk sorgente e selezionando il tipo di database desiderato (-n per MongoDB, -s per SQLite3)
   - attendere il termine del salvataggio nel database locale
   - aprire il server, digitando su console all'interno della cartella nodeServer:

       ``node server.js``
   - dopodiché aprire il browser e cercare il seguente indirizzo:

       localhost:3000

   è ora possibile effettuare le query di interesse, ma va inizializzato il database locale.

   Per avere le specie con i dati genomici (versione MongoDB):

``python3 genbank.py --file "nomefile".gbk -n``

   Per la taxonomy:

``python3 genbank.py --find -m --email`` (inserire email per accesso su ncbi, ammesso che si abbia accesso)



NOTA: a ogni riavvio del sistema, è necessario riattivare mongodb (``systemctl start mongod``) e il server (``node server.js``)






# Concetti utili:
MONGODB conserva i dati implicitamente in una memoria virtuale. Per trasportarli da un sistema
a un altro è necessario utilizzare il comando mongodump per generare una cartella contenente
il db di interesse e sudo mongorestore sul file bson generato dal mongodump una volta importata la cartella generata.
 


# La rappresentazione FASTA e FASTQ
FASTA conserva soltanto la Sequenza di nucleotidi o amminoacidi, codificando ogni gene in singole lettere
per indice di posizione. Nella rappresentazione in Genbank, troviamo tale dato nel file JSON che salviamo
in locale, alla voce translation per ogni Coding Sequence sotto ogni Specie, secondo la seguente gerarchia:
 
# SPECIE
##   FEATURES
###       CDS
####           /translation="LSLAVGTTITLASYHWLL[...]""

FASTQ è un "quality score" che associa alla sequenza, per ogni indice di posizione, un valore 
qualitativo codificato in ASCII. Un esempio a seguire:
@SRR64[...]       Name Sequence
CCTCGTCTA[...]    DNA Sequence
+SRR64[...]       Quality address
BBBBBFFFF[...]    Quality Score



# Sequenziamento genetico
E' il processo di determinazione dell'ordine dei nucleotidi (Adenina, Citosina, Guanina e Timina) che 
costituiscono il frammento di DNA in analisi. Le tecniche principali di Squenziamento sono
Sanger e NGS (Illumina ne é un esempio).



# Allineamento genetico
E' il processo di confronto di due o più sequenze di DNA o proteine per identificare regioni identiche
o simili per individuare eventuali relazioni funzionali, strutturali o filogenetiche.
Le tecniche di allineamento prevedono il confronto globale e locale.
Un'applicazione algoritmica di allineamento locale è data da Smith Waterman.
Un'applicazione algoritmica di allineamento globale é data da Needleman-Wunsch.



![Algae Project Struct](algaeStruct.png "Struct of the Project")


Json Files saved:

`dataStruct.json`:
```json
    {
        "struct":[
            "Name": "LOCUS",
            "Id": "VERSION",
            "Seq_Hex": "SHA256HEX",
            "Description": "DEFINITION",
            "Features": [
                {
                    "Type": "typeFeatures",
                    "Others": "others info of the feature",
                    "Location": [
                        0,
                        0
                    ]
                }
            ]
        ],
        "hex":[
            {
                "Seq_Hex":"SHA256HEX",
                "Seq_Raw":"ORIGIN"
            }
        ]
    }
    //example
    {
        "struct":[
            {
            "Name": "OQ443078",
            "Id": "OQ443078.1",
            "Seq_Hex": "de5f363dd78bc0f2a534c55783739a672ce136c67b34767a052d8ecd7e4deb5b",
            "Description": "Bacillus velezensis strain 3(JS) 16S ribosomal RNA gene, partial sequence",
            "Features": [
                {
                    "Type": "source",
                    "organism": "Bacillus velezensis",
                    "mol_type": "genomic DNA",
                    "strain": "3(JS)",
                    "db_xref": "taxon:492670",
                    "country": "India: Vikramgad",
                    "lat_lon": "19.80 N 73.09 E",
                    "collected_by": "Vir Acharya",
                    "identified_by": "Raunak Giri",
                    "Location": [
                        0,
                        936
                    ]
                },
                {
                    "Type": "rRNA",
                    "product": "16S ribosomal RNA",
                    "Location": [
                        0,
                        936
                    ]
                }
            ]
        }
        ],
        "hex":[
            {
                "Seq_Hex":"de5f363dd78bc0f2a534c55783739a672ce136c67b34767a052d8ecd7e4deb5b",
                "Seq_Raw":"ACCTGCCTGTAAGACTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTCTGAACCGCATGGTTCAGACATAAAAGGTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTGCCGTTCAAATAGGGCGGCACCTTGACGGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCCTAGAGATAGGACGTCCCCTTCGGGGGCAGAGTGACAGGTGGTGC"
            }
        ]
    }
```

`dataSource.json`
```json
    {
        "OrganismName":[
            "VERSION"
        ]
    }
    //example
    {
        "Danio rerio": [
            "NM_001327985.1",
            "NM_212741.1",
            "NM_001128675.2",
            "NM_001081554.1",
            "NM_199618.1",
            "NM_001303262.1",
            "NM_001137555.3",
            "NM_001020798.1"
        ]
    }
```

`algaeSource.json`
```json
    {
        "OrganismName":[
            "VERSION"
        ]
    }
    //example
    {
        "Bracteacoccus aggregatus": [
            "MZ090013.1",
            "MZ067570.1",
            "MH205944.1",
            "MH703758.1",
            "MH703740.1"
        ]
    }
```

`microAlgaeSource.json`
```json
    {
        "OrganismName":[
            "VERSION"
        ]
    }
    //example
    {
        "Ankistrodesmus fusiformis": [
            "OM683277.1",
            "OM683275.1"
        ]
    }
```
