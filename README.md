# MicroAlgae Parsing & Searching

Credits
[`federico-rosatelli`](https://github.com/federico-rosatelli) [`Mat`](https://github.com/AxnNxs) [`Loriv3`](https://github.com/Loriv3) [`Samsey`](https://github.com/Samseys) [`Calli`](https://github.com/BboyCaligola)


# Requisiti Hardware:
Questo progetto ospiterà una grande quantità di dati, dovuta principalmente a n sequenze FASTA 
e/o FASTQ di una vasta varietà di organismi. Pertanto si consiglia di disporre di:

1. PC ad-hoc connesso costantemente a internet, dove verrà trasferito il DB e il Server virtuale
2. Archiviazione da 1+ TB per ospitare tutti i dati



# Requisiti Software:

Anaconda:             "https://www.anaconda.com"

Biopython:            "https://biopython.org/wiki/Documentation"

Hashlib:              "https://docs.python.org/3/library/hashlib.html"

Argparse:             "https://docs.python.org/3/library/argparse.html"

Pymongo:              "https://pymongo.readthedocs.io/en/stable/"

MongoDB CE "Jammy":   "https://www.mongodb.com/docs/manual/tutorial/install-mongodb-on-ubuntu/"

NVM:                  "https://www.linode.com/docs/guides/how-to-install-use-node-version-manager-nvm/"

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

EBI:                  "https://www.ebi.ac.uk/"           DB europeo




# Altri link e Repository:

Drive:                "https://drive.google.com/drive/folders/19RXRHEb-7-O9gaUjXz5ho-Q2_HsbKlEW"

Github:               "https://github.com/federico-rosatelli/biologia"



# Struttura del Progetto:
L'insieme di questi 4 moduli hanno dato vita al progetto qui in opera:

- Database: MongoDB
- Parser:   Bioparse.py
- Backend:  GoServer
- Frontend: VueJS



# Guida al primo utilizzo
Per la prima installazione su un qualsiasi PC, seguire i seguenti passaggi:

[TODO] Docker image [TODO] 
- Scaricare/clonare il progetto tramite git alla repository "https://github.com/federico-rosatelli/biologia"
- Installare tutti i Requisiti Software, settando, se necessario, le variabili di ambiente a livello globale
- Inizializzare MongoDB. Dare i seguenti comandi da terminale per evitare problemi di permessi di scrittura:
  
    ``sudo chown -R mongodb:mongodb /var/lib/mongodb``          Setting di permessi
        
    ``sudo systemctl start mongod``                             Avviare MongoDB
        
    ``sudo systemctl enable mongod.service``                    Per rendere automatico l'avvio
        
- [DEPRECATED] ~~runnare lo script genbank.py, specificando un file .gbk sorgente e selezionando il tipo di database desiderato (-n per MongoDB, -s per SQLite3)~~ [COMMENT]: # La maggior parte dei dati presenti nel Database del progetto sono stati scaricati manualmente a causa delle restrizioni sui
download da parte delle piattaforme dei Riferimenti Ufficiali. Tali metodi sono elencati in parsingMethods.py. In seguito, si é proceduto a eliminare le voci che non erano interesse del progetto,come per esempio il ramo tassonomico dei Bacteria, intervenendo a livello di DB interno con i comandi forniti da MongoDB, effettuando query dei dati tramite regex e ranking tassonomico.

- MongoDB conserva i dati implicitamente in memoria. Effettuare il download dei dump del DB presenti al link [INSERIRELINK!], decomprimere gli archivi e lanciare, per ogni cartella decompressa:

    ``mongorestore --db=Biologia <path_to_folder_extracted>``
    
    Da questo momento é possibile consultare il DB da terminale.


- Scaricare i dati mancanti di Genome insieme alle analisi prodotte dal software BUSCO dal link [INSERIRELINK!]
- Scaricare i dati mancanti di TSA ed SRA eseguendo lo script [INSERIRELINK!]



- [DEPRECATED] ~~Aprire il server, digitando su console all'interno della cartella nodeServer:~~
    ~~``npm install --save express`` per il primo avvio;~~
    ~~``node server.js``~~ [COMMENT]: # Il backend è stato sostituito da una versione più completa e consistente scritta in GoLang.

- Inizializzare il server, digitando, all'interno della cartella GOServer:

    ``go build ./cmd/webapi/``                                  Per avviare il backend

    ``run npm dev``                                             Per avviare il frontend







- Ora é possibile consultare il DB da terminale. Si riportano alcuni comandi utili allo scopo:
    
    ``mongosh``                                                 Per utilizzare MongoDB
    ``show dbs``                                                Per visualizzare i DB attivi
    ``use Biologia``                                            Per entrare nel DB creato durante il mongorestore
    ``show collections``                                        Per visualizzare le collezioni di elementi del DB
    ``db.<collection_name>.find()``                             Query che ritorna tutte le occorrenze se vuota
    ``db.<collection_name>.findOne()``                          Query che ritorna la prima occorrenza dell'input



- Ora é possibile consultare il DB da interfaccia web al seguente indirizzo:

    ``http://localhost:5173``

 
 
 
# Le collections trattate:
- nucleotide_data    contiene tutti i dati delle microalghe che sono riscontrabili su NCBI alla voce Nucleotide;
- nucleotide_basic   contiene i dati, allegeriti soltanto a nome e NCBI_ID, per questioni di performance quando si effettuano query di conteggio;
- taxonomy_data      contiene tutti i dati delle microalghe che sono riscontrabili su NCBI alla sezione Taxonomy;
- taxonomy_tree      contiene i link di Lineage per singola specie.



# Dati esterni al DB:
Per motivi di performance e trattamento di grandi moli di dati, si é deciso di:

- Per la collection nucleotide_data, le specie a cui fanno riferimento più di 10.000 entry su NCBI sono state tagliate.
Ciò perché NCBI limita fortemente la velocità di download e a livello prettamente numerico il peso del DB locale sarebbe
salito di un ordine di grandezza. Le voci sono incomplete se al numero di occorrenze viene accodato un "+" nella view
di nucleotide consultabile da interfaccia web. Le voci mancanti possono essere inserite, per singola specie, utilizzando
le funzioni efetch ed esearch di biopython e incrementando il parametro retstart al precedente retmax (il massimo è appunto
10.000). Un ulteriore requisito é essere registrati su NCBI e inserire le proprie credenziali per fare le consultazioni:

    ``Entrez.email = <email_registered_on_NCBI>``
    ``Entrez.api_key = "cc030996838fc52dd1a2653fad76bf5fe408"``



# La rappresentazione FASTA e FASTQ
FASTA conserva soltanto la Sequenza di nucleotidi o amminoacidi, codificando ogni gene in singole lettere
per indice di posizione. Nella rappresentazione in Genbank, troviamo tale dato nel file JSON che salviamo
in locale, alla voce translation per ogni Coding Sequence sotto ogni Specie, secondo la seguente gerarchia:
 
- SPECIE
    - FEATURES
        - CDS
            - /translation="LSLAVGTTITLASYHWLL[...]""

FASTQ è un "quality score" che associa alla sequenza, per ogni indice di posizione, un valore 
qualitativo codificato in ASCII. Un esempio a seguire:

- @SRR64[...]       Name Sequence
- CCTCGTCTA[...]    DNA Sequence
- +SRR64[...]       Quality address
- BBBBBFFFF[...]    Quality Score



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
