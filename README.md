# MicroAlgae Parsing & Searching

`federico-rosatelli & AxnNxs & Loriv3`

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
