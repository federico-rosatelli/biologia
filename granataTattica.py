########################################################################################
# Libraries & API
########################################################################################

from pymongo import MongoClient
from Bio import SeqIO, Entrez
from time import ctime, perf_counter
from Bio.Seq import Seq
from matplotlib import pyplot as plt
#from validate_email import validate_email
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
from xml.dom import minidom
import xml.etree.ElementTree as ET



########################################################################################
# Init and Declaration of Global variables
########################################################################################

client = MongoClient('localhost', 27017)
db = client["Biologia"]
collection_taxonomy_data = db["taxonomy_data"]
collection_nucleotide_data = db["nucleotide_data"]
collection_taxonomy_tree = db["taxonomy_tree"]
collection_table_basic = db["table_basic"]
collection_protein_data = db["protein_data"]
collection_SRA_data = db["sequences_data"]
# collection_Biosample_data = db ["sequences_data"]
# collection_SRA_data = db ["sequences_data"]
ignore_names = ["environmental samples"]                    # taxonomy ID to be avoided in collection_taxonomy_data
Entrez.api_key = "cc030996838fc52dd1a2653fad76bf5fe408"
Entrez.email = "russo.1864451@studenti.uniroma1.it"         
datasTaxon = [(2961,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2961"),
            (145388,"https://www.ncbi.nlm.nih.gov/genome/?term=txid145388"),
            (72520,"https://www.ncbi.nlm.nih.gov/genome/?term=txid72520"),
            (70448,"https://www.ncbi.nlm.nih.gov/genome/?term=txid70448"),
            (44056,"https://www.ncbi.nlm.nih.gov/genome/?term=txid44056"),
            (41880,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41880"),
            (41875,"https://www.ncbi.nlm.nih.gov/genome/?term=txid41875"),
            (36894,"https://www.ncbi.nlm.nih.gov/genome/?term=txid36894"),
            (35677,"https://www.ncbi.nlm.nih.gov/genome/?term=txid35677"),
            (3077,"https://www.ncbi.nlm.nih.gov/genome/?term=txid3077"),
            (1764295,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1764295"),
            (1650286,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1650286"),
            (1486918,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1486918"),
            (1093141,"https://www.ncbi.nlm.nih.gov/genome/?term=txid1093141"),
            (574566,"https://www.ncbi.nlm.nih.gov/genome/?term=txid574566"),
            (554065,"https://www.ncbi.nlm.nih.gov/genome/?term=txid554065"),
            (426638,"https://www.ncbi.nlm.nih.gov/genome/?term=txid426638"),
            (265525,"https://www.ncbi.nlm.nih.gov/genome/?term=txid265525"),
            (257627,"https://www.ncbi.nlm.nih.gov/genome/?term=txid257627"),
            (2009235,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2009235"),
            (2730355,"https://www.ncbi.nlm.nih.gov/genome/?term=txid2730355")]



########################################################################################
# Methods
########################################################################################

def SRAFind(name:str):                   #case 15
    '''Function that merge data found on NCBI between Genome and Nucleotide section and
    print its structore in json format'''
    handle = Entrez.esearch(db='sra', term=name + "[Organism] AND \"biomol rna\"[Properties] AND \"library layout paired\"[Properties]", retmode='text', retmax=10000)
    record = Entrez.read(handle, validate=False)
    handle.close()
    handle = Entrez.efetch(db="sra", id=record["IdList"], retmode='xml')
    #read = Entrez.read(handle, validate=False)
    handle = handle.read().decode(encoding='utf-8')
    #print(handle)
    root = ET.fromstring(handle)
    seq_list = []
    print(f"Stiamo facendo {name}")
    for child in root:
        dict_seq = funzioneRicorsiva(child)
        seq_list.append(dict_seq)
    #collection_SRA_data.insert_many(seq_list)
    with open('provariscorsiva.json', "w") as jswr:
        json.dump(seq_list, jswr, indent=4)        
    print(len(seq_list))
    return


def funzioneRicorsiva(bla) -> dict:
    response = {}
    for key, value in bla.attrib.items():
        response[key] = value
    itern = 0
    for child in bla:
        if child.tag in response:
            child.tag = child.tag + f"_{itern}"
            itern += 1
        response[child.tag] = funzioneRicorsiva(child)
        if child.text:
            response[child.tag]["value"] = child.text
    return response



########################################################################################
# Execution code unit
########################################################################################

# print("Bomba in buca!")
# count = 0
# for item in datasTaxon:
#     count += 1
#     txid = f"txid{item[0]}"
#     SRAFind(txid)
#     print(f"{count}")

#SRAFind("txid3077")



def csvWriteProf():
    print("MEEE SOOO MBRIAAACATO")
    filter = {
        "STUDY.center_name": "BioProject",
        "SAMPLE.IDENTIFIERS.EXTERNAL_ID.namespace": "BioSample",  
    }

    proj = {
        "_id":0,
        "TaxId":"$SAMPLE.SAMPLE_NAME.TAXON_ID.value",
        "BioProject":"$STUDY.alias",
        "BioSample":"$SAMPLE.IDENTIFIERS.EXTERNAL_ID.value",
        "ExperimentCode":"$EXPERIMENT.accession",
        "SRR":"$RUN_SET.RUN.accession",
        "Abstract":"$STUDY.DESCRIPTOR.STUDY_ABSTRACT.value",
        "LinkSRA":"$RUN_SET.RUN.SRAFiles"
    }
    #bson.D{{Key: "$project", Value: bson.M{"QtyNucleotides": bson.M{"$size": "$Nucleotides"}, "QtyProteins": bson.M{"$size": "$Proteins"}, "QtyProducts": bson.M{"$size": "$Products"}, "_id": 0, "ScientificName": 1, "TaxId": 1, "Genomes": 1}}}
    #{ "$project": { "employee": "$name", "salary": "$salary" }}

    aggregate = [
        {"$match": filter},
        {"$project":proj}
    ]

    all_data = collection_SRA_data.aggregate(aggregate)
    datas = []
    for data in all_data:
        links = []
        for key in data["LinkSRA"]:
            if "url" in data["LinkSRA"][key]["Alternatives"]:
                links.append(data["LinkSRA"][key]["Alternatives"]["url"])
        data["LinkSRA"] = ' '.join(links)
        datas.append(data)
    with open('bioprojects.csv',"w") as csvW:
        writer = csv.writer(csvW,delimiter="\t")
        for row in datas:
            if "Abstract" in row:
                rws = [row["TaxId"],row["BioProject"] or '', row["BioSample"] or '', row["ExperimentCode"] or '', row["SRR"] or '', row["LinkSRA"] or '', row["Abstract"] or '',]
            else:
                rws = [row["TaxId"],row["BioProject"] or '', row["BioSample"] or '', row["ExperimentCode"] or '', row["SRR"] or '', row["LinkSRA"] or '', "" or '']
            writer.writerow(rws)
    # with open('provaaggr.json', "w") as jswr:
    #     json.dump(datas, jswr, indent=4)    
    

#csvWriteProf()
#SRAFind("txid35677")


########################################################################################
# AREA SCRATCH
########################################################################################

# LIBRARY_CONSTRUCTION_PROTOCOL

	# <STUDY center_name="NCGR" alias="NCGR_marinellite" accession="SRP042159">
	# 	<IDENTIFIERS>
	# 		<PRIMARY_ID>SRP042159</PRIMARY_ID>
	# 		<EXTERNAL_ID namespace="BioProject">PRJNA248394</EXTERNAL_ID>
	# 		<SUBMITTER_ID namespace="NCGR">NCGR_marinellite</SUBMITTER_ID>
	# 	</IDENTIFIERS>
	# 	<DESCRIPTOR>
	# 		<STUDY_TITLE>Marine Microbial Eukaryote Transcriptome Sequencing Project</STUDY_TITLE>
	# 		<STUDY_TYPE existing_study_type="Other"/>
	# 		<STUDY_ABSTRACT>The Marine Microbial Eukaryote Transcriptome Sequencing Project was a collaboration between the National Center for Genome Resources, the Gordon and Betty Moore Foundation's Marine Microbiology Initiative, and the international marine microbial eukaryote research community to sequence the transcriptomes of approximately 700 samples from hundreds of diverse organisms.   Marine microbial eukaryotes are found in all major eukaryotic branches of the evolutionary tree. They perform a diverse range of functions in marine ecosystems, including photosynthesis, predation, and parasitism. Despite the great abundance of microeukaryotes in the ocean, their importance as absorbers of carbon dioxide and their critical contribution to marine food webs (among many ecological roles), the gene content of these microorganisms is only beginning to be explored because their genomes can be structurally complex and can be many gigabases in size. This program was intended to increase the research community's baseline of scientific knowledge by creating catalogs of genes that suggest how these organisms thrive in diverse marine habitats and how they influence marine ecosystems, biogeochemical cycles, and the composition of the atmosphere.  Transcriptome datasets are complex and should be approached with awareness, especially in light of the very deep Illumina sequencing used to generate these data. While many researchers attempted to provide axenic and uni-algal total RNA extracts for sequencing, it did not always turn out to be that way. The deep Illumina sequencing may have occasionally picked up very low levels of non-target RNA, and low levels of bacteria may have been present but not known to the laboratory that provided the sample.  For additional project information, details about the varied samples, and data, including transcriptome assemblies and annotations, please see: http://marinemicroeukaryotes.org/ and http://camera.calit2.net/mmetsp/.</STUDY_ABSTRACT>
	# 		<CENTER_PROJECT_NAME>Marine Microbial Eukaryote Transcriptome Sequencing Project</CENTER_PROJECT_NAME>
	# 	</DESCRIPTOR>
	# </STUDY>
# 



# values da prendere:
    # For Bioproject: save PRJNA Code
        # For Biosample: save SRS Code
            # Save Abstract through SRP: (example) Protein models were generated for the VTC4 protein of Chlamydomonas reinhardtii and putative VTC4 homologues in Chlorella vulgaris, Desmodesmus cf. armatus, and Gonium pectorale. RNA sequencing, transcriptome assembly, and protein prediction were used to identify the putative protein sequences of C. vulgaris and D. armatus, while the sequence for G. pectorale was found through a search of NCBI. Biochemical assays were used to compare polyphosphate synthesis kinetics between the algae.
            # SRX Code
                # SRA link
                # SRA code
                # Paper/Comment/Data


#<SUBMITTER_ID namespace="SUB5109492" label="1">DLS-AmphidiniumU_S2_L001_R1_001.fastq</SUBMITTER_ID>
        #  tag             attrib              attrib          text


####### example

# <foo>
#    <bar>
#       <type foobar="1"/>
#       <type foobar="2"/>
#    </bar>
# </foo>



####### Code

# import xml.etree.ElementTree as ET
# root = ET.parse('thefile.xml').getroot()

# for type_tag in root.findall('bar/type'):
#     value = type_tag.get('foobar')
#     print(value)



###### Output

# 1
# 2
