from Bio import SeqIO
import json

def parsing(file_name):

    record = SeqIO.parse(file_name, "genbank")
    gene_array = []
    for gene in record:
        json_gene = {}
        json_gene["Name"] = gene.name
        json_gene["Id"] = gene.id
        json_gene["Seq"] = gene.seq
        gene_array.append(json_gene)

print(gene_array)