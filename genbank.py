from Bio import SeqIO
import json

NAME = "file.gbk"
def parsing(file_name):

    record = SeqIO.parse(file_name, "genbank")
    gene_array = []
    for gene in record:
        json_gene = {}
        json_gene["Name"] = gene.name
        json_gene["Id"] = gene.id
        json_gene["Seq"] = str(gene.seq)
        json_gene["Description"] = gene.description
        json_gene["Features"] = {}
        for feature in gene.features:
            json_gene["Features"][feature.type] = {}
            for qual in feature.qualifiers:
                json_gene["Features"][feature.type][qual] = "".join(feature.qualifiers.get(qual))
            json_gene["Features"][feature.type]["Location"] = [int(feature.location.start),int(feature.location.end)]
        gene_array.append(json_gene)
    return gene_array

data = parsing(NAME)

with open("datastruct.json","w") as js:
    json.dump(data,js,indent=4)
