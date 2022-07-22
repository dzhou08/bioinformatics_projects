from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio import Entrez

Entrez.email = 'danielzh08@gmail.com'

# open genes.txt file containing one gene per line
with open('genes.txt') as gene_file:
    for gene_type in gene_file:
        
        # serach genbank for gene ID
        handle = Entrez.esearch(db = "nuccore", term = gene_type)
        record = Entrez.read(handle)
        gi_list = record["IdList"]
        gi_str = ",".join(gi_list)

        print(gi_str)

        # get genebank file in gb format
        handle = Entrez.efetch(db = "nuccore", id = gi_str, rettype = "gb", retmode = "text")
        records = SeqIO.parse(handle, "gb")

        for record in records:
            #record = SeqIO.read("sequence.gb", "gb")
            for f in record.features:

                # find all CDS types 
                # do nucleotide translation and comparison
                if f.type == "CDS":
                    gene_seq = Seq(record.seq)
                    gene_sub_seq = f.location.extract(gene_seq)
                    calculated_translation_string = gene_sub_seq.translate()
                    # + "*" because extract() function always returns a string ending with "*"
                    translation_string = f.qualifiers["translation"][0] + "*"
                    if calculated_translation_string == translation_string:
                        print(f'{f.qualifiers["protein_id"][0]} Translation results match')
                    else:
                        print(f'{f.qualifiers["protein_id"][0]} Translation results do not match')