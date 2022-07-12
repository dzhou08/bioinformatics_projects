from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature,FeatureLocation
record = SeqIO.read("sequence.gb", "gb")
for f in record.features:
    if f.type == "CDS":
        gene_seq = Seq(record.seq)
        gene_sub_seq= f.location.extract(gene_seq)
        calculated_translation_string = gene_sub_seq.translate()
        for item in f.qualifiers:
            if item == "translation":
                # + "*" because extract() function always returns a string ending with "*"
                translation_string = f.qualifiers[item][0] + "*"
        if calculated_translation_string == translation_string:
            print("vamos")
        break