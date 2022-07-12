import csv
import string

#reading gene file
gene_file = open('FGFR3_Gene_FASTA.csv')
csvreader = csv.reader(gene_file)
DNA = ""
for row in csvreader:
    for item in row:
        DNA += item
gene_file.close()
print(DNA)
pre_RNA = ''
for i in DNA:
    if i == "t":
        pre_RNA += "A"
    elif i == "g":
        pre_RNA += "C"
    elif i == "c":
        pre_RNA += "G"
    elif i == "a":
        pre_RNA += "U"
    


#reading exon file
exon_file = open('exon_positions.csv')
csvreader = csv.reader(exon_file)
exon_positions = []
for row in csvreader:
    for item in row:
        exon_positions.append(item)
exon_file.close()

#retrieve exon positions in the full DNA sequence
mature_RNA = ""
for exon_position in exon_positions:
    position = exon_position.split("..")
    mature_RNA += DNA[int(position[0]) : int(position[1])]
    print(position[0],position[1])
codon_str = ""
for n in range(int(len(mature_RNA)/3)):
    codon_string = mature_RNA[3*n : 3*(n+1)]
    codon_str += codon_dict[codon_string]

test = ''
for i in mature_RNA:
    if i == "U":
        test += "T"
    elif i == "G":
        test += "G"
    elif i == "C":
        test += "C"
    elif i == "A":
        test += "A"
print(test)
print(codon_str)
if codon_str == "MGAPACALALCVAVAIVAGASSESLGTEQRVVGRAAEVPGPEPGQQEQLVFGSGDAVELSCPPPGGGPMGPTVWVKDGTGLVPSERVLVGPQRLQVLNASHEDSGAYSCRQRLTQRVLCHFSVRVTDAPSSGDDEDGEDEAEDTGVDTGAPYWTRPERMDKKLLAVPAANTVRFRCPAAGNPTPSISWLKNGREFRGEHRIGGIKLRHQQWSLVMESVVPSDRGNYTCVVENKFGSIRQTYTLDVLERSPHRPILQAGLPANQTAVLGSDVEFHCKVYSDAQPHIQWLKHVEVNGSKVGPDGTPYVTVLKTAGANTTDKELEVLSLHNVTFEDAGEYTCLAGNSIGFSHHSAWLVVLPAEEELVEADEAGSVYAGILSYGVGFFLFILVVAAVTLCRLRSPPKKGLGSPTVHKISRFPLKRQQVSLESNASMSSNTPLVRIARLSSGEGPTLANVSELELPADPKWELSRARLTLGKPLGEGCFGQVVMAEAIGIDKDRAAKPVTVAVKMLKDDATDKDLSDLVSEMEMMKMIGKHKNIINLLGACTQGGPLYVLVEYAAKGNLREFLRARRPPGLDYSFDTCKPPEEQLTFKDLVSCAYQVARGMEYLASQKCIHRDLAARNVLVTEDNVMKIADFGLARDVHNLDYYKKTTNLVLWGPALGDLHAGGLPVPRHPCGGALQAAEGGPPHGQARQLHTRPVHDHAGVLACRALPEAHLQAAGGGPGPCPYRDVHRRVPGPVGAFRAVLPGWPGHPQLQLLRGRLRVCPRPAAPGPTQQWGLADVKGHWSPTM":
    print("SIUUUUUUUU")