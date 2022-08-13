import pandas as pd
import re

# read csv file
data = pd.read_csv('TCGA_mutprop_cancer_race_info.csv')

# retrieve all distinct cancer names (set())
data_header = data.columns
cancer_types = set()
for col in data_header:
    try:
        cancer_name = re.search("TCGA-(.+?)\(", col).group(1)
        cancer_types.add(cancer_name)
    except AttributeError:
        pass

# initialize the count list
asian_gen_cancer_dict = {}
black_gen_cancer_dict = {}
for index, row in data.iterrows():
    gene_name = row['gene_name']
    asian_gen_cancer_dict[gene_name] = 0
    black_gen_cancer_dict[gene_name] = 0

for cancer_type in cancer_types:
    # find all (usually three) columns that are associated with a cancer
    # e.g. When cancer_type = "ACC", the result is ["TCGA-ACC(Asian)(2)","TCGA-ACC(Black Or African American)(1)","TCGA-ACC(White)(78)"]
    cancer_type_cols = [col for col in data.columns if cancer_type in col]
    # find the column of the particular race
    for cancer_type_col in cancer_type_cols:
        if "Asian" in cancer_type_col:
            asian_col_name = cancer_type_col
        elif "Black" in cancer_type_col:
            black_col_name = cancer_type_col
        elif "White" in cancer_type_col:
            white_col_name = cancer_type_col

    # identify if mutation is unique to race
    for index, row in data.iterrows():
        gene_name = row['gene_name']
        # increase number of unique mutation occurences for Asian and Black race
        if row[asian_col_name] != 0 and row[white_col_name] == 0:
            asian_gen_cancer_dict[gene_name] += 1
        if row[black_col_name] != 0 and row[white_col_name] == 0:
            black_gen_cancer_dict[gene_name] += 1

# print end result for every gene
for index, row in data.iterrows():
    gene_name = row['gene_name']
    print(f'{gene_name}: Unique Asian mutation: {asian_gen_cancer_dict[gene_name]}; Unique Black mutation: {black_gen_cancer_dict[gene_name]}')