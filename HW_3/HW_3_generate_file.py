from enum import unique
import pandas as pd
import re

# read csv file
data = pd.read_csv('TCGA_mutprop_cancer_race_info.csv')

# retrieve all distinct cancer names (set())
data_header = data.columns
cancer_types = set()
for col_header in data_header:
    try:
        cancer_name = re.search("TCGA-(.+?)\(", col_header).group(1)
        cancer_types.add(cancer_name)
    except AttributeError:
        pass

# initialize the count list
asian_gen_cancer_dict = {}
black_gen_cancer_dict = {}
asian_gen_cancer_name_dict = {}
black_gen_cancer_name_dict = {}
for index, row in data.iterrows():
    gene_name = row['gene_name']
    asian_gen_cancer_dict[gene_name] = 0
    black_gen_cancer_dict[gene_name] = 0
    asian_gen_cancer_name_dict[gene_name] = []
    black_gen_cancer_name_dict[gene_name] = []

for cancer_type in cancer_types:

    # find all (usually three) columns that are associated with a cancer
    # e.g. When cancer_type = "ACC", 
    # the result is ["TCGA-ACC(Asian)(2)","TCGA-ACC(Black Or African American)(1)","TCGA-ACC(White)(78)"]
    cancer_type_cols = [col_header for col_header in data.columns if cancer_type in col_header]
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
            asian_gen_cancer_name_dict[gene_name].append(cancer_type)
        if row[black_col_name] != 0 and row[white_col_name] == 0:
            black_gen_cancer_dict[gene_name] += 1
            black_gen_cancer_name_dict[gene_name].append(cancer_type)

    #print(f'{cancer_type}: Unique_asian_mutation: {unique_asian_mutation}; Unique_black_mutation: {unique_black_mutation}')
# print end result for every gene

column_names = ["gene_name",
                "Unique Asian Mutation",
                "Unique Asian Cancers",
                "Unique Black Mutation",
                "Unique Black Cancers"]

output_df = pd.DataFrame(columns = column_names)
count=0
for index, row in data.iterrows():
    count+=1
    gene_name = row['gene_name']
    asian_gene_mutation = asian_gen_cancer_dict[gene_name]
    cancers_for_asians = ", ".join(asian_gen_cancer_name_dict[gene_name])
    black_gene_mutation = black_gen_cancer_dict[gene_name]
    cancers_for_blacks = ", ".join(black_gen_cancer_name_dict[gene_name])
    print(f'{index} {count} {gene_name} *** {cancers_for_asians}')

    new_df = pd.DataFrame({'gene_name' : [gene_name], 
               'Unique Asian Mutation' : [asian_gene_mutation],
               'Unique Asian Cancers' : [cancers_for_asians],
               'Unique Black Mutation' : [black_gene_mutation],
               'Unique Black Cancers' : [cancers_for_blacks]})

    output_df = pd.concat([output_df, new_df], ignore_index = True, axis = 0)

# write to csv file
output_df.to_csv("./output.csv", index = False)