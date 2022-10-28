from enum import unique
import pandas as pd
import re
import numpy as np

# read csv file
data = pd.read_csv('TCGA_mutprop_cancer_race_info.csv')

# retrieve all distinct cancer names (set())

for column in data:
    try:
        cancer_name = re.search("TCGA-(.+?)\(", column).group(1)

        
        new_df = data[['gene_name', column]]
        new_df = new_df.sort_values(by=[column], ascending = False)

        if 'Asian' in column:
            new_df.to_csv("./output/" + cancer_name + "_Asian" + ".csv", index = False)

        elif 'Black' in column:
            new_df.to_csv("./output/" + cancer_name + "_Black" + ".csv", index = False)

    except AttributeError:
        pass


















'''
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
    asian_col_name, black_col_name, white_col_name = '','',''
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
        if asian_col_name != '' and white_col_name != '' \
            and row[asian_col_name] != 0 and row[white_col_name] == 0:
            asian_gen_cancer_dict[gene_name] += 1
            asian_gen_cancer_name_dict[gene_name].append(cancer_type)
        if black_col_name != '' and white_col_name != '' \
            and row[black_col_name] != 0 and row[white_col_name] == 0:
            black_gen_cancer_dict[gene_name] += 1
            black_gen_cancer_name_dict[gene_name].append(cancer_type)

    #print(f'{cancer_type}: Unique_asian_mutation: {unique_asian_mutation}; Unique_black_mutation: {unique_black_mutation}')
# print end result for every gene

column_names = ["gene_name",
                "unique_asian_mutation",
                "unique_asian_cancers",
                "unique_black_mutation",
                "unique_black_cancers"]

output_df = pd.DataFrame(columns = column_names)
count = 0
for index, row in data.iterrows():
    count += 1
    gene_name = row['gene_name']
    asian_gene_mutation = asian_gen_cancer_dict[gene_name]
    cancers_for_asians = ", ".join(asian_gen_cancer_name_dict[gene_name])
    black_gene_mutation = black_gen_cancer_dict[gene_name]
    cancers_for_blacks = ", ".join(black_gen_cancer_name_dict[gene_name])
    print(f'{index} {count} {gene_name} *** {cancers_for_asians}')

    new_df = pd.DataFrame({'gene_name' : [gene_name], 
               'unique_asian_mutation' : [asian_gene_mutation],
               'unique_asian_cancers' : [cancers_for_asians],
               'unique_black_mutation' : [black_gene_mutation],
               'unique_black_cancers' : [cancers_for_blacks]})

    output_df = pd.concat([output_df, new_df], ignore_index = True, axis = 0)

# write to csv file
output_df.to_csv("./output.csv", index = False)


venn_diagram_column_names = ['cancer_name',
                            'count_of_unique_asian_genes',
                            'uniquely_mutated_genes_for_asians',
                            'count_of_unique_black_genes',
                            'uniquely_mutated_genes_for_blacks']

v_df = pd.DataFrame(columns = venn_diagram_column_names)
v_df = v_df.assign(cancer_name = list(cancer_types),
                    count_of_unique_asian_genes = [0 for x in range(len(cancer_types))],
                    uniquely_mutated_genes_for_asians = [[] for x in range(len(cancer_types))],
                    count_of_unique_black_genes = [0 for x in range(len(cancer_types))],
                    uniquely_mutated_genes_for_blacks = [[] for x in range(len(cancer_types))])
venn_diagram_output_data = pd.read_csv('output.csv')
for index, row in venn_diagram_output_data.iterrows():
    gene_name = row['gene_name']
    asian_cancer_mutations_str = row['unique_asian_cancers']
    if asian_cancer_mutations_str == asian_cancer_mutations_str:
        for cancer in asian_cancer_mutations_str.split(', '):
            m = v_df['cancer_name'].eq(cancer)
            
            v_df['uniquely_mutated_genes_for_asians'] = v_df['uniquely_mutated_genes_for_asians']\
                                                        .mask(m, 
                                                            v_df['uniquely_mutated_genes_for_asians']
                                                          .apply(lambda x: list(x) + [str(gene_name).strip()]))

            v_df['count_of_unique_asian_genes'] = v_df['count_of_unique_asian_genes']\
                                                        .mask(m, 
                                                            v_df['count_of_unique_asian_genes']
                                                            .apply(lambda x: x + 1))
    
    black_cancer_mutations_str = row['unique_black_cancers']
    if black_cancer_mutations_str == black_cancer_mutations_str:
        for cancer in black_cancer_mutations_str.split(', '):
            m = v_df['cancer_name'].eq(cancer)

            v_df['uniquely_mutated_genes_for_blacks'] = v_df['uniquely_mutated_genes_for_blacks']\
                                                        .mask(m, 
                                                            v_df['uniquely_mutated_genes_for_blacks']
                                                            .apply(lambda x: list(x) + [str(gene_name).strip()]))

            v_df['count_of_unique_black_genes'] = v_df['count_of_unique_black_genes']\
                                                        .mask(m, 
                                                            v_df['count_of_unique_black_genes']
                                                            .apply(lambda x: x + 1))

print(v_df)
v_df.to_csv("cancer_gene.csv")
'''