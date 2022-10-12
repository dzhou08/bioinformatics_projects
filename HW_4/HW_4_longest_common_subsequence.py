from enum import unique
from operator import itemgetter
import pandas as pd
import re
import numpy as np
import itertools as it
 
input_cancer_names = input("What are the cancers? ")
race = input("What's the race (Options: Asians or Blacks)?: ")

data = pd.read_csv('cancer_gene.csv')
cancer_names = input_cancer_names.split(', ')

cancer_dict = {}

for cancer in cancer_names:
    value_list = data.loc[data['cancer_name'] == cancer, 'uniquely_mutated_genes_for_' + race.lower()].iloc[0].replace("[", "").replace("]", "").replace(" ", "").replace("'", "").split(',')
    #print(value_list) 
    cancer_dict[cancer] = set(value_list)

venn_diagram_column_names = ['cancer_combinations',
                            'count_of_common_genes',
                            'common_genes']

venn_df = pd.DataFrame(columns = venn_diagram_column_names)

for i in range(len(cancer_names), 0, -1):
    for item in list(it.combinations(cancer_names, i)):
        count = 0
        for j in list(item):
            if count > 0:
                return_value = return_value.intersection(cancer_dict[j])
            else:
                return_value = cancer_dict[j]
            count += 1
        
        new_df = pd.DataFrame({'cancer_combinations' : [str(list(item))],
                   'count_of_common_genes' : [len(return_value)], 
                   'common_genes' : [return_value]})

        venn_df = pd.concat([venn_df,new_df], ignore_index = True, axis = 0)
venn_df.to_csv('venn_' + race.lower() + '.csv', index = False)