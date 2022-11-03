from enum import unique
import pandas as pd
import re
import numpy as np

def read_race_gene_file(race):
    data = pd.read_csv(cancer + '_' + race + '_Vs_White.csv', index_col = 0)
    # only keep the p_value < 0.05 rows
    data = data[data['p value'] < 0.05].rename(columns={"p value": "p value " + race})
    return data
    
cancers = input("What are the cancers? (All uppercase) ").split()

for cancer in cancers:
    # read Asian gene file
    asian_data = read_race_gene_file("Asian")

    # read Black gene file
    black_data = read_race_gene_file("Black")
    
    # only keep the common genes for both races
    int_df = pd.merge(asian_data, black_data, how ='inner', on =['Gene', 'Gene'])
    int_df.to_csv(cancer + '_asian_and_black_common_genes.csv', index = False)


    
        