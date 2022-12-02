from enum import unique
import pandas as pd
import re
import numpy as np

def read_race_gene_file(cancer,race):
    try:
        file_name = 'input_data/' + cancer + '_' + race + '_Vs_White.csv'
        data = pd.read_csv(file_name, index_col = 0)
        print(file_name)
        # only keep the p_value < 0.05 rows
        data = data[data['p value'] < 0.05].rename(columns = {"p value": "p value " + race})
        # pd.merge(asian_data, black_data, how ='inner', on =['Gene', 'Gene'])

        mutationData = pd.read_csv('TCGA_mutprop_cancer_race_info.csv', index_col = 0)
        colName = 'TCGA-' + cancer + "(" + race
        
        data_header = mutationData.columns
        for col in data_header:
            try:
                
                if col.startswith(colName):
                    mutationDataSubset = mutationData[['gene_name', col]]

                    #data.merge(mutationDataSubset, left_on='Gene', right_on='gene_name', )
                    inner_merged = pd.merge(data, mutationDataSubset, left_on = ["Gene"], right_on = ["gene_name"])
                    inner_merged = inner_merged[['Gene', 'p value ' + race, col]]
                    inner_merged = inner_merged.sort_values([col,"Gene"], ascending = [False,True])
                    inner_merged = inner_merged.rename(columns = {col: "mutation rates"})
                    inner_merged.to_csv('output_data/' + cancer + '_' + race + '_genes.csv', index = False) 
                    return

            except KeyError:
                pass

    except FileNotFoundError:
        pass

    # return data
    
cancers = input("What are the cancers? (All uppercase) ")

if cancers == "ALL":
    data = pd.read_csv('TCGA_mutprop_cancer_race_info.csv')
    data_header = data.columns
    cancers = set()
    for col_header in data_header:
        try:
            cancer_name = re.search("TCGA-(.+?)\(", col_header).group(1)
            cancers.add(cancer_name)
        except AttributeError:
            pass
else:
    cancers = cancers.split()

for cancer in cancers:
    # read Asian gene file
    read_race_gene_file(cancer, "Asian")

    # read Black gene file
    read_race_gene_file(cancer, "Black")

    
        