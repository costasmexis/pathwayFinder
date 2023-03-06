import pandas as pd
import numpy as np
import networkx as nx
import io
from tqdm import tqdm

from Bio.KEGG import REST

def to_df(x):
    return pd.read_csv(io.StringIO(x), sep='\t', header=None)
    
compounds_list = REST.kegg_list('compound').read()
compounds_list = to_df(compounds_list)

compounds_list.rename(columns={0:'id', 1:'name'}, inplace=True)
compounds_list['mol_weight'] = np.nan
compounds_list['formula'] = np.nan

# save to csv
compounds_list.to_csv('data/compounds_list_KEGG.csv', index=None)

for row in tqdm(range(len(compounds_list))):
    
    request = REST.kegg_get('cpd:'+compounds_list.iloc[row]['id']).read()
    
    for line in request.split('\n'):
        if line.startswith('MOL_WEIGHT'):
            mol_weight = float(line.split()[1])
            compounds_list.iloc[row]['mol_weight'] = mol_weight
        if line.startswith('FORMULA'):
            formula = line.split()[1]
            compounds_list.iloc[row]['formula'] = formula
            
# save to csv
compounds_list.to_csv('data/compounds_list_KEGG.csv', index=None)
