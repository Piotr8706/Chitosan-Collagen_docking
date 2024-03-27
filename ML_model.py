import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json


"""
This script uses machine learning to predict biding energy between chitosan and
collagen type II.  Chitosan DD and collagen HD can be any number from 0 to 1.
Data is extrated using main.py script
 """  
#dataset = pd.read_json('./Files/Hydrogen_Bonds.json')
with open('./Files/Hydrogen_Bonds.json') as json_data:
    data = json.load(json_data)
#pd.DataFrame.from_dict(data, orient='index').T.set_index('index')
print(data['1']['125'])