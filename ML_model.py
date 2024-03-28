import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split


"""
This script uses machine learning to predict biding energy between chitosan and
collagen type II.  Chitosan DD and collagen HD can be any number from 0 to 1.
Data is extrated using main.py script
 """  
#dataset = pd.read_json('./Files/Hydrogen_Bonds.json')
df = pd.read_csv('./Files/All_data_combined.csv')
print(df.describe())
print(df.info())
print(df.isnull().values.any())
sns.countplot(df['Binding_Energy'])
plt.plot
plt.show()
x = df.drop(columns='Binding_Energy')
y = df['Binding_Energy']
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=0)
print(x_train.shape)
