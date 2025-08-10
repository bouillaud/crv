import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

df = pd.read_csv("hits.csv")

# plt.scatter(df['y'][df['edep']>1e-3], df['z'][df['edep']>1e-3])
# plt.scatter(df['y'][df['edep']>1e-3], df['z'][df['edep']>1e-3])
# plt.scatter(df['y'][df['edep']>1e-3], df['z'][df['edep']>1e-3])

plt.scatter(df['y'], df['z'], s=8+np.log10(df['edep']), c=df['x'])

plt.show() 