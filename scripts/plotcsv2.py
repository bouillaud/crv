import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 

df = pd.read_csv("hits.csv")

# plt.scatter(df['y'][df['edep']>1e-3], df['z'][df['edep']>1e-3])

pid_type = []
pid_counts = {}
for p in df['pid']:
	if p not in pid_type:
		pid_type.append(p)
		pid_counts[p] = 0
	else:
		pid_counts[p] += 1

print(pid_type)
print(pid_counts)

plt.scatter(df['y'], df['z'], s=8+np.log10(df['edep']), c=df['x'])
plt.show() 

plt.scatter(df['y'][df['pid']==pid_type[0]], df['z'][df['pid']==pid_type[0]], s=8+np.log10(df['edep'][df['pid']==pid_type[0]]), c=df['x'][df['pid']==pid_type[0]])
plt.show()