import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.decomposition import PCA

import pandas as pd

def plot_3d(df):
	df = df.drop(['doc_name'], axis=1)
	X = df.drop(['type'],axis=1)
	y = df['type']
	
	pca = PCA(n_components=3)
	X_r = pca.fit(X).transform(X)
	lw=2

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	colors = ['navy', 'darkorange']
	for color, i, target_name in zip(colors, ['burst', 'irregular', 'tonic'], ['burst', 'irregular', 'tonic']):
		ax.scatter(X_r[np.array(y == i), 0], X_r[np.array(y == i), 1], X_r[np.array(y == i), 2], color=color, alpha=.8, lw=lw,
					label=target_name)
	
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	plt.show()

def plot_2d(df):
	X = df[['cv', 'burst_mean']].values
	y = df['type']
	
	lw=2
	
	colors = ['navy', 'darkorange', 'green']
	for color, i, target_name in zip(colors, ['burst', 'irregular', 'tonic'], ['burst', 'irregular', 'tonic']):
		plt.scatter(X[np.array(y == i), 0], X[np.array(y == i), 1], color=color, alpha=.8, lw=lw,
					label=target_name)
	
	ax = plt.axes()
	
	ax.set_xlabel('cv')
	ax.set_ylabel('burst')
	
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	plt.show()

df=pd.read_csv('all_res.csv', sep=',')

plot_3d(df)	