from sklearn import datasets
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt
from matplotlib import figure
import pandas as pd
import numpy as np

#Input: normalized and with no nan values dataframe, title for plot
def tsne_measurements (num_data, tl):
	y = num_data.loc[:, 'asthma']
	#Fit and transform with a TSNE
	tsne = TSNE(perplexity=15, init='pca')
	#fit tsne in the numeric data
	X = num_data.loc[:,num_data.columns != 'asthma'] 
	X = X.loc[:,X.columns != 'currentsmoking']
	X = X.loc[:,X.columns != 'ics.use']
	
	#run tsne with the features of the fabia's biclusters
	#lung
	#X = X.loc[:,['fvc.post','fvc','tlcoc.sb','fev1.post','fev1','amp.fev1.base','fef75','fef75.post',
	#'fev1.perc.pred.post','fvc.perc.pred.post','fvc.perc.pred','fef25','fef25.post','fef50']]
	#'fvc.post','fvc','tlcoc.sb','fev1.post','fev1','amp.fev1.base'
	#'fev1.perc.pred.post','fvc.perc.pred.post','fvc.perc.pred'
	#'fef75','fef75.post'
	#'fef25','fef25.post','fef50'

	#sputum
	#X = X.loc[:,['sp.viab.v2','sp.squamous.v2','sp.viab.v1','sp.squamous.v1']]

	#blood
	#X = X.loc[:,['ld.lym.perc','ld.seg.perc','ld.mon.perc','ld.eoc.perc']]
	
	X_2d = tsne.fit_transform(X)

	#2 classes: healthy, asthma	
	y_data = pd.Series()	
	counter = 0
	for i in y:
		if i == 'healthy':
			y_data.loc[counter]= 'healthy'
		else:
			y_data.loc[counter]='asthma'
		counter = counter+1

	target_ids = ['healthy','asthma']
	#Visualize the data
	colors = 'g', 'r'
	create_plot(target_ids,colors,tl,X_2d,y_data)

	#4 classes: current asthma, clinical remission, complete remission, healthy
	y_data = y
	target_ids = ['current asthma','clinical remission','complete remission','healthy']
	#Visualize the data
	colors = 'r', 'b','y','g'
	create_plot(target_ids,colors,tl,X_2d,y_data)

	#2 classes, smoker, nonsmoker
	y_data = num_data.loc[:, 'currentsmoking']
	target_ids = ['smoker','nonsmoker']
	#Visualize the data
	colors = 'r', 'g'
	create_plot(target_ids,colors,tl,X_2d,y_data)

	#2 classes, ics use, no ics use
	y_data = num_data.loc[:, 'ics.use']	

	'''
	#remove the datapoints that correspond to healthy people (we are only interested in patients that receive medicine or not)
	y.index = range(y.count())
	y_data.index = range(y_data.count()) #we set row indices again, because they had missing numbers

	index = 0
	y_list = []
	indices_healthy = []
	for i in y:
		if i!='healthy':
			y_list.append(y_data[index])
		else:
			indices_healthy.append(index)
		index = index + 1
	y_data_new = pd.DataFrame({'ics.use':y_list}, columns=['ics.use'])
	'''

	target_ids = ['ics use','no ics use']
	#Visualize the data
	colors = 'r', 'g'
	create_plot(target_ids,colors,tl,X_2d,y_data)


#Input: target ids, colors for data points, title, the data to be plotted
def create_plot(target_ids, colors,tl,X_2d,y_data):
	plt.figure(figsize=(10,5))
	for i, c, label in zip(target_ids, colors, target_ids):
		plt.scatter(X_2d[y_data == i, 0], X_2d[y_data == i, 1], c=c, label=label)
	plt.legend()
	plt.title(tl)
	plt.show()

#Input: initial dataframe, specific dataframe with  measurements
#Output: normalized specific dataframe with no nan values and with label column
def remove_nan(df,X):
	row_count = len(df.index)
	data = pd.DataFrame()

	for column in X:
		nan_count = X[column].isna().sum() #count of nan values of a column
		if nan_count < row_count*0.35: #if the nans are less than 35% of the column values, then keep the column
			data = data.append(X[column])	#with append, matrix  becomes transpose

	data = data.append(df.loc[:, 'asthma']) #add label columns (asthma, smoking, medicine)
	data = data.append(df.loc[:, 'currentsmoking'])

	#fill 'ics.use' nan values for healthy people (set to no ics use)
	med = df.loc[:, 'ics.use']
	med = med.fillna('no ics use')

	data = data.append(med)

	num_data = pd.DataFrame()
	index = 1
	list_index = []
	for column in data:
		nan_count = data[column].isna().sum() #count of nan values of a column
		if nan_count == 0:	#keep only the rows with no nan values
			num_data = num_data.append(data[column]) #with append, matrix becomes transpose again, so we have the initial format again
			list_index.append(index)

		index = index+1

	#normalize numeric data
	num_data_tmp = num_data.loc[:,num_data.columns != 'asthma']
	num_data_tmp = num_data_tmp.loc[:,num_data_tmp.columns != 'currentsmoking']
	df_proc = normalization(num_data_tmp.loc[:,num_data_tmp.columns != 'ics.use'])	
	#add asthma column label in the normalized data

	lbl1 = (num_data.loc[:,'asthma']).tolist()
	lbl2 = (num_data.loc[:,'currentsmoking']).tolist()
	lbl3 = (num_data.loc[:,'ics.use']).tolist()
	lbl_df = pd.DataFrame({'asthma':lbl1, 'currentsmoking':lbl2, 'ics.use':lbl3}, columns=['asthma','currentsmoking','ics.use']) #for correct reindexing
	df_proc = df_proc.join(lbl_df)
	
	'''
	with open ('row_biopsy.txt','w') as f:	#write in a txt file the indices for rows and cols to be kept (for not empty values)
		for item in list_index:
			f.write("%s\n" % item)
	with open ('col_biopsy.txt','w') as f:
		for item in list(df_proc):
			f.write("%s\n" % item)
	
	count_asthma = 0
	count_clinical = 0
	count_complete = 0
	count_healthy = 0
	for i in lbl1:
		if i == 'current asthma':
			count_asthma += 1
		elif i == 'clinical remission':
			count_clinical += 1
		elif i == 'complete remission':
			count_complete += 1
		else:
			count_healthy += 1
	print count_asthma
	print count_clinical
	print count_complete
	print count_healthy
	'''
	return df_proc	

#z-score normalization
#Input: dataframe to be normalized
#Output: normalized dataframe
def normalization(df):
	names = df.columns
	scaler = StandardScaler()
	scaled_df = scaler.fit_transform(df)
	scaled_df = pd.DataFrame(scaled_df, columns=names)
	return scaled_df

#preprocessing for the files GE.txt, METH.txt, microrna4_default_aggr_normalized_log.txt
#Input: the filename, a separator to read the file, the number of classes, the plot title
def preprocessing(filename,separator,tl):
	#processing GE.txt
	df_ge = pd.read_csv(filename, sep=separator)
	df_ge = df_ge.transpose() #transpose in order to have each patient in each row
	df_ge = df_ge.iloc[1:]	#drop the id column for now
	ids = df_ge.index

	rna_list = []
	for i in ids:	#assign to every patient his asthma label
		index = df.loc[df['rnaseq.id']==i]	#(Individual subjects can be linked to the rnaseq.id in the Database_biopten_v2_6.csv.) 
		value = index['asthma'].tolist()
		rna_list.append(value[0])

	rna_df = pd.DataFrame({'asthma':rna_list})	
	df_ge.index = range(len(df_ge))	#rename row char labels to integers
	df_ge = df_ge.join(rna_df)	#add asthma column to the gene dataframe

	rna_list = []
	for i in ids:	#assign to every patient his smoking label
		index = df.loc[df['rnaseq.id']==i]	#(Individual subjects can be linked to the rnaseq.id in the Database_biopten_v2_6.csv.) 
		value = index['currentsmoking'].tolist()
		rna_list.append(value[0])

	rna_df = pd.DataFrame({'currentsmoking':rna_list})	
	df_ge.index = range(len(df_ge))	#rename row char labels to integers
	df_ge = df_ge.join(rna_df)	#add currentsmoking column to the gene dataframe

	rna_list = []
	for i in ids:	#assign to every patient his ics.use label
		index = df.loc[df['rnaseq.id']==i]	#(Individual subjects can be linked to the rnaseq.id in the Database_biopten_v2_6.csv.) 
		value = index['ics.use'].tolist()
		rna_list.append(value[0])

	rna_df = pd.DataFrame({'ics.use':rna_list})	
	df_ge.index = range(len(df_ge))	#rename row char labels to integers
	df_ge = df_ge.join(rna_df)	#add currentsmoking column to the gene dataframe

	#call function for tsne
	tsne_measurements(df_ge,tl)


if __name__ == "__main__":
	#read csv to dataframe
	df = pd.read_csv('Database_biopten_v2_6.csv', sep=';')
	
	#blood measurements data
	tsne_measurements(remove_nan(df,df.loc[:, 'lc.ige':'triglycerids']),'Blood')

	#select sputum measurements data
	tsne_measurements(remove_nan(df,df.loc[:, 'sp.vol.v1':'sp.baso.abs.ml.v2']),'Sputum')

	#select lung function measurements data
	X_1 = df.loc[:,'ivc.pred':'pc20.meth']  #select dataframe slices that correspond to lung function data
	X_2 = df.loc[:,'amp.vc.base':'pc20.hist']
	X_lung = X_1.join(X_2)
	tsne_measurements(remove_nan(df,X_lung),'Lung Function')

	#select biopsy measurements data
	X_1 = (df.loc[:,'bm.thickness']).to_frame()	#select dataframe slices that correspond to biopsy data
	X_2 = df.loc[:,'mucus.area':'height.epithelium']
	X_3 = df.loc[:,'intact.perc':'other.perc']
	X_4 = df.loc[:,'np57':'cd20.infiltrate.whole.coupe']
	X_5 = df.loc[:,'aa1.cells.muscle':'foxp3']
	X_bio = (((X_1.join(X_2)).join(X_3)).join(X_4)).join(X_5)
	tsne_measurements(remove_nan(df,X_bio),'Biopsy')

	#perform preprocessing and tsne for GE.txt, METH.txt and microrna4_default_aggr_normalized_log.txt files
	preprocessing('GE.txt','\t','Gene expression')
	preprocessing('METH.txt','\t','DNA methylation')
	preprocessing('microrna4_default_aggr_normalized_log.txt',' ','miRNA expression')

	#lung function combined with sputum
	X_1 = df.loc[:, 'sp.vol.v1':'sp.baso.abs.ml.v2']
	X_new = X_1.join(X_lung)
	tsne_measurements(remove_nan(df,X_new),'Lung Function & Sputum')

	#lung function with blood
	X_2 = df.loc[:, 'lc.ige':'triglycerids']
	X_new = X_2.join(X_lung)
	tsne_measurements(remove_nan(df,X_new),'Lung Function & Blood')

	#blood with sputum
	X_new = X_1.join(X_2)
	tsne_measurements(remove_nan(df,X_new),'Sputum & Blood')

	#lung, blood, sputum
	X_new = X_new.join(X_lung)
	tsne_measurements(remove_nan(df,X_new),'Lung Function & Blood & Sputum')

