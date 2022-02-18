import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.chdir(r'C:\Users\lga420\Desktop\python_test\GAPDH_REV')
data = pd.read_excel("20190814_GAPDH_R_BCA_py.xlsx")
#dtaclean: nan and coluns that have values for empty wells (H?)
del data['Temperature(°C)']
del data['Time(hh:mm:ss)']
data_clean=data.dropna(axis=1, how='all')
print(data_clean)
total_col=data_clean.shape[1]
samples=4
data_list=[]
data_list_mean=[]
data_list_background=[]
data_list_mean_background=[]
data_dilution=[]

dilutions=[0,5/240,10/240,15/240,20/240,25/240,30/240,40/240]
len_d=len(dilutions)
protein_BCA=[i * 2 for i in dilutions]

for x in range(0, samples):
    vec=list(range(x,total_col,samples))
    data_sample=data_clean.iloc[:,vec]
    data_sample=data_sample[data_sample.columns[data_sample.min() > 0.05]] #select data >0.05 eliminate empty wells
    data_sample_mean=data_sample.mean()
    data_sample_background = data_sample.subtract(data_sample_mean[0])
    data_sample_mean_background = data_sample_mean.subtract(data_sample_mean[0])
    data_list.append(data_sample)
    data_list_mean.append(data_sample_mean)
    data_list_background.append(data_sample_background)
    data_list_mean_background.append(data_sample_mean_background)

#organize data in a dataframe
bsa=pd.DataFrame(data_list_mean_background[0])
f1=pd.DataFrame(data_list_mean_background[1])
f2=pd.DataFrame(data_list_mean_background[2])
f4=pd.DataFrame(data_list_mean_background[3])
a=pd.concat([bsa, f1, f2, f4])
b=dilutions*samples
bsa_d=bsa.T/dilutions
a.reset_index(inplace=True)
a.insert(1, 'dilution', b, True)
a=pd.concat([a,pd.DataFrame(protein_BCA)], ignore_index=True, axis=1)
d=a.T
new_index=['coordenates','dilution','absorbance','[protein]']
d.index=(new_index)
a=d.T
f=4 #number of points considered
plt.title('BCA')
plt.xlabel('BCA mg/ml')
plt.ylabel('Abs 562')
x=(a.iloc[0:(f),3].values).astype('float') #abs
y=(a.iloc[0:(f),2].values).astype('float') #protein
plt.scatter(x,y) #discarded last point (bending)
#plt.show()
p=np.poly1d(np.polyfit(x,y,1)) #fit
print(a)
#8--> /8.273 - 0.09101
#6--> /8.87 - 0.03473
#5--> /8.975 - 0.02886
#4--> /9.729 - 0.002528
#3--> /9.916 - 0.007745
sample_1=(a.iloc[len_d:(len_d*2), 2])/9.729 - 0.002528
sample_2=(a.iloc[(len_d*2):(len_d*3), 2])/9.729 - 0.002528
sample_3=(a.iloc[(len_d*3):(len_d*4), 2])/9.729 - 0.002528
a.loc[len_d:(len_d*2),'[protein]'] = sample_1
a.loc[(len_d*2):(len_d*3),'[protein]'] = sample_2
a.loc[(len_d*3):(len_d*4),'[protein]'] = sample_3
a.replace(0,np.nan, inplace=True)
corrected_values=a.loc[:,'[protein]']/a.loc[:,'dilution']
a.insert(1, 'total [protein]', corrected_values, True)
a.fillna(0, inplace=True)
f2=1
f4=2

protein_f2=(a.iloc[(len_d+1):(len_d*2), 1])*f2
protein_f4=(a.iloc[(len_d*2+1):(len_d*3), 1])*f4
protein_mean_f2=protein_f2.mean()
protein_mean_f4=protein_f4.mean()
total_protein=(protein_mean_f2+protein_mean_f4)/2
print(total_protein) #mg/ml
print(protein_f2)
