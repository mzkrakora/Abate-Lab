#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import statements

import numpy as np
from datascience import *

import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
import warnings
warnings.simplefilter('ignore', FutureWarning)

import csv


# In[2]:


#creating tables

# with open("names.csv", newline='') as csvfile:
#     namereader = csv.reader(csvfile, delimiter = ' ', quotechar = '|')
#     name=[]
#     for row in namereader:
        
#         print(', '.join(row))


# In[3]:


umi_corrected_phage = Table.read_table("umi_corrected_phage.csv")
len(umi_corrected_phage)
umi_corrected_phage= umi_corrected_phage.drop("Unnamed: 44")
umi_corrected_phage


# In[ ]:


len(list(umi_corrected_phage[0]))
#Goal is to plot count by antibody and then have error bars there as well
umi_corrected_phage.bar()


# In[10]:


with open('umi_corrected_phage.csv', 'r') as f:
    data = [line.strip().split(',') for line in f]


# In[11]:


data2 = np.array(data)[:,:-1].astype(float)


# In[15]:


data2[data2>=1]=1
#duplicate counts normalized to >= 1


# In[17]:


plt.hist(data2.sum(axis=1), bins=44)
plt.yscale("log")


# In[ ]:


#the max phage in each cell cannot be more than 10% of the total phage population

 
names = ["ANPEP.01","ANPEP.02","ANTXR1.01","ANTXR1.04","ANTXR1.10","EPHA7.01","EPHA7.02","EPHA7.03","EPHA7.04","FGFR2.02","FGFR4.01","FGFR4.03","FGFR4.04","FLT3.01","FLT3.02","FLT3.03","GFP.01","GFP.02","GFP.03","GFP.04","GFP.05","ICAM1.10","ICAM1.13","ICAM1.20","INSR.01","NCR3LG1.01","NCR3LG1.06","NCR3LG1.12","NCR3LG1.14","PDGFRB.01","PDGFRB.02","PDGFRB.03","PDL1.02","PDL1.03","RET.01","RET.02","ROR1.02","ROR1.09","ROR1.11","ROR1.15","TRAlpha.01","ZNF18.01","ZNF2.01", "X44"]
GFP = names[16:21]
GFP

