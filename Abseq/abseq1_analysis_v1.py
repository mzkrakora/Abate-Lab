#!/usr/bin/env python
# coding: utf-8

# # Analysis of AbSeq1 Experiment
# Ben Demaree 9.21.2018
# 
# *note: must have run abseq pipeline beforehand with data exported to sqlite database*

# In[1]:


# imports
from __future__ import division
import itertools
import numpy as np
import pandas as pd
import seaborn as sns; sns.set(style="white", color_codes=True)
from sqlalchemy import create_engine
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.family'] = 'Helvetica'
plt.rcParams['axes.unicode_minus'] = False
from sklearn.manifold import TSNE
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.decomposition import PCA


# ### 1. Protein data analysis
# Compare distributions of CD19, CD30, CD45 across all cells.

# In[ ]:


# begin by importing parsed data from sqlite database
db_path = '/home/bdemaree/missionbio/abseq1/dbs/abseq1.db'
# create DB interface
db = create_engine('sqlite:///' + db_path)
# import ab data
ab_df = pd.read_sql_query('''select * from umi_counts''', db)
# print list of column names (different methods of counting umis)
print list(ab_df) 


# In[100]:


# transform ab data into counts by cell
ab_counts = pd.pivot_table(ab_df, values='adjacency', index=['cell_barcode'], columns=['ab_barcode'],aggfunc=np.sum).fillna(0)
#ab_counts.head
ab_combos = [list(x) for x in itertools.combinations(list(ab_counts), 2)]


# plot scatters of ab combinations - **all cells**

# In[101]:


# plot scatters of all ab combinations - all cells
min_total_count = 20  # minimum sum of all ab counts for a given cell
ax_max = 400
# plot scatters of all ab combinations - all cells
for c in ab_combos:
    # restrict to valid cells as identified by the MB pipeline
    sns.jointplot(x=ab_counts[ab_counts.sum(axis=1) > min_total_count][c[0]], 
                 y=ab_counts[ab_counts.sum(axis=1) > min_total_count][c[1]],
                 space=0,
                 joint_kws={'alpha': 0.1},
                 marginal_kws=dict(bins=np.linspace(0, ax_max, 20)),
                 xlim=(0,ax_max),
                 ylim=(0,ax_max))


# plot scatters of ab combinations - **valid cells only**

# In[102]:


# load mission bio cell tsv into dataframe and use to plot valid cells only
cell_tsv_path = '/home/bdemaree/missionbio/abseq1/mb_data/abseq1_180921021949.barcode.cell.distribution.merged.tsv'
cell_tsv = pd.read_table(cell_tsv_path, sep='\t', header=0, index_col=0)
#cell_tsv.head
#cell_tsv.shape


# In[103]:


# plot scatters of all ab combinations - valid cells
ax_max = 100
for c in ab_combos:
    # restrict to valid cells as identified by the MB pipeline
    sns.jointplot(x=ab_counts[ab_counts.index.isin(cell_tsv.index)][c[0]], 
                 y=ab_counts[ab_counts.index.isin(cell_tsv.index)][c[1]],
                 space=0,
                 joint_kws={'alpha': 0.1},
                 marginal_kws=dict(bins=np.linspace(0, ax_max, 20)),
                 xlim=(0,ax_max),
                 ylim=(0,ax_max))


# In[2]:


good_ab_counts = ab_counts[ab_counts.index.isin(cell_tsv.index)].values.sum() # ab count from good cells
total_ab_counts = ab_counts.values.sum()  # ab count from all cells
frac_good = good_ab_counts/total_ab_counts  # fraction of ab counts from good cells
print 'Fraction of ab counts belonging to valid cells is %0.4f.' % frac_good


# In[105]:


# plot ab ratios (change 0 to 1 - total count must be greater than 10)
# normalize each antibody by its average count
min_total_count = 10
ab_counts_norm = ab_counts/ab_counts.mean(axis=0)
ab_counts_norm[ab_counts_norm == 0] += 1e-6
valid_rows = ab_counts_norm[ab_counts.index.isin(cell_tsv.index) & np.asarray(ab_counts.sum(axis=1).index > min_total_count)]
for c in ab_combos:
    x=np.divide(valid_rows[c[0]],
              valid_rows[c[1]])
    plt.figure()
    plt.hist(x, bins=np.linspace(0, 10, 40))
    plt.xlabel('count %s/count %s' % (c[0], c[1]))


# ### 2. SNP Analysis
# Use genotypying data to identify cell populations.

# In[106]:


# note: run loom_to_tsv.py beforehand
# import snp tsv into pandas
snp_tsv_path = '/home/bdemaree/missionbio/abseq1/mb_data/abseq1_180921021949.cells.tsv'
snp_tsv = pd.read_table(snp_tsv_path, sep='\t', header=0, index_col=0)
snp_tsv.columns = np.sort(cell_tsv.index.values)  # loom file barcodes are alphabetized


# In[108]:


# run a pca on the snp data
pca = PCA(n_components=10)  # 10 dimensions
principal_components = pca.fit_transform(np.transpose(snp_tsv))
#principal_components.shape


# In[138]:


# plot the pca results
plt.scatter(principal_components[:,0], principal_components[:,1], alpha=0.1)


# In[110]:


# run tsne on the dataset
tsne = TSNE(metric='correlation', perplexity=50, random_state=1, verbose=1)
tsne_results = tsne.fit_transform(np.transpose(snp_tsv))             # perform tsne computation


# In[139]:


# plot the tsne results
plt.scatter(tsne_results[:, 0], tsne_results[:, 1])  # create scatter plot


# In[140]:


# hierarchical clustering of tsne

# generate the linkage matrix
Z = linkage(tsne_results, 'ward')
max_d = 500      # cluster cutoff distance

clusters = fcluster(Z, max_d, criterion='distance')     # extract flat clusters
plt.figure(figsize=(10, 8))
x = tsne_results[:, 0]
y = tsne_results[:, 1]
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.scatter(x, y, c=clusters, cmap='Set2')  # plot points with cluster dependent colors

# add cluster labels to graph
for cl in set(clusters):
    pts = [ind for ind in range(len(clusters)) if clusters[ind] == cl]      # find cells in this cluster
    x_vals = [x[i] for i in pts]
    y_vals = [y[i] for i in pts]
    center = (np.mean(x_vals), np.mean(y_vals))     # calculate center
    t = plt.annotate(str(cl), xy=center, fontsize=10)
    t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor=None))


# In[128]:


# make a dictionary of clustering data
# key: cluster_num, value: list of cell barcodes
cluster_dict = {}
for c in list(set(clusters)):
    cluster_dict[c] = [snp_tsv.columns.values[i] for i in range(len(clusters)) if clusters[i] == c]
    
# cluster 1 is Raji, cluster 2 is K562


# In[129]:


# does the tsne contain true groups of different cell types?

ax_max = 100

# replot abs using clustering from tsne or pca plots

# simple slicing
# pc_df = pd.DataFrame(data=principal_components, index=snp_tsv.columns, columns=['x', 'y'])
# tsne_df = pd.DataFrame(data=tsne_results, index=snp_tsv.columns, columns=['x', 'y'])
#pca_cells = pc_df[pc_df.x > 20].index
#tsne_cells = tsne_df[tsne_df.y < -30].index
# for c in ab_combos:
#     # restrict to valid cells as identified by the MB pipeline
#     sns.jointplot(x=ab_counts.loc[tsne_cells][c[0]],
#                  y=ab_counts.loc[tsne_cells][c[1]],
#                  space=0,
#                  joint_kws={'alpha': 0.1},
#                  marginal_kws=dict(bins=np.linspace(0, ax_max, 20)),
#                  xlim=(0,ax_max),
#                  ylim=(0,ax_max))

# all clusters
plt.figure(figsize=(20, 50))
sub_count = 1
for c in ab_combos:
    plt.subplot(int('13%s' % sub_count))
    plt.scatter(ab_counts.loc[snp_tsv.columns.values][c[0]],
                ab_counts.loc[snp_tsv.columns.values][c[1]],
                cmap='Set2',
                c=clusters,
                alpha=0.1)
    plt.xlim([0, ax_max])
    plt.ylim([0, ax_max])
    plt.xlabel(c[0])
    plt.ylabel(c[1])
    plt.title('All tSNE Clusters')
    plt.gca().set_aspect('equal')
    sub_count += 1

# clusters 1 and 2
colors = list(iter(plt.cm.Set2(np.linspace(0,1,len(set(clusters))))))
plt.figure(figsize=(20, 50))
sub_count = 1
for c in ab_combos:
    plt.subplot(int('13%s' % sub_count))
    plt.scatter(ab_counts.loc[cluster_dict[1] + cluster_dict[2]][c[0]],
                ab_counts.loc[cluster_dict[1] + cluster_dict[2]][c[1]],
                c=[colors[0]]*len(cluster_dict[1]) + [colors[1]]*len(cluster_dict[2]),
                alpha=0.1)
    plt.xlim([0, ax_max])
    plt.ylim([0, ax_max])
    plt.xlabel(c[0])
    plt.ylabel(c[1])
    plt.title('tSNE Clusters 1 (Raji) and 2 (K562)')
    plt.gca().set_aspect('equal')
    sub_count += 1
    
# combine some clusters we think are K562 and Raji
# plt.figure(figsize=(20, 50))
# sub_count = 1
# k562_clusters = cluster_dict[2] + cluster_dict[3]
# raji_clusters = cluster_dict[1] + cluster_dict[4]
# for c in ab_combos:
#     plt.subplot(int('13%s' % sub_count))
#     plt.scatter(ab_counts.loc[k562_clusters + raji_clusters][c[0]],
#                 ab_counts.loc[k562_clusters + raji_clusters][c[1]],
#                 c=[colors[0]]*len(k562_clusters) + [colors[1]]*len(raji_clusters),
#                 alpha=0.1)
#     plt.xlim([0, ax_max])
#     plt.ylim([0, ax_max])
#     plt.xlabel(c[0])
#     plt.ylabel(c[1])
#     plt.gca().set_aspect('equal')
#     sub_count += 1


# In[122]:


# perform hierarchical clustering of all genotyping data

# generate the linkage matrix
Z = linkage(np.transpose(snp_tsv), 'ward')


# In[133]:


# from https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('cell index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')

# plot dendrogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('cell')
plt.ylabel('distance')

max_d = 600      # cluster cutoff distance

# generate dendrogram
fancy_dendrogram(
    Z,
    leaf_rotation=90.,          # rotates the x axis labels
    leaf_font_size=8.,          # font size for the x axis labels
    truncate_mode='level',      # only show last p levels
    p=8,                        # levels to show
    max_d=max_d                 # plot cutoff line
)


# In[135]:


# replot the ab scatters
clusters = fcluster(Z, max_d, criterion='distance')     # extract flat clusters
print 'number of clusters: %d' % len(set(clusters))
# two clusters
plt.figure(figsize=(20, 50))
sub_count = 1
for c in ab_combos:
    plt.subplot(int('13%s' % sub_count))
    plt.scatter(ab_counts.loc[snp_tsv.columns.values][c[0]],
                ab_counts.loc[snp_tsv.columns.values][c[1]],
                cmap='Set2',
                c=clusters,
                alpha=0.1)
    plt.xlim([0, ax_max])
    plt.ylim([0, ax_max])
    plt.xlabel(c[0])
    plt.ylabel(c[1])
    plt.title('hierarchical clustering')
    plt.gca().set_aspect('equal')
    sub_count += 1

