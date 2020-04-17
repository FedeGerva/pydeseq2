# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.2
#   kernelspec:
#     display_name: pyDESeq2_env
#     language: python
#     name: pydeseq2_env
# ---

# ## DESeq2 in python - pydeseq2

# Nowadays, it is becoming more clear the fundamental role of flexibility between programming languages in bioinformatics analysis. In fact, there are new packages that allow to call functions written in one language from another language. One of the most know package is [rpy2](https://rpy.sourceforge.io/rpy2/doc-dev/html/overview.html) which allow to use R function and call R object directly in a Python environment. These packages are of major importance because they allow the users to reuse existing code and perform analysis integrating different scripting languages.
#
# The bioinformatics community is growing and evolving day by day searching for new solutions that guarantee reproducibility, scalability and flexibility. A huge amount of bioinformatics tools developed in past are written in R language. However, in the last few years, python is emerged as a new promising languages for bioinformatics analysis, since it reduces the time and resources of the analysis. One of the most used bioinformatics packages for the NGS analysis (including RNA-seq, ChIP-seq, ATAC-seq, Hi-C etc) is [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). This package is part of the [Bioconductor project](https://www.bioconductor.org/) and it’s written in R language.
#
# Here I created an example of the DESeq2 analysis integrating R with python through rpy2 package.
#
# The toy dataset, the methods and plots are a reimplementation of [“Analyzing RNA-seq data with DESeq2”](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat)
#
# The script was run on [JupyterHub]( https://jupyter.org/hub) using a [conda]( https://docs.conda.io/en/latest/) environment(see below to have further information). Conda v4.5.11.
#
# The aims of this work are:
# * Generate a reproducible conda environment to run DESeq2 in Python 
# * Create a tutorial for DESeq2 analysis in Python using rpy2 package. This tutorial has been tested on RNA-seq and ATAC-seq data.
# * Replicate the analysis reported in DESeq2 [vignette]( https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
# * Reimplement DESeq2 canonical plots in python using [Matplotlib](https://matplotlib.org/) and [Seaborn](https://seaborn.pydata.org/)
# * Reduce the running time to produce DESeq2 object (such as dds and rld)
# * Create a mini-package to run all the commands

# ### Create conda environment
#
#
# Here we will use conda to create an environment for this script. The requirements and specific versions to successfully run this script are listed below. 

# +
# # Create a new conda environment
# git clone https://github.com/FedeGerva/pydeseq2
# conda create --name pyDESeq2_env --file ~/pydeseq2/conda_env/conda_requirements.txt

# #Open conda env
# conda activate pyDESeq2_env

# pip install -r ~/pydeseq2/conda_env/1.4.2/pip_requirements.txt #impo rpy2 > 3.2.5 (see https://github.com/rpy2/rpy2/issues/631)

# #install pydeseq2 package
# pip install -e ~/pydeseq2

# R --vanilla <<code
# chooseCRANmirror(graphics=FALSE, ind=88) #Italy GARR

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
    
# BiocManager::install("GenomeInfoDbData") #version 1.2.2
# BiocManager::install("pasilla") #version 1.14.0
# BiocManager::install("vsn") #version 3.54.0
# code

# ipython kernel install --user --name pyDESeq2_env --display-name "pyDESeq2_env"
# -

# ### Import Python and R packages
#
# The first step is to import all the packages and libraries that we need to work with Python and R. Below you can see two chunks, one to import python packages and rpy2 package, the second one to import all the R packages and library. To use rpy2/R you need to import first rpy2 packages and then all the R libraries. To import R vectors and functions you can use the [R object package](https://rpy.sourceforge.io/rpy2/doc-dev/html/robjects.html#module-rpy2.robjects), please notice that “.” is syntactically valid for R objects but not for python’s. For this reason, “.” of R function are substituted with “_”. The following library are the one essential to run the script, keep in consideration the possibility to add more.
#

# #### Import python modules

# +
# from numpy.random import multinomial, random
import numpy as np
from numpy import *
from matplotlib.ticker import FuncFormatter
from matplotlib import collections  as mc
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import pandas as pd
import os.path
import os

##if you're on JupyterHub you may need to specify the path to R - do it before import any rpy2 modules
jupiter_dir="/jupyterminiconda3" #name of the directory in which conda is saved
from pathlib import Path
home = str(Path.home())
os.environ['R_HOME']=os.path.join(home + jupiter_dir + "/envs/pyDESeq2_env/lib/R/")

##Import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
robjects.numpy2ri.activate()
# %load_ext rpy2.ipython
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula

##Import pydeseq2 package
import pydeseq2
# -

# #### Import R packages

# +
##Import R packages
deseq = importr('DESeq2') #import deseq2
pasilla = importr('pasilla') #import pasilla
vsn = importr('vsn') #import vsn

##Import R functions
assay = robjects.r['assay']
assays=robjects.r['assays']
coef=robjects.r['coef']
colData=robjects.r['colData']
colnames=robjects.r['colnames']
counts=robjects.r['counts']
dispersions=robjects.r['dispersions']
mcols=robjects.r['mcols']
metadata = robjects.r['metadata']
names=robjects.r['names']
normTransform=robjects.r['normTransform']
rownames=robjects.r['rownames']
sign = robjects.r['sign']
sizeFactors=robjects.r['sizeFactors']
system_file = robjects.r['system.file']
substr=robjects.r['substr']
summary = robjects.r['summary']
# -

##Initialize/setting working directory
wd=os.path.join(home, "demo_pyDESeq2")  #declare working dir
os.chdir(wd) #set working directory

# ### Set options
#
# The following chunks are made to set all the visualization options that you want to add. I also create a dictionary with the annotation variable and colours to apply in each plot.

##Display df option
display_max_rows=pd.get_option('display.max_rows')
display_max_columns=pd.get_option('display.max_columns')
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', display_max_rows)
pd.set_option('display.max_columns', display_max_columns)
##Plot options
plt.style.use('default')
sns.set_style("white")

##Create a dict with all the condition and the colors that you want to use - you will need to choose the colors and attributes to include inside the dict
TCcmap={'untreated:single-read': '#66c2a5','untreated:paired-end': '#fc8d62', 'treated:single-read':'#8da0cb', 'treated:paired-end':'#e78ac3', #pca
        'treated': '#9a009a','untreated': '#009900', 
        'single-read':'#110AEA','paired-end':'#cc0000',
       } 


# ### Import count matrix and samples information
#
# The toy dataset chosen to run this tutorial is the count data from [Pasilla](http://bioconductor.org/packages/release/data/experiment/html/pasilla.html) package since it was the one chose by DESeq2 tutorial, making it easier to compare and follow.

# +
##definition of functions to import objects. These functions are not in the package to be more flexible
def importCts(pathCts):
    """imort count matrix as pandas df"""
    cts=pd.read_table(pathCts[0], sep="\t", header=0, index_col=0, dtype="str")
    cts=cts.apply(pd.to_numeric)
    return cts

def importColdata(pathColdata, pathCts):
    """imort sample annotation as pandas df"""
    coldata=pd.read_csv(pathColdata[0],index_col=0) #import coldata information
    coldata.index=coldata.index.str.replace(r'fb$', '')
    coldata=coldata.loc[importCts(pathCts).columns]
    return coldata


# -

cts=importCts(system_file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=True))
coldata=importColdata(system_file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=True), system_file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=True))

##Print head of the count table
pd.DataFrame(cts).head()

pd.DataFrame(coldata)

# +
coldata.index=coldata.index.str.replace(r'fb$', '') #remove fb to have the same names in coldata and cts

coldata=coldata.loc[cts.columns] #reorder coldata according to cts order

coldata.index.equals(cts.columns) #check that the order and names are the same

# +
##set name of the experiment
name='pyDESeq2_demo'

##set design formula
condition = Formula("~condition")
# sampleinfo.head()

##set name to save rds files
file_dds = os.path.join(wd, "dds_"+ name+".rds")
file_rld = os.path.join(wd, "rld_"+ name+ ".rds")
file_vsd = os.path.join(wd, "vst_"+ name+ ".rds")
# -

##Pre-filtering - not mandatory
keep=np.sum(cts,axis=1) >= 10
cts=cts[keep]

##convert cts and coldata from pandas to r dataframe - you can also pass numpy array and skip this process
r_cts, r_coldata=pydeseq2.pyConvertPandas(cts, coldata)

# ### Create or load dds/rld/vsd files
#
# In the next session, the ` dds `, ` rld `, ` vsd ` are created if they don’t exist otherwise are simply loaded. To create these objects, we need to prepare the count table(`cts`) and ` coldata` in the proper format and to set the design formula.

dds=pydeseq2.pyDESeqDataSetFromMatrix(r_cts, r_coldata, design=condition)
rld=pydeseq2.pyrlogTransformation(dds)
vsd=pydeseq2.pyrVarStabTionformation(dds)

print(dds)

# ### Size factor estimation
#
# The size factor is calculated using _sizeFactors_ function of DESeq2 and the correlation between size factor and library size is shown in the associated plot.

print('size factors: ')
print(sizeFactors(dds))

# +
#Create df with libsize info and annotation
libSize_df=pydeseq2.pyCreateLibSizedf(dds, coldata=coldata)

#Plot
pydeseq2.pyPlotLibSizeFact(libSize_df, save='libsize.png')
# -

# ### Dispersion plot and fitting alternatives
#
# The _plotDispEsts_ and _pymeanSdPlot_ applied using rpy2 automatically generate and save the plot. For this reason, the plot is saved as png and it is shown in the notebook using the module display of IPhyton function.

pydeseq2.pyPlotDispEsts(dds, wd=wd)

# ### Effects of transformations on the variance

# +
##this gives log2(n + 1)
ntd=normTransform(dds)
ntd_df=pd.DataFrame(np.matrix(assay(ntd)), index=rownames(dds), columns=colnames(dds))

dds_obj=[ntd, dds, rld]

for i in dds_obj:
    pydeseq2.pymeanSdPlot(i, wd=wd)
# -

# ### PCA plot and Clusterin
#
# PCA and clustering plot are made combining R and python functions. 

pydeseq2.pyPlotPCA(rld, intgroup_name=['type', 'condition'])

# ### Clustering

pydeseq2.pyPlotClustering(rld, coldata, intgroup_name=['type'] )

# ### Differential analysis
#
# The analysis is executed as described in the reference manual. The contrast must be passed as np.array. The following plots were produced extracting the data and re implementing the plot in python.

##Inspect the comparison already present
print(deseq.resultsNames(dds))

# +
contrast=np.array(['condition','treated', 'untreated']) #set the constrast

res, res_df=pydeseq2.pyresults(dds, contrast=contrast, alpha=0.1, save=True)
# -

print(res) #print the res table

res_df['padjxFC']=round(-np.log10(res_df['padj'])*np.sign(res_df['log2FoldChange']), 3) #add column with the padj * sign of log2FC - not mandatory
res_df.head()

# +
##Extract results
n_padj=0.1 #set padj treshold

print(res_df.sort_values('padj', ascending=True).head())

print()

print('summary of the res:')
summary(res)

print()

print('Total number of differential interval')
print(res_df[res_df['padj'] <= n_padj].dropna().count()[1])
# -

pydeseq2.pyPlotMeannormCount_lg10pval(res_df)

# ### Metadata

# +
##Metadata - alpha set
print('alpha: ')
print(metadata(res)[4])

##Percentage of genes set to NA
##BaseMean threshold of genes set to NA 
print('filterThreshold: ')
print(metadata(res)[0])
# -

pydeseq2.pyPlotFilterNumRej(res)

pydeseq2.pyPlotHistpvalue(res, res_df)

# ### MAplot

pydeseq2.pyPlotMA(res_df, n_padj=0.1, ylim=4)

# ### Volcano plot

pydeseq2.pyPlotVolcano(res_df,  n_padj=0.1, ylim_max=40, xlim=2)

# +
n_log2FC=0

##Retrive number of gene/interval using the treshold set for the vulcano plot
diff_up=res_df[(res_df['log2FoldChange'] >= n_log2FC) & (res_df['padj'] <= n_padj)]
print('#of differentially upregulated')
print(shape(diff_up)[0])

diff_down=res_df[(res_df['log2FoldChange'] <= -n_log2FC) & (res_df['padj'] <= n_padj)]
print('#of differentially downregulated')
print(shape(diff_down)[0])
# -

# ### Boxplot of differentially upregulated/downregulated of genes-intervals - normalized counts

norm_filt_melt[['gene_id','variable', 'value','condition', 'type']]

# +
norm_count=counts(dds, normalized=True) #extract normalized count matrix 
norm_count=pd.DataFrame(np.matrix(norm_count), index=rownames(dds), columns=colnames(dds))

norm_filt=norm_count.loc[diff_up.index,] #names of the interval 2keep - if you want to have downregulated you just need to change de variable
col2keep=norm_filt.columns #names of the samples

norm_filt['gene_id']=norm_filt.index
norm_filt_melt=pd.melt(norm_filt, id_vars='gene_id', value_vars=col2keep)
coldata['variable']=coldata.index
norm_filt_melt=pd.merge(norm_filt_melt, coldata,  on='variable')

fig, (boxPlot_counts) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

##plot
norm_filt_melt[['value','condition']].boxplot(by=['condition'], ax=boxPlot_counts, grid=False, notch=True, showfliers=False, patch_artist=True)

##axis
boxPlot_counts.set_title('')
boxPlot_counts.get_figure().suptitle("")
boxPlot_counts.set_xlabel('category', size=10)
boxPlot_counts.set_ylabel('normalized counts', size=10)
boxPlot_counts.tick_params(axis='both', which='major', labelsize=5)
# -

# ### Boxplot of differentially upregulated/downregulated of genes-intervals - mean of normalized count

# +
norm_filt_melt_mean=norm_filt_melt.groupby(['variable', 'condition'], as_index=False)['value'].mean()

fig, (scatter_meanNorm) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

##plot
scatter_meanNorm.scatter(norm_filt_melt_mean.condition, norm_filt_melt_mean['value'], s=20, facecolors='none', 
                         edgecolors='black', linewidth=0.5, c=[TCcmap.get(x,"No_key") for x in norm_filt_melt_mean['condition']])

##axis
scatter_meanNorm.set_xlabel('condition', size=10)
scatter_meanNorm.set_ylabel('mean of normalized count', size=10)
scatter_meanNorm.tick_params(axis='both', which='major', labelsize=5)
# -

# ### Heatmap of the count matrix

# +
select=norm_count.mean(axis=1).sort_values(ascending=False)[0:20].index.to_numpy() ##select genes for the heatmap
col_ord= ['treated1','treated2','treated3', 'untreated1','untreated2','untreated3','untreated4'] ##columns order
ntd_df_sel=ntd_df.loc[select][col_ord] ##select row and reorder columns

##plot
g=sns.clustermap(ntd_df_sel, row_cluster=False, col_cluster=False, col_colors=[[TCcmap.get(x,"No_key") for x in coldata.reindex(np.array(ntd_df_sel.columns))['type']],
                 [TCcmap.get(x,"No_key") for x in coldata.reindex(np.array(ntd_df_sel.columns))['condition']]], 
                 cmap='coolwarm', yticklabels=True, xticklabels=True, linewidth=0, annot_kws={"size": 10})

##axis
for label in coldata.reindex(np.array(ntd_df_sel.columns))['type'].unique():
    g.ax_col_dendrogram.bar(0, 0, color=TCcmap[label],
                            label=label, linewidth=0)
    
for label in coldata.reindex(np.array(ntd_df_sel.columns))['condition'].unique():
    g.ax_col_dendrogram.bar(0, 0, color=TCcmap[label],
                            label=label, linewidth=0)

g.ax_col_dendrogram.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.8), frameon=False) ##legend label position
g.cax.set_position([1.02, .30, .03, .2]) ##colorbar position
# -

# ### Other Plots

# #### Denisty plot normalized count

# +
fig, (dens_norm) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

sns.kdeplot(norm_filt_melt.loc[norm_filt_melt['condition']=='treated',]['value'].values, shade=True, color="r", ax=dens_norm, linewidth=1, cut=0)
sns.kdeplot(norm_filt_melt.loc[norm_filt_melt['condition']=='untreated',]['value'].values, shade=True, color="b", ax=dens_norm, linewidth=1, cut=0)
dens_norm.set_xlim(0, 10000)
# -

# ### Investigate the dds

tuple(assays(dds).slotnames())

print(dispersions(dds))

print(coef(dds))

print(substr(names(mcols(dds)),1,10) )

# +
#Work with S4 object
#S4 object: https://rpy2.github.io/doc/v2.9.x/html/generated_rst/s4class.html

from rpy2.robjects.methods import RS4
class DESeqDataSet(RS4):
    pass

dds_myclass = DESeqDataSet(assays(dds))
dds_myclass
# -

type(dds_myclass)

tuple(dds_myclass.slotnames())


