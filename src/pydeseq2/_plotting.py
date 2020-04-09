# -*- coding: utf-8 -*-
def pyPlotLibSizeFact(df, color='black', alpha=0.7, s=30, linewidth=0.5, Dictcmap=None, save=None, dpi=300,  **kwargs):
    """
    Plot a scatter with the correlation between sizeFactor and library size
    
    Parameters
    -------
    input:
    df
        pd.dataFrame with the information os size factor and library size
    color
        name of one column of the coldata. For now is possible to use just 1 condition
    alpha
        The alpha blending value. Value between 0 and 1. Default: 0.7 
    s
        scalar or array-like, shape (n, ), optional. Defautlt: 30
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    Dictcmap
        Dict object with the association between samples and color code
    save
        path to save plot. Default:None
    """
    
    ##library
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.style.use('default')
    
    ##plot
    colorPlot=np.array([color])

    ##plot
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300, gridspec_kw={'width_ratios':[1]})

    if colorPlot != 'black':
        ax.scatter(df['libSize'], df['sizeFactors'], color=[Dictcmap.get(x,"No_key") for x in np.array(df[colorPlot[0]])],
            alpha=alpha, s=s, linewidth=linewidth, **kwargs)

    else:
        ax.scatter(df['libSize'], df['sizeFactors'], color='black',
            alpha=alpha, s=s, linewidth=linewidth, **kwargs)

    ##axis
    ax.set_xlabel('library size', size=10)
    ax.set_ylabel('size factor', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)

    ##add lm info
    m, b = np.polyfit(df['libSize'], df['sizeFactors'], 1) #calculate lm
    ax.plot(df['libSize'], m*df['libSize'] + b, color='red',  linewidth=1) #add line
    
    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()

def pyPlotDispEsts(dds, wd='.', name='dispersion'):
    """
    Dispersion plot
    
    Parameters
    -------
    DESeqDataSet
        DESeqDataSet, the Default: 'dds'. Other options can be for example rld, vst
    wd
        working directory on which it saves the plot
    name
        name to save file. The default plot is saved in png to avoid errors.
    """
    
    import os
    from IPython.display import Image, display
    import rpy2.robjects as robjects
    options=robjects.r['options']
    png=robjects.r['png']
    plotDispEsts=robjects.r['plotDispEsts']
    dev_off=robjects.r['dev.off']
    
    options(device='png')
    png(os.path.join(wd, name+'.png'))
    plotDispEsts(dds)
    dev_off()

    display(Image(filename=os.path.join(wd, name+'.png')))

def pymeanSdPlot(dds, wd='.', name='meanSdPlot'):
    """
    Dispersion plot
    
    Parameters
    -------
    DESeqDataSet
        DESeqDataSet, the Default: 'dds'. Other options can be for example ntd,rld, vst
    wd
        working directory on which it saves the plot
    name
        name to save file. The default plot is saved in png to avoid errors.
    """
    import os
    from IPython.display import Image, display
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    options=robjects.r['options']
    png=robjects.r['png']
    vsn = importr('vsn') #import vsn
    assay = robjects.r['assay']
    meanSdPlot=robjects.r['meanSdPlot']
    dev_off=robjects.r['dev.off']
    
    options(device='png')
    png(os.path.join(wd, name+'.png'))
    meanSdPlot(assay(dds))
    dev_off()

    display(Image(filename=os.path.join(wd, name+'.png')))

def pyPlotPCA(rld, intgroup_name=None, n_top=500, ncompx=1, ncompy=2, alpha=0.7, s=50, linewidth=0.5, Dictcmap=None, save=None, dpi=300, **kwargs):
    """
    Plot PCA implemented in python
    
    Parameters
    -------
    input:
    DESeqDataSet
         DESeqDataSet, the Default: 'rld'.
    introup
        one or a list of condition used to color the PCA
    n_top
        #genes used to perform the PCA plot. Default: 500
    ncompx
        number of component shown on x axis. Default: 1 (PC1)
    ncompy
        number of component shown on y axis. Default: 2 (PC2)
    alpha
        The alpha blending value. Value between 0 and 1. Default: 0.7 
    s
        scalar or array-like, shape (n, ), optional. Defautlt: 50
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    Dictcmap
        Dict object with the association between samples and color code
    save
        path to save plot. Default: None
    """

    import numpy as np
    import matplotlib.patches as mpatches
    from matplotlib.patches import Patch
    import matplotlib.pyplot as plt
    import pandas as pd
    from numpy import transpose
    import rpy2.robjects as robjects
    assay = robjects.r['assay']
    rownames=robjects.r['rownames']
    colnames=robjects.r['colnames']
    rowVars=robjects.r['rowVars']
    prcomp=robjects.r['prcomp']
    colData=robjects.r['colData']
    as_data_frame=robjects.r['as.data.frame']
    
    plt.style.use('default')

    n_top=n_top
    
    if intgroup_name == None:
        raise Exception("You need to set a color condition in intgroup_name field")
    else:
        intgroup=np.array(intgroup_name)

    table=pd.DataFrame(np.matrix(assay(rld)), index=rownames(rld), columns=colnames(rld))
    rv=np.array(rowVars(assay(rld)))
    select=np.argsort(rv)[::-1][range(min(n_top,len(rv)))]

    pca=prcomp(transpose(np.array(table.iloc[select,])))
    percentVar=np.array(pca[0])**2 / sum(np.array(pca[0])**2) ##the percentVar in R is saved as attr of the pca_data

    if (len(intgroup)) > 1:
        intgroup_df=pd.DataFrame(as_data_frame(colData(rld)), index=rownames(colData(rld)), columns=colnames(colData(rld)))[intgroup]
        group=intgroup_df[intgroup].apply(lambda x: ':'.join(x[x.notnull()]), axis = 1)

    else:
        group=np.array(pd.DataFrame(as_data_frame(colData(rld)), index=rownames(colData(rld)), columns=colnames(colData(rld)))[intgroup[0]])

    d={'PC1': np.array(pca[4])[:,ncompx-1], 'PC2':np.array(pca[4])[:,ncompy-1], 'group': group, 'name':colnames(rld)}
    df = pd.DataFrame(data=d)

    ##PlotPCA
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)
    ax.scatter(df['PC1'], df['PC2'], color=[Dictcmap.get(x,"No_key") for x in np.array(df['group'])], 
               alpha=alpha, s=s, linewidth=linewidth, **kwargs)

    ##axis
    #retrieve values from color dictionary and attribute it to corresponding labels
    leg_color=list(set([Dictcmap.get(x,"No_key") for x in np.array(df['group'])]))
    inv_Dictcmap = {v: k for k, v in Dictcmap.items()}
    leg_group=[inv_Dictcmap[k] for k in list(set([Dictcmap.get(x,"No_key") for x in np.array(df['group'])]))]
    sub_Dictcmap=dict((k, Dictcmap[k]) for k in leg_group)

    leg_el = [mpatches.Patch(facecolor = leg_color, edgecolor = "black", label = leg_group, alpha = 0.4) for leg_group, leg_color in sub_Dictcmap.items()]
    ax.legend(handles = leg_el, loc='upper center', bbox_to_anchor=(0.3, 0.8, 0.5, 0.5), frameon=False, fontsize='xx-small')

    ax.set_xlabel(str(round(percentVar[ncompx-1]*100, 2))+ '% PC'+ str(ncompx), size=10)
    ax.set_ylabel(str(round(percentVar[ncompy-1]*100, 2))+ '% PC'+ str(ncompy), size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)
    
    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()

def pyPlotClustering(rld, coldata, intgroup_name=None, Dictcmap=None, save=None, dpi=300, **kwargs):
    """
    Plot dist matrix implemented in python
    
    Parameters
    -------
    input:
    DESeqDataSet
         DESeqDataSet, the Default: 'rld'.
    coldata
        coldata object with samples info
    introup
        one condition used to color the clustering anno
    Dictcmap
        Dict object with the association between samples and color code
    save
        path to save plot. Default: None
    """
    
    import pandas as pd
    import warnings
    import matplotlib.pyplot as plt
    from numpy import transpose
    import numpy as np
    from scipy.spatial import distance_matrix
    import seaborn as sns; sns.set(color_codes=True)
    import rpy2.robjects as robjects
    colData=robjects.r['colData']
    assay = robjects.r['assay']
    rownames=robjects.r['rownames']
    colnames=robjects.r['colnames']
    
    plt.style.use('default')
 
    if intgroup_name == None:
        warnings.warn("Warning: You need to set a color condition in intgroup_name field")
    else:
        intgroup=np.array(intgroup_name)[0]
    
    ##extract rld info
    rld_t_table=pd.DataFrame(transpose(assay(rld)), index=colnames(rld), columns=rownames(rld))
    dist_mtx=pd.DataFrame(distance_matrix(rld_t_table.values, rld_t_table.values), index=colnames(rld), columns=colnames(rld))
 
    if Dictcmap:    
        ##plot
        g=sns.clustermap(dist_mtx, cmap="mako", row_colors=[Dictcmap.get(x,"No_key") for x in coldata.reindex(np.array(dist_mtx.index))[intgroup]],
                      yticklabels=True, xticklabels=True, **kwargs)

        ##axis
        for label in coldata.reindex(np.array(dist_mtx.index))[intgroup].unique():
            g.ax_col_dendrogram.bar(0, 0, color=Dictcmap[label],
                                    label=label, linewidth=0)
        
        g.ax_col_dendrogram.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.4), frameon=False) ##legend label position
 
    else:
        g=sns.clustermap(dist_mtx, cmap="mako",
                  yticklabels=True, xticklabels=True)

    g.cax.set_position([1.02, .35, .03, .2]) ##colorbar position
    
    if save:
        g.savefig(save, dpi=dpi)

def pyPlotMeannormCount_lg10pval(df, alpha=0.7, s=10, linewidth=0.5, save=None, dpi=300, **kwargs):
    """
    Plot −log10(pvalues) for all the geans over the normalized mean counts
    
    Parameters
    -------
    input:
    df
        pd.dataFrame containing the res
    alpha
        The alpha blending value. Value between 0 and 1. Default: 0.7 
    s
        scalar or array-like, shape (n, ), optional. Defautlt: 10
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    save
        path to save plot. Default: None
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.style.use('default')

    ##Scatter plot with basemean and -log10 pvalue
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

    ax.scatter(df['baseMean']+1, -np.log10(df['pvalue']), color='black', edgecolor='black',
             alpha=alpha, s=s, linewidth=linewidth, **kwargs)

    ##axis
    ax.set_xlabel('mean of normalized counts', size=10)
    ax.set_ylabel('-log10(pvalue)', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.set_ylim(0,30)
    ax.set_xscale('log')
    
    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()

def pyPlotFilterNumRej(res, s=10, linewidth=0.5, save=None, dpi=300, **kwargs):
    """
    Plot number of rejections (adjusted p value less than a significance level), over the quantiles of a filter statistic (the mean of normalized counts) 
    
    Parameters
    -------
    input:
    res
        res object
    s
        scalar or array-like, shape (n, ), optional. Defautlt:10
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    save
        path to save plot. Default: None
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    import rpy2.robjects as robjects
    metadata = robjects.r['metadata']
    
    plt.style.use('default')
    
    ##Scatter plot #reject and quantiles of filters
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

    ##plot
    ax.scatter(metadata(res)[2]['theta'], metadata(res)[2]['numRej'], facecolors='none', edgecolors='black', s=s, linewidth=linewidth, **kwargs)
    ax.plot(metadata(res)[3][0], metadata(res)[3][1], color='red', lw=linewidth)

    ##axis
    ax.set_xlabel('quantiles of filter', size=10)
    ax.set_ylabel('#rejections', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)

    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()

def pyPlotHistpvalue(res, res_df, linewidth=0.5, save=None, dpi=300, **kwargs):
    """
    Plot −log10(pvalues) for all the geans over the normalized mean counts
    
    Parameters
    -------
    input:
    res
        res object    
    res_df
        pd.dataFrame containing the res
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    save
        path to save plot. Default: None
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.patches as mpatches
    from matplotlib.patches import Patch
    import rpy2.robjects as robjects
    metadata = robjects.r['metadata']
    
    plt.style.use('default')
    
    ##Scatter plot #reject and quantiles of filters
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

    ##Histogram of p values for all tests. The area shaded in blue indicates the subset of those that pass the filtering, the area in khaki those that do not pass
    use = res_df[res_df['baseMean'] > metadata(res)[0]].dropna()

    ##plot
    ax.hist(res_df[res_df.isin(use)].dropna()['pvalue'], color='powderblue', edgecolor='black', linewidth=linewidth, bins=np.arange(0, 1, 0.05), **kwargs)
    ax.hist(res_df[~res_df.isin(use)].dropna(subset=['pvalue'])['pvalue'], color='khaki', edgecolor='black', linewidth=linewidth, bins=np.arange(0, 1, 0.05), **kwargs)

    ##axis
    ax.set_xlabel('pValue', size=10)
    ax.set_ylabel('frequency', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)

    ##Create custom legend
    custom_leg = [Patch(facecolor='powderblue', edgecolor='black', linewidth=0.5, label='Color Patch'),
                  Patch(facecolor='khaki', edgecolor='black', linewidth=0.5, label='Color Patch')]

    ax.legend(custom_leg, ['passed', 'not passed'], fontsize='xx-small', loc='upper right')

    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()


def pyPlotVolcano(df, n_padj=0.05, n_log2FC=0, alpha=0.7, s=3, linewidth=0.5, ylim_max=None, xlim=None, a_color="#D90416", b_color="#033E8C", c_color="grey", save=None, dpi=300, **kwargs):
    """
    volcanoPlot
    
    Parameters
    -------
    input:
    df
        pd.dataFrame containing the res
    n_padj
        padj threshold to apply. Default:0.05.
    n_log2FC
        log2FC threshold to apply. Default:0
    alpha
        The alpha blending value. Value between 0 and 1. Default: 0.7 
    s
        scalar or array-like, shape (n, ), optional. Defautlt:3
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    ylim_max
        ylim to set. Default:None
    xlim
        xlim to set. Default:None
    a_color, b_color, c_color
        color for the volcano plot. Default are: 'A':"#D90416", 'B':"#033E8C",  'C':"grey"
    save
        path to save plot. Default:None
    """

    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.style.use('default')
        
    Vmap={'A':a_color, 'B':b_color,  'C':c_color}
    
    volcano=df
    volcano['threshold']=volcano.apply(lambda row: 'A' if (row.log2FoldChange > n_log2FC) and (row.padj < n_padj) else ('B' if (row.log2FoldChange < -n_log2FC) & (row.padj < n_padj) else 'C'), axis=1)

    ##plot
    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

    ax.scatter(volcano['log2FoldChange'], -np.log10(volcano['padj']), alpha=alpha, color=[Vmap.get(x,"No_key") for x in volcano['threshold']],
             s=s, linewidth=linewidth, **kwargs) 
    
    ##axis
    ax.set_xlabel('log2 fold change', size=10)
    ax.set_ylabel('-log10 padj', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.axhline(-np.log10(n_padj), ls='--', color='black', lw=.5, alpha=0.5 )
    
    if ylim_max:
        ax.set_ylim(0,ylim_max)
    
    if xlim:
        ax.set_xlim(-xlim,xlim)
        
    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()   


def pyPlotMA(df, n_padj=0.05, alpha=0.7, s=3, linewidth=0.5, ylim=None, save=None, dpi=300, **kwargs):
    """
    MAplot
    
    Parameters
    -------
    input:
    df
        pd.dataFrame containing the res
    n_padj
        padj threshold to fiter. Default: 0.05.
    alpha
        The alpha blending value. Value between 0 and 1. Default: 0.7 
    s
        scalar or array-like, shape (n, ), optional. Defautlt: 3
    linewidth
        scalar or array-lik. The linewidth of the marker edge. Default: 0.5
    ylim
        ylim to set. Default: None
    save
        path to save plot. Default: None
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.style.use('default')

    fig, (ax) = plt.subplots(1,1,figsize=(3, 3),dpi=300)

    col = np.where(df['padj'] < n_padj,'#cd0000', '#525252')

    ##plot
    ax.scatter(df['baseMean'], df['log2FoldChange'], alpha=alpha, color=col,
             s=s, linewidth=linewidth, **kwargs)

    ##axis
    ax.set_xlabel('mean of normalized count', size=10)
    ax.set_ylabel('log2 fold change', size=10)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.axhline(0, color='red', lw=1, alpha=0.5 )
    ax.set_xscale('log')
    
    if ylim:
        ax.set_ylim(-ylim,ylim)

    if save:
        plt.savefig(save, dpi=dpi)
    
    plt.show()
