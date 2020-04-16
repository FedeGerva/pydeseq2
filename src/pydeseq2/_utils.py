def pyDESeqDataSetFromMatrix(cts, coldata, design=None, file_dds='dds.rds', **kwargs):
    """
    Create and import dds object

    Parameters
    -------
    input:
    cts
        count matrix
    coldata
        dataframe with sample information
    design
        design formula
    file_dds
        This field is mandatory. Default: "dds.rds"
    **kwargs
        Argument passed to DESeqDataSetFromMatrix function of DESeq2
    """
    
    import os
    import sys
    import numpy as np
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    deseq = importr('DESeq2') #import deseq2   
    readRDS = robjects.r['readRDS']
    saveRDS = robjects.r['saveRDS']
    
    if isinstance(cts, (rpy2.robjects.vectors.DataFrame, np.ndarray, np.generic)):

        if not os.path.isfile(file_dds):
            print('dds does not exist, it will be created. Keep calm and wait.')
            dds = deseq.DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=design, **kwargs)
            dds = deseq.DESeq(dds)
            saveRDS(object=dds, file=file_dds)
        else:
            print('dds exists and it will be loaded')
            dds=readRDS(file=file_dds)
    
    else:
        raise Exception("count matrix is not an object of class rpy2.robjects.vectors.DataFrame, np.ndarray or np.generic")

    return dds


def pyrlogTransformation(dds, file_rld='rld.rds', blind=True, fitType = "parametric"):
    """
    Apply Regularized Log' Transformation to dds

    Parameters
    -------
    input:
    dds
        DESeqDataSet object
    file_rld
        Name of the file to be saved. Default: "rld.rds"
    blind
        logical, whether to blind the transformation to the experimental design. Default=True 
    fitType
        Default:"parametric"
    """
    import os
    import sys
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.rinterface import SexpS4
    deseq = importr('DESeq2') #import deseq2   
    readRDS = robjects.r['readRDS']
    saveRDS = robjects.r['saveRDS']

    if isinstance(dds, SexpS4):
        
        if not os.path.isfile(file_rld):
            print('rld does not exist, it will be created. Keep calm and wait.')
            rld = deseq.rlog(dds, blind=blind, fitType=fitType)
            saveRDS(object=rld, file=file_rld)
        else:
            print('rld exists and it will be loaded')
            rld=readRDS(file=file_rld)
    else:
        raise Exception("dds is not an object of class rpy2.robjects.methods.RS4")

    return rld


def pyrVarStabTionformation(dds, file_vsd='vsd.rds', blind=True, nsub = 1000, fitType = "parametric"):
    """
    Apply Variance Stabilizing Transformation on dds

    Parameters
    -------
    input:
    dds
        DESeqDataSet object
    file_vsd
        Name of the file to be saved. Default: "vsd.rds"
    blind
        logical, whether to blind the transformation to the experimental design. Default=True 
    nsub
        the number of genes to subset to. Default: 1000
    fitType
        Default:"parametric" 
    """
    import os
    import sys
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.rinterface import SexpS4
    deseq = importr('DESeq2') #import deseq2   
    readRDS = robjects.r['readRDS']
    saveRDS = robjects.r['saveRDS']
    
    if isinstance(dds, SexpS4):

        if not os.path.isfile(file_vsd):
            print('vsd does not exist, it will be created. Keep calm and wait.')
            vsd = deseq.vst(dds, blind=blind, nsub=nsub, fitType=fitType)
            saveRDS(object=vsd, file=file_vsd)
        else:
            print('vsd exists and it will be loaded')
            vsd=readRDS(file=file_vsd)
    
    else:
        raise Exception("dds is not an object of class rpy2.robjects.methods.RS4")

    return vsd


def pyCreateLibSizedf(dds, coldata):
    """
    Create dataframe with size factors and coldata info

    Parameters
    -------
    input:
    dds
        DESeqDataSet object
    coldata
        dataframe with samples annotation
    """
    import numpy as np
    import pandas as pd
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.rinterface import SexpS4
    counts=robjects.r['counts']
    colnames=robjects.r['colnames']
    rownames=robjects.r['rownames']
    sizeFactors=robjects.r['sizeFactors']

    if isinstance(dds, SexpS4):
    
        ##dds row count matrix
        dds_counts=pd.DataFrame(np.matrix(counts(dds)), columns = colnames(dds), index=rownames(dds))

        ##df with info
        data= {'sizeFactors': sizeFactors(dds),
               'libSize': dds_counts.sum().round(2).div(10^6)
              }
        df = pd.DataFrame(data, columns = ['sizeFactors', 'libSize'], index=colnames(dds))
        df=df.merge(coldata, left_index=True, right_index=True)
        
    else:
        raise Exception("dds is not an object of class rpy2.robjects.methods.RS4")

    return df


def pyresults(dds, save=None, wd='.', name='res_table', index_label='gene_id', **kwargs):
    """
    Extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values. 

    Parameters
    -------
    input:
    dds
        DESeqDataSet object
    save
        If True, save res dataframe as csv.
    name
        Name of the csv file.
    wd
        working directory on which it saves the table
    index_label
        name of the index labels. Default: 'gene_id'
    **kwargs
        Argument passed to results function of DESeq2
        
    Returns
    -------
    Returns two object, res (R object) and res_df (pandas DataFrame) with the result table.
    """
    
    import os
    import sys
    import pandas as pd
    import rpy2
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.rinterface import SexpS4
    robjects.numpy2ri.activate()
    deseq = importr('DESeq2') #import deseq2 
    to_dataframe = robjects.r('function(x) data.frame(x)')
    rownames=robjects.r['rownames']
    
    if isinstance(dds, SexpS4):
        
        res=deseq.results(dds, **kwargs)
        res_df=to_dataframe(res)
        res_df=pd.DataFrame(res_df, columns=['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'], index=rownames(dds))

    else:
        raise Exception("dds is not an object of class rpy2.robjects.methods.RS4")
    
    if save:
        res_df.to_csv(os.path.join(wd, name+'.csv'), index_label="gene_id")        
    
    return(res, res_df)
