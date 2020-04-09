========
pydeseq2
========
Description of function to run DESeq2 in python

Install
-------
In order to install the development package use `pip`::

#    git clone https://github.com/FedeGerva/pydeseq2.git 
#    cd pydeseq2
#    pip install -e ~/pydeseq2

Then you can import the library as usual::

    import pydeseq2

===========

Nowadays, it is becoming more clear the fundamental role of flexibility between programming languages in bioinformatics analysis. In fact, there are new packages that allow to call functions written in one language from another language. One of the most know package is rpy2(https://rpy.sourceforge.io/rpy2/doc-dev/html/overview.html) which allow to use R function and call R object directly in a Python environment. These packages are of major importance because they allow the users to reuse existing code and perform analysis integrating different scripting languages.

The bioinformatics community is growing and evolving day by day searching for new solutions that guarantee reproducibility, scalability and flexibility. A huge amount of bioinformatics tools developed in past are written in R language. However, in the last few years, python is emerged as a new promising languages for bioinformatics analysis, since it reduces the time and resources of the analysis. One of the most used bioinformatics packages for the NGS analysis (including RNA-seq, ChIP-seq, ATAC-seq, Hi-C etc) is DESeq2(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). This package is part of the [Bioconductor project](https://www.bioconductor.org/) and it’s written in R language.

Here I created an example of the DESeq2 analysis integrating R with python through rpy2 package.

The toy dataset, the methods and plots are a reimplementation of “Analyzing RNA-seq data with DESeq2” (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat)

Note
====

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
