#pydeseq2 package
# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

__author__ = """Federica Gervasoni"""
__email__ = 'federica.gervasoni@unimi.it'
__version__ = '0.1.0'

#"""Utils"""
from ._utils import pyConvertPandas, pyDESeqDataSetFromMatrix, pyrlogTransformation, pyrVarStabTionformation, pyCreateLibSizedf, pyresults

#""" Plotting """
from ._plotting import pyPlotLibSizeFact
from ._plotting import pyPlotDispEsts, pymeanSdPlot
from ._plotting import pyPlotPCA, pyPlotClustering, pyPlotMeannormCount_lg10pval, pyPlotFilterNumRej, pyPlotHistpvalue
from ._plotting import pyPlotVolcano, pyPlotMA

