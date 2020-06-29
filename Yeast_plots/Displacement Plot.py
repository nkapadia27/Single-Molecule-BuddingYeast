import scipy.io as sio
import matplotlib.pyplot as mtl
import seaborn as sns
import pandas as pd
#from pandas import Series
#from pandas import DataFrame
import numpy as npy
from scipy.stats import expon
from math import exp
from statsmodels.distributions.empirical_distribution import ECDF
Histone_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Displacements\YTK1434_Displacements.mat")
Histone_dat = Histone_load["cum_displacements"]
Histone_dat = npy.array(Histone_dat)

#ecdf = ECDF(Histone_dat, side= 'right')
#x_CDF = npy.linspace(min(Histone_dat), max(Histone_dat))


#ECDF_values = ecdf(x_CDF)
#print(ECDF_values)
fig, axs = mtl.subplots(nrows = 1, ncols = 2)

sns.distplot(Histone_dat, color="green", axlabel= 'Step Size (pixels)',kde_kws=dict(shade=True),  kde = True, hist = False, ax = axs[0])
sns.distplot(Histone_dat, color = 'red', axlabel= 'Step Size (pixels)', hist = False, kde_kws=dict(cumulative=True), ax = axs[1])
mtl.savefig(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Displacements\Displacements.ps")
mtl.show()