import scipy.io as sio
import matplotlib.pyplot as mtl
import seaborn as sns
import pandas as pd
#from pandas import Series
#from pandas import DataFrame
import numpy as npy
from scipy.stats import expon
from math import exp
trunc_pt = 3
def Trunc_exp (data):
    mu_est = npy.mean(data)
    trunc_est = (mu_est - trunc_pt)
    return (trunc_est)
#Load CDC45
#CDC_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\20s\CDC45.mat")
#CDC_dat = CDC_load["On_time_combined"]
#CDC_dat = npy.array(CDC_dat)
#Load MCM4
#MCM4_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\1s\MCM4.mat")
#MCM4_dat = MCM4_load["On_time_combined"]
#MCM4_dat = npy.array(MCM4_dat)

#Load Pol2
#POL2_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\1s\POL2.mat")
#POL2_dat = POL2_load["On_time_combined"]
#POL2_dat = npy.array(POL2_dat)
#Load DPB4
CLB5_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\1s\CLB5_revised.mat")
CLB5_dat = CLB5_load["On_time_combined"]
CLB5_dat = npy.array(CLB5_dat)

#Load Pol3
#POL3_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\20s\POL3.mat")
#POL3_dat = POL3_load["On_time_combined"]
#POL3_dat = npy.array(POL3_dat)
#Load Pol32
#POL32_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\20s\POL32.mat")
#POL32_dat = POL32_load["On_time_combined"]
#POL32_dat = npy.array(POL32_dat)

#Load POL12
#POL12_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\8s\POL12.mat")
#POL12_dat = POL12_load["On_time_combined"]
#POL12_dat = npy.array(POL12_dat)

#Load CTF4
#CTF4_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\8s\CTF4.mat")
#CTF4_dat = CTF4_load["On_time_combined"]
#CTF4_dat = npy.array(CTF4_dat)


#H3 Histone
Histone_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\1s\YTK1434.mat")
Histone_dat = Histone_load["On_time_combined"]
Histone_dat = npy.array(Histone_dat)

#CDC_dat = pd.DataFrame(CDC_dat)
#MCM4_dat = pd.DataFrame(MCM4_dat)
#POL2_dat = pd.DataFrame(POL2_dat)
CLB5_dat = pd.DataFrame(CLB5_dat)
#POL3_dat = pd.DataFrame(POL3_dat)
#POL32_dat = pd.DataFrame(POL32_dat)
#POL12_dat = pd.DataFrame(POL12_dat)
#CTF4_dat = pd.DataFrame(CTF4_dat)
Histone_dat = pd.DataFrame(Histone_dat)



Data_total = pd.concat([CLB5_dat, Histone_dat], axis = 1)
Data_total = npy.array(Data_total)
Data_frame = pd.DataFrame(Data_total, columns=['CLB5','Histone H3'])
Data_frame = Data_frame.melt(var_name='Protein', value_name='Track Durations')

fig, ax = mtl.subplots()

sns.violinplot(x='Protein', y='Track Durations', data = Data_frame, cut = 0, inner = None)
#fig2,ax2 = mtl.subplots()
sns.barplot(x='Protein', y='Track Durations', data = Data_frame, estimator=Trunc_exp)
mtl.savefig(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Track Durations Combined\1s\Violin Comparison_CLB5.ps")
#ax2 = ax.twinx()


mtl.show()
print("Finished")

