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
CTF4_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Diffusion Plots\CTF4_Diffusion_total_rev.mat")
CTF4_dat = CTF4_load["D_values_tot_rev"]
CTF4_dat = npy.array(CTF4_dat)

#Load Pol3
POL32_load = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Diffusion Plots\POL32_Diffusion_total_rev.mat")
POL32_dat = POL32_load["D_values_tot_rev"]
POL32_dat = npy.array(POL32_dat)
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
Histone_load_tr = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Diffusion Plots\YTK1434_Diffusion_total_rev_tr.mat")
Histone_dat_tr = Histone_load_tr["D_values_tot_rev"]
Histone_dat_tr = npy.array(Histone_dat_tr)

Histone_load_exp = sio.loadmat(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Diffusion Plots\YTK1434_Diffusion_total_rev_exp.mat")
Histone_dat_exp = Histone_load_exp["D_values_tot_rev"]
Histone_dat_exp = npy.array(Histone_dat_exp)
#CDC_dat = pd.DataFrame(CDC_dat)
#MCM4_dat = pd.DataFrame(MCM4_dat)
#POL2_dat = pd.DataFrame(POL2_dat)
POL32_dat = pd.DataFrame(POL32_dat)
#POL3_dat = pd.DataFrame(POL3_dat)
#POL32_dat = pd.DataFrame(POL32_dat)
#POL12_dat = pd.DataFrame(POL12_dat)
#CTF4_dat = pd.DataFrame(CTF4_dat)
Histone_dat_tr = pd.DataFrame(Histone_dat_tr)
Histone_dat_exp = pd.DataFrame(Histone_dat_exp)



Data_total = pd.concat([POL32_dat, Histone_dat_tr, Histone_dat_exp], axis = 1)
Data_total = npy.array(Data_total)
Data_frame = pd.DataFrame(Data_total, columns=['POL32 Prediction','Histone H3 Training', 'Histone H3 Prediction'])
#Data_frame = Data_frame.melt(var_name='Diffusion Coefficients (um^2/s)', value_name='PDF')
print(Data_frame)
#fig, ax = mtl.subplots()
sns.distplot( POL32_dat, color="blue",label = 'POL32 Prediction', kde = True, hist = False)
sns.distplot( Histone_dat_tr, color="red", label = 'Histone H3 Training', kde = True, hist = False)
sns.distplot( Histone_dat_exp, color="green", label = 'Histone H3 Prediction', kde = True, hist = False)
mtl.legend()

#sns.violinplot(x='Protein', y='Track Durations', data = Data_frame, cut = 0, inner = None)
#fig2,ax2 = mtl.subplots()
#sns.barplot(x='Protein', y='Track Durations', data = Data_frame, estimator=Trunc_exp)
mtl.savefig(r"C:\Users\Reyes-Lamothe Lab\Pictures\Microscopy Data\Nitin\PALM Yeast\Diffusion Plots\Diffusion Comparison.ps")
#ax2 = ax.twinx()


mtl.show()
print("Finished")

