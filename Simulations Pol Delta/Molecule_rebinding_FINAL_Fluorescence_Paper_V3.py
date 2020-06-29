import pandas as pd
import numpy as npy
import random as rnd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.stats import expon
import math as mt
from scipy.io import savemat
#Version works!
#Allows for moving replisomes
#units are 10nm and seconds

#This version calculates track duration while including gaps

size_window = 400
 #setup a window

center_sphere = size_window/2
num_simulations = 40
nuclear_rad = 65
replisome_rad = 15
#num_replisomes = 400

time_steps = 4000
dt = 0.05
Tbind = 20
Tbleach = 23
#D = 50000 #estimate this based on size

#D_rep = 5

P_bind = expon.cdf(x=dt, scale=Tbind)
P_bleach = expon.cdf(x=dt, scale=Tbleach)
#P_BF_BF = 1 - P_bind - P_bleach
P_BD = 1 - P_bind
#P_UF_UF = 1 - P_bleach
P_F = 1 - P_bleach
#P_BD = 1 - P_bind

BF_state = 0
UF_state = 1
BD_state = 2
UD_state = 3
#replisomes_array = npy.linspace(300, 600, 2, endpoint=True, dtype = int)
#print(replisomes_array)

excess_copy_array = npy.linspace(3000, 6000, 2, endpoint=True,dtype=int)
print(excess_copy_array)
sim_array_frac = npy.arange(0.25, 1.00, 0.50)

print(sim_array_frac)
sim_array_D_rep = npy.logspace(1, 2, 2,  endpoint=True)
sim_array_D_rep = sim_array_D_rep*5
print(sim_array_D_rep)
sim_array_D = npy.logspace(3, 4, 2, endpoint = True)
sim_array_D = sim_array_D *5
print(sim_array_D)
#sim_array_len = len(sim_array)
#for num_replisomes in replisomes_array:
#for Excess_copies in excess_copy_array:
alpha_d = 0.5
for bound_frac in sim_array_frac:
    for D_rep in sim_array_D_rep:
        for D in sim_array_D:
            #bound_frac = m
            #D = p

            Excess_copies = 1600
            num_replisomes = 300
            copy_num = num_replisomes + Excess_copies
            var = 2 * D * dt
            var_rep = 2 * D_rep * (dt**alpha_d)
            print("num_replisomes",num_replisomes)
            print("Excess Copies", Excess_copies)
            print("Bound Fraction", bound_frac)
            print("D",D)
            print("D_rep",D_rep)
            Timelapse_name =  "Rebind_V3" + "_" + "Rep" + str(num_replisomes) + "_" + "NucRad" + str(nuclear_rad) + "_" + "Excess" + \
                              str(Excess_copies) + "_" + "D_POL" + str(D) + "_" + "D_rep" + str(D_rep) + "_" + \
                              "BD_Frac"+ str(bound_frac) + "dt" + str(dt) + "alpha" + str(alpha_d)
            #BD_state = 0


            #Transition Probability Matrix
            #Transition_mat = npy.array([[P_BF_BF, P_bind, P_bleach, 0], [0, P_UF_UF, 0, P_bleach], [0, 0, P_BD, P_bind], [0, 0, 0, 1]])
            Transition_mat_D = npy.array([P_BD, P_bind])
            Transition_mat_F = npy.array([P_F, P_bleach])
            #print(Transition_mat_D)
            #print(Transition_mat_F)

            Multi_spots = npy.zeros((num_simulations,1))
            for k in range(num_simulations):
                frame_window = npy.zeros((size_window, size_window, size_window))
                replisome_positions = npy.zeros((num_replisomes, 3, time_steps))
                replisome_state = npy.zeros((num_replisomes,1))
                #Initialize with fraction replisomes being bound
                bound_replisomes = npy.int_(num_replisomes * bound_frac)
                replisome_state[:bound_replisomes] = 1
                BF_num_times = npy.zeros((time_steps, 1))
                BF_num_times[0] = 1
                state_track = npy.zeros((time_steps, num_replisomes))
                tracks_array = npy.zeros((copy_num, 3, time_steps))
                delta_array = npy.zeros((copy_num, 3, 2))
                bd_frequencies = npy.zeros((num_replisomes, 1))

                # R^2 = (x-x0)^2 + (y-y0)^2 + (z-zo)^2
                #MIght have to fix this issue
                for i in range(bound_replisomes):
                        xrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                        yrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                        zrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                        while ((xrep-center_sphere)**2 + (yrep-center_sphere)**2 + (zrep-center_sphere)**2 - replisome_rad**2 >= 0) or \
                                frame_window[xrep, yrep, zrep] == 2:
                            xrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                            yrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                            zrep = rnd.randrange(center_sphere-replisome_rad, center_sphere + replisome_rad, 1)
                        replisome_positions[i, :, 0] = [xrep, yrep, zrep]
                        frame_window[xrep, yrep, zrep] = 2
                for i in range(bound_replisomes, num_replisomes):
                    xrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                    yrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                    zrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                    while ((xrep - center_sphere) ** 2 + (yrep - center_sphere) ** 2 + (
                            zrep - center_sphere) ** 2 - replisome_rad ** 2 >= 0) or \
                            frame_window[xrep, yrep, zrep] == 2 or frame_window[xrep, yrep, zrep] == 1:
                        xrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                        yrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                        zrep = rnd.randrange(center_sphere - replisome_rad, center_sphere + replisome_rad, 1)
                    replisome_positions[i, :, 0] = [xrep, yrep, zrep]
                    frame_window[xrep, yrep, zrep] = 1


                delta_array[:bound_replisomes, :, 0] = replisome_positions[:bound_replisomes,:,0]
                delta_array[1:bound_replisomes, 0, 1] = BD_state

                delta_array[0:bound_replisomes, 1, 1] = npy.int_(npy.arange(1,bound_replisomes + 1))
                #print(npy.int_(npy.arange(1,num_replisomes + 1)))
                delta_array[bound_replisomes:copy_num, 1, 1] = 0
                # Pick one bound molecule to be fluorescent
                delta_array [0, 0, 1] = BF_state
                state_track[0, :bound_replisomes] = 1





                for i in range(bound_replisomes, copy_num):
                    x_pos = rnd.randrange(center_sphere-nuclear_rad, center_sphere + nuclear_rad, 1)
                    y_pos = rnd.randrange(center_sphere-nuclear_rad, center_sphere + nuclear_rad, 1)
                    z_pos = rnd.randrange(center_sphere-nuclear_rad, center_sphere + nuclear_rad, 1)
                    while (x_pos - center_sphere) ** 2 + (y_pos - center_sphere) ** 2 + (z_pos - center_sphere) ** 2 - \
                           nuclear_rad ** 2 >= 0: #and frame_window[x_pos, y_pos, z_pos] == 2:

                        x_pos = rnd.randrange(center_sphere - nuclear_rad, center_sphere + nuclear_rad, 1)
                        y_pos = rnd.randrange(center_sphere - nuclear_rad, center_sphere + nuclear_rad, 1)
                        z_pos = rnd.randrange(center_sphere - nuclear_rad, center_sphere + nuclear_rad, 1)
                    delta_array[i, :, 0] = [x_pos, y_pos, z_pos]
                    delta_array[i, 0, 1] = UD_state
                #delta_array_3 = delta_array[:, 0, 1]



                #ax.scatter(delta_array[:, 0, 0], delta_array[:, 1, 0], delta_array[:, 2, 0], c='g')



                #print([delta_array[:, :, 0]])
                tracks_array[:, :, 0] = delta_array[:, :, 0]
                delta_array_2 = delta_array[:, :, 0]
                print("Initialization Complete")
                for i in range(time_steps):
                    print("Time Step", i)
                    if i == 0:
                        #fig = plt.figure()
                        #ax = fig.add_subplot(111, projection='3d')
                        #ax.scatter(replisome_positions[:, 0, i], replisome_positions[:, 1, i], replisome_positions[:, 2, i], c='b')
                        filename_save = Timelapse_name + "_" + str(i) + ".png"
                        BF_molecules = npy.argwhere(delta_array[:, 0, 1] == BF_state)
                        UF_molecules = npy.argwhere(delta_array[:, 0, 1] == UF_state)
                        BD_molecules = npy.argwhere(delta_array[:, 0, 1] == BD_state)
                        UD_molecules = npy.argwhere(delta_array[:, 0, 1] == UD_state)
                        #ax.scatter(tracks_array[BF_molecules, 0, i], tracks_array[BF_molecules, 1, i],
                                  # tracks_array[BF_molecules, 2, i], c='g')
                        #ax.scatter(tracks_array[UF_molecules, 0, i], tracks_array[UF_molecules, 1, i],
                                   #tracks_array[UF_molecules, 2, i], c='r')
                        #ax.scatter(tracks_array[BD_molecules, 0, i], tracks_array[BD_molecules, 1, i],
                                   #tracks_array[BD_molecules, 2, i], c='c')
                        #ax.scatter(tracks_array[UD_molecules, 0, i], tracks_array[UD_molecules, 1, i],
                                   #tracks_array[UD_molecules, 2, i], c='y')
                        #plt.xlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                        #plt.xlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                        #ax.set_zlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                        # plt.show()
                        #plt.savefig(filename_save)
                        #plt.close()
                        continue
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111, projection='3d')

                    for j in range(num_replisomes):
                        x_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))
                        y_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))
                        z_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))
                        replisome_positions[j, 0, i] = replisome_positions[j, 0, i - 1] + x_step_rep
                        replisome_positions[j, 1, i] = replisome_positions[j, 1, i - 1] + y_step_rep
                        replisome_positions[j, 2, i] = replisome_positions[j, 2, i - 1] + z_step_rep
                        x_coord_rep = npy.int_(replisome_positions[j, 0, i])
                        y_coord_rep = npy.int_(replisome_positions[j, 1, i])
                        z_coord_rep = npy.int_(replisome_positions[j, 2, i])
                        while (x_coord_rep>size_window or y_coord_rep>size_window or z_coord_rep> size_window) or \
                                ((replisome_positions[j, 0, i] - center_sphere) ** 2 + (replisome_positions[j, 1, i] - center_sphere) ** 2 + \
                                (replisome_positions[j, 2, i] - center_sphere) ** 2 - nuclear_rad ** 2 >= 0) or \
                                frame_window[x_coord_rep, y_coord_rep, z_coord_rep] == 1 or \
                                frame_window [x_coord_rep,y_coord_rep, z_coord_rep] == 2:
                            x_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))
                            y_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))
                            z_step_rep = rnd.normalvariate(0, mt.sqrt(var_rep))

                            replisome_positions[j, 0, i] = replisome_positions[j, 0, i - 1] + x_step_rep
                            replisome_positions[j, 1, i] = replisome_positions[j, 1, i - 1] + y_step_rep
                            replisome_positions[j, 2, i] = replisome_positions[j, 2, i - 1] + z_step_rep
                            x_coord_rep = npy.int_(replisome_positions[j, 0, i])
                            y_coord_rep = npy.int_(replisome_positions[j, 1, i])
                            z_coord_rep = npy.int_(replisome_positions[j, 2, i])
                        x_coord_rep_rel = npy.int_(replisome_positions[j, 0, i-1])
                        y_coord_rep_rel = npy.int_(replisome_positions[j, 1, i-1])
                        z_coord_rep_rel = npy.int_(replisome_positions[j, 2, i-1])
                        frame_window[x_coord_rep_rel, y_coord_rep_rel, z_coord_rep_rel] = 0
                        if replisome_state[j] == 0:
                            frame_window[x_coord_rep, y_coord_rep, z_coord_rep] = 1
                        elif replisome_state[j] == 1:
                            frame_window[x_coord_rep, y_coord_rep, z_coord_rep] = 2
                    #test_crit = replisome_positions[3,:,1]
                    #rep_find = npy.argwhere(replisome_positions[:, :, 1] == test_crit)
                    #print(rep_find)
                        #MAYBE CHANGE TRACK POSITION HERE


                    #ax.scatter(replisome_positions[:, 0, i], replisome_positions[:, 1, i], replisome_positions[:, 2, i], c='b')
                    #print("Timestep",i)
                    for j in range (copy_num):
                        #print("Copy Number", j, delta_array[j, 0,  1 ])
                        if delta_array[j, 0, 1] == BF_state:
                            BF_num_times[i] = 1
                            Transition_p_D = Transition_mat_D
                            rand_multi_select_D = npy.random.multinomial(1, Transition_p_D)
                            element_multi_D = npy.nonzero(rand_multi_select_D)
                            state_select_D = npy.sum(element_multi_D)
                            if state_select_D == 0:
                                Transition_p_F = Transition_mat_F
                                rand_multi_select_F = npy.random.multinomial(1, Transition_p_F)
                                element_multi_F = npy.nonzero(rand_multi_select_F)
                                state_select_F = npy.sum(element_multi_F)
                                if state_select_F == 0:
                                    delta_array[j, 0, 1] = BF_state
                                    #state_track[i, j] = BF_state
                                    #delta_array[j, 1, 1] = replisome
                                    replisome_ID = npy.int_(delta_array[j, 1, 1])
                                    #print(replisome_ID)
                                    tracks_array[j, :, i] = replisome_positions[replisome_ID - 1, :, i]
                                    state_track[i, replisome_ID-1] = 1
                                elif state_select_F == 1:
                                    delta_array[j, 0, 1] = BD_state
                                    #state_track[i, j] = BD_state
                                    replisome_ID = npy.int_(delta_array[j, 1, 1])
                                    tracks_array[j, :, i] = replisome_positions[replisome_ID - 1, :, i]
                                    state_track[i, replisome_ID - 1] = 1

                            elif state_select_D == 1:
                                #delta_array[j, 0 , 1] = 1
                                replisome_ID = npy.int_(delta_array[j, 1, 1])
                                replisome_state[replisome_ID - 1] = 0
                                state_track[i, replisome_ID - 1] = 0
                                x_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 0, i])
                                y_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 1, i])
                                z_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 2, i])
                                #rep_rel_find = npy.int_(npy.argwhere(replisome_positions[:,:,i] == [x_coord_rel, y_coord_rel, z_coord_rel]))
                                frame_window[x_coord_rel, y_coord_rel, z_coord_rel] = 1

                                x_step = rnd.normalvariate(0, mt.sqrt(var))
                                y_step = rnd.normalvariate(0, mt.sqrt(var))
                                z_step = rnd.normalvariate(0, mt.sqrt(var))
                                tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                #print(tracks_array[j, 0, i], x_step)
                                tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                #print(tracks_array[j, 1, i], y_step)
                                tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step
                                #print("X", tracks_array[j, 0, i-1], x_step)
                                #print("Y", tracks_array[j, 1, i-1], y_step)
                                #print("Z", tracks_array[j, 2, i-1], z_step)
                                #print(tracks_array[j, 2, i], z_step)
                                x_coord = npy.int_(tracks_array[j, 0, i])
                                y_coord = npy.int_(tracks_array[j, 1, i])
                                z_coord = npy.int_(tracks_array[j, 2, i])
                                while (tracks_array[j, 0, i] - center_sphere) ** 2 + (tracks_array[j, 1, i] - center_sphere) ** 2 + \
                                        (tracks_array[j, 2, i] - center_sphere) ** 2 - nuclear_rad ** 2 >= 0:
                                    x_step = rnd.normalvariate(0, mt.sqrt(var))
                                    y_step = rnd.normalvariate(0, mt.sqrt(var))
                                    z_step = rnd.normalvariate(0, mt.sqrt(var))

                                    tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                    tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                    tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step

                                x_coord = npy.int_(tracks_array[j, 0, i])
                                y_coord = npy.int_(tracks_array[j, 1, i])
                                z_coord = npy.int_(tracks_array[j, 2, i])
                                #print("X", tracks_array[j, 0, i - 1], x_step)
                                #print("Y", tracks_array[j, 1, i-1], y_step)
                                #print("Z", tracks_array[j, 2, i-1], z_step)
                                Transition_p_F = Transition_mat_F
                                rand_multi_select_F = npy.random.multinomial(1, Transition_p_F)
                                element_multi_F = npy.nonzero(rand_multi_select_F)
                                state_select_F = npy.sum(element_multi_F)
                                if state_select_F == 0:
                                    if frame_window[x_coord, y_coord, z_coord] == 1:
                                        delta_array[j, 0, 1] = BF_state
                                        #state_track[i, j] = BF_state
                                        rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                        rep_find = rep_find[0,0]
                                        print(rep_find)
                                        delta_array[j,1,1] = rep_find + 1
                                        frame_window[x_coord, y_coord, z_coord] = 2
                                        replisome_state[rep_find] = 1
                                        state_track[i, rep_find] = 1

                                    elif frame_window[x_coord, y_coord, z_coord] == 0 \
                                            or frame_window[x_coord, y_coord, z_coord] == 2:
                                        delta_array[j, 0, 1] = UF_state
                                        #state_track[i, j] = UF_state
                                        #rep_find = npy.argwhere(
                                            #npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                        #rep_find = rep_find[0, 0]
                                        delta_array[j, 1, 1] = 0
                                        #replisome_state[replisome_ID - 1] = 0
                                        #state_track[i, replisome_ID - 1] = 0
                                    #delta_array[j, 0, 1]] = 0
                                    #tracks_array[j, :, i] = tracks_array[j, :, i - 1]
                                elif state_select_F == 1:
                                    if frame_window[x_coord, y_coord, z_coord] == 1:
                                        delta_array[j, 0, 1] = BD_state
                                        #state_track[i, j] = BD_state
                                        rep_find = npy.argwhere(
                                            npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                        rep_find = rep_find[0, 0]
                                        delta_array[j, 1, 1] = rep_find + 1
                                        frame_window[x_coord, y_coord, z_coord] = 2
                                        replisome_state[rep_find] = 1
                                        state_track[i, rep_find] = 1

                                    elif frame_window[x_coord, y_coord, z_coord] == 0 \
                                            or frame_window[x_coord, y_coord, z_coord] == 2:
                                        delta_array[j, 0, 1] = UD_state
                                        #state_track[i, j] = UD_state
                                        #rep_find = npy.argwhere(
                                            #npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                        #rep_find = rep_find[0, 0]
                                        delta_array[j, 1, 1] = 0
                                        #replisome_state[replisome_ID-1] = 0
                                        #state_track[i, replisome_ID - 1] = 0

                            #elif state_select_D == 2:
                                #delta_array[j, 0, 1] = 2
                                #tracks_array[j, :, i] = tracks_array[j, :, i-1]
                            #print("Finished 1")



                        elif delta_array[j, 0, 1] == UF_state:
                            #Transition_p_D = Transition_mat_D
                            #rand_multi_select_D = npy.random.multinomial(1, Transition_p_D)
                            #element_multi_D = npy.nonzero(rand_multi_select_D)
                            #state_select_D = npy.sum(element_multi_D)
                            #if state_select_D == 1:
                                #delta_array[j, 0, 1] = 1
                            x_step = rnd.normalvariate(0, mt.sqrt(var))
                            y_step = rnd.normalvariate(0, mt.sqrt(var))
                            z_step = rnd.normalvariate(0, mt.sqrt(var))
                            tracks_array[j, 0, i] = tracks_array[j, 0, i-1] + x_step
                            tracks_array[j, 1, i] = tracks_array[j, 1, i-1] + y_step
                            tracks_array[j, 2, i] = tracks_array[j, 2, i-1] + z_step
                                #print("X", tracks_array[j, 0, i-1], x_step)
                                #print("Y", tracks_array[j, 1, i-1], y_step)
                                #print("Z", tracks_array[j, 2, i-1], z_step)
                            x_coord = npy.int_(tracks_array[j, 0, i])
                            y_coord = npy.int_(tracks_array[j, 1, i])
                            z_coord = npy.int_(tracks_array[j, 2, i])
                            while (tracks_array[j, 0, i] - center_sphere)**2 + (tracks_array[j, 1, i] - center_sphere)**2 + \
                                    (tracks_array[j, 2, i] - center_sphere)**2 - nuclear_rad ** 2 >= 0:
                                x_step = rnd.normalvariate(0, mt.sqrt(var))
                                y_step = rnd.normalvariate(0, mt.sqrt(var))
                                z_step = rnd.normalvariate(0, mt.sqrt(var))
                                tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step
                                    #print("X", tracks_array[j, 0, i - 1], x_step)
                                    #print("Y", tracks_array[j, 1, i - 1], y_step)
                                   # print("Z", tracks_array[j, 2, i - 1], z_step)

                            x_coord = npy.int_(tracks_array[j, 0, i])
                            y_coord = npy.int_(tracks_array[j, 1, i])
                            z_coord = npy.int_(tracks_array[j, 2, i])
                            Transition_p_F = Transition_mat_F
                            rand_multi_select_F = npy.random.multinomial(1, Transition_p_F)
                            element_multi_F = npy.nonzero(rand_multi_select_F)
                            state_select_F = npy.sum(element_multi_F)
                            if state_select_F == 0:
                                if frame_window[x_coord, y_coord, z_coord] == 1:
                                    delta_array[j, 0, 1] = BF_state
                                    #state_track[i, j] = BF_state
                                    rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    rep_find = rep_find[0, 0]
                                    delta_array[j, 1, 1] = rep_find + 1
                                    frame_window[x_coord, y_coord, z_coord] = 2
                                    replisome_state[rep_find] = 1
                                    state_track[i, rep_find] = 1
                                elif frame_window[x_coord, y_coord, z_coord] == 0 or frame_window[x_coord, y_coord, z_coord] == 2:
                                    delta_array[j, 0, 1] = UF_state

                                    #rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    #rep_find = rep_find[0, 0]
                                    #delta_array[j, 1, 1] = 0
                                    #replisome_state[replisome_ID] = 0
                            elif state_select_F == 1:
                                if frame_window[x_coord, y_coord, z_coord] == 1:
                                    delta_array[j, 0, 1] = BD_state
                                    #state_track[i, j] = BD_state
                                    rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    rep_find = rep_find[0, 0]
                                    delta_array[j, 1, 1] = rep_find + 1
                                    frame_window[x_coord, y_coord, z_coord] = 2
                                    replisome_state[rep_find] = 1
                                    state_track[i, rep_find] = 1
                                elif frame_window[x_coord, y_coord, z_coord] == 0 or frame_window[x_coord, y_coord, z_coord] == 2:
                                    delta_array[j, 0, 1] = UD_state

                                    #rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    #rep_find = rep_find[0, 0]
                                    #delta_array[j, 1, 1] = 0
                                    #replisome_state[replisome_ID] = 0
                            #print("Finished 2")


                        elif delta_array[j, 0, 1] == BD_state:

                            Transition_p_D = Transition_mat_D
                            rand_multi_select_D = npy.random.multinomial(1, Transition_p_D)
                            element_multi_D = npy.nonzero(rand_multi_select_D)
                            state_select_D = npy.sum(element_multi_D)
                            if state_select_D == 0:
                                delta_array[j, 0, 1] = BD_state

                                replisome_ID = npy.int_(delta_array[j, 1, 1])
                                tracks_array[j, :, i] = replisome_positions[replisome_ID - 1, :, i]
                                state_track[i, replisome_ID - 1] = 1
                            elif state_select_D == 1:
                                #delta_array[j, 0, 1] = 3
                                replisome_ID = npy.int_(delta_array[j, 1, 1])
                                state_track[i, replisome_ID - 1] = 0
                                replisome_state[replisome_ID - 1] = 0
                                x_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 0, i])
                                y_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 1, i])
                                z_coord_rel = npy.int_(replisome_positions[replisome_ID - 1, 2, i])
                                frame_window[x_coord_rel, y_coord_rel, z_coord_rel] = 1
                           
                                x_step = rnd.normalvariate(0, mt.sqrt(var))
                                y_step = rnd.normalvariate(0, mt.sqrt(var))
                                z_step = rnd.normalvariate(0, mt.sqrt(var))
                                tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step
                                #print("X", tracks_array[j, 0, i-1], x_step)
                                #print("Y", tracks_array[j, 1, i-1], y_step)
                                #print("Z", tracks_array[j, 2, i-1], z_step)
                                x_coord = npy.int_(tracks_array[j, 0, i])
                                y_coord = npy.int_(tracks_array[j, 1, i])
                                z_coord = npy.int_(tracks_array[j, 2, i])
                                while (tracks_array[j, 0, i] - center_sphere) ** 2 + (tracks_array[j, 1, i] - center_sphere) ** 2 + \
                                        (tracks_array[j, 2, i] - center_sphere) ** 2 - nuclear_rad ** 2 >= 0:
                                    x_step = rnd.normalvariate(0, mt.sqrt(var))
                                    y_step = rnd.normalvariate(0, mt.sqrt(var))
                                    z_step = rnd.normalvariate(0, mt.sqrt(var))
                                    tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                    tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                    tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step

                                x_coord = npy.int_(tracks_array[j, 0, i])
                                y_coord = npy.int_(tracks_array[j, 1, i])
                                z_coord = npy.int_(tracks_array[j, 2, i])


                                if frame_window[x_coord, y_coord, z_coord] == 1:
                                    delta_array[j, 0, 1] = BD_state

                                    rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    rep_find = rep_find[0, 0]
                                    #print(rep_find)
                                    delta_array[j, 1, 1] = rep_find + 1
                                    frame_window[x_coord, y_coord, z_coord] = 2
                                    replisome_state[rep_find] = 1
                                    state_track[i, rep_find] = 1
                                elif frame_window[x_coord, y_coord, z_coord] == 0 or frame_window[x_coord, y_coord, z_coord] == 2:
                                    delta_array[j, 0, 1] = UD_state

                                    #rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                    #rep_find = rep_find[0, 0]
                                    delta_array[j, 1, 1] = 0
                                    #replisome_state[replisome_ID-1] = 0
                                    #state_track[i, replisome_ID - 1] = 0
                            #print("Finished 3")

                        elif delta_array[j, 0, 1] == UD_state:
                            x_step = rnd.normalvariate(0, mt.sqrt(var))
                            y_step = rnd.normalvariate(0, mt.sqrt(var))
                            z_step = rnd.normalvariate(0, mt.sqrt(var))
                            tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                            tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                            tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step
                            #print("X", tracks_array[j, 0, i-1], x_step)
                            #print("Y", tracks_array[j, 1, i-1], y_step)
                            #print("Z", tracks_array[j, 2, i-1], z_step)
                            x_coord = npy.int_(tracks_array[j, 0, i])
                            y_coord = npy.int_(tracks_array[j, 1, i])
                            z_coord = npy.int_(tracks_array[j, 2, i])
                            while (tracks_array[j, 0, i] - center_sphere) ** 2 + (tracks_array[j, 1, i] - center_sphere) ** 2 + \
                                    (tracks_array[j, 2, i] - center_sphere) ** 2 - nuclear_rad ** 2 >= 0:
                                x_step = rnd.normalvariate(0, mt.sqrt(var))
                                y_step = rnd.normalvariate(0, mt.sqrt(var))
                                z_step = rnd.normalvariate(0, mt.sqrt(var))
                                tracks_array[j, 0, i] = tracks_array[j, 0, i - 1] + x_step
                                tracks_array[j, 1, i] = tracks_array[j, 1, i - 1] + y_step
                                tracks_array[j, 2, i] = tracks_array[j, 2, i - 1] + z_step

                            x_coord = npy.int_(tracks_array[j, 0, i])
                            y_coord = npy.int_(tracks_array[j, 1, i])
                            z_coord = npy.int_(tracks_array[j, 2, i])
                            if frame_window[x_coord, y_coord, z_coord] == 1:
                                delta_array[j, 0, 1] = BD_state

                                rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                rep_find = rep_find[0, 0]
                                #print(rep_find)
                                delta_array[j, 1, 1] = rep_find + 1
                                frame_window[x_coord, y_coord, z_coord] = 2
                                replisome_state[rep_find] = 1
                                state_track[i, rep_find] = 1

                            elif frame_window[x_coord, y_coord, z_coord] == 0 or frame_window[x_coord, y_coord, z_coord] == 2:
                                delta_array[j, 0, 1] = UD_state

                                #rep_find = npy.argwhere(npy.int_(replisome_positions[:, :, i]) == [x_coord, y_coord, z_coord])
                                #rep_find = rep_find[0, 0]
                                #delta_array[j, 1, 1] = 0
                                #replisome_state[rep_find] = 0
                            #print("Finished 4")
                    BF_molecules = npy.argwhere(delta_array[:,0,1] == BF_state)
                    UF_molecules = npy.argwhere(delta_array[:,0,1] == UF_state)
                    BD_molecules = npy.argwhere(delta_array[:,0,1] == BD_state)
                    UD_molecules = npy.argwhere(delta_array[:,0,1] == UD_state)
                    filename_save = Timelapse_name + "_" + str(i) + ".png"
                    #ax.scatter(tracks_array[BF_molecules, 0, i], tracks_array[BF_molecules, 1, i],
                              # tracks_array[BF_molecules, 2, i], c='g')
                    #ax.scatter(tracks_array[UF_molecules, 0, i], tracks_array[UF_molecules, 1, i],
                              # tracks_array[UF_molecules, 2, i], c='r')
                    #ax.scatter(tracks_array[BD_molecules, 0, i], tracks_array[BD_molecules, 1, i],
                             #  tracks_array[BD_molecules, 2, i], c='c')
                    #ax.scatter(tracks_array[UD_molecules, 0, i], tracks_array[UD_molecules, 1, i],
                               #tracks_array[UD_molecules, 2, i], c='y')


                    #plt.xlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                    #plt.xlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                    #ax.set_zlim((center_sphere - nuclear_rad, center_sphere + nuclear_rad))
                    #plt.show()
                    #plt.savefig(filename_save)
                    #plt.close()
                #BF_freq = npy.count_nonzero(BF_num_times,0)
                BF_non_zero_ele = npy.nonzero(BF_num_times) #find nonzero elements
                #print(BF_non_zero_ele)
                BF_end_ele = npy.amax(BF_non_zero_ele)
                BF_freq = BF_end_ele + 1 #number of time steps including gaps
                #print(BF_end_ele)#find the end of the non zero elements
                rev_BF  = BF_num_times[:BF_end_ele + 1,:]
                diff_array = npy.diff(rev_BF, axis=0)
                find_neg = npy.argwhere(diff_array[:, 0] == -1)
                find_pos = npy.argwhere(diff_array[:, 0] == 1)
                num_ele = len(find_neg)
                gaps_array = npy.zeros((num_ele, 1))
                for j in range(num_ele):
                    gap = find_pos[j] - find_neg[j]
                    gaps_array[j] = gap
                #print("REVISED",rev_BF)
                #print(BF_freq)

                for j in range(num_replisomes):
                    freq_bd = npy.count_nonzero(state_track[:,j] ,0)

                    #print(freq_bd)
                    bd_frequencies[j] = freq_bd/time_steps
                bd_frequencies2 = npy.array(bd_frequencies)
                mat_save_name = Timelapse_name + "_" + str(k) + ".mat"
                mat_gaps_name = Timelapse_name + "_" "Gaps" + "_" + str(k) + ".mat"
                savemat(mat_save_name, {'bound':bd_frequencies2} )
                savemat(mat_gaps_name, {'Gaps': gaps_array})
                #print(type(bd_frequencies))
                #print(bd_frequencies)
                #bd_frequencies.astype(int)
                #sns.distplot(bd_frequencies2)
                #plt.show()


                #fig2 = plt.figure()
                #ax1 = fig2.add_subplot(111)

                print("Simulation Finished",k)
                Multi_spots[k] = BF_freq

                #print(BF_num_times)
            Fluorescent_save_name = Timelapse_name + "Fluor" + ".mat"
            savemat(Fluorescent_save_name, {'Fluorescent_Track_Durations':Multi_spots} )
        #sns.distplot(Multi_spots, fit= expon)

        #plt.show()
        #plt.savefig("HistogramFIT.png")
        #print(expon.fit(Multi_spots))
            #print(BF_num_times[0:,1])










