# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 13:48:42 2019

@author: archi
"""
import numpy as np 
import matplotlib.pyplot as plt 
from sparwars import *
from simplebeam import *
from mpl_toolkits.mplot3d import Axes3D

#input values                       (MUST MATCH THE SPARWARS and SIMPLE BEAM VALUES)
radius = 4.46
taper = 1
chord_length = 1
inc_angle = 5
twist = 10
disc_steps = 10
nseg = disc_steps -1 
skin_thickness = list(0.04*np.ones(nseg))
V_flight = 50
rpm = 41 * 60/(2*np.pi)
rho = 0.5
CL = 1.2
W_aircraft = 2500
LDratio = 15

#material properities           NEED TO BE CHANGED 
den = 1.58e3
E = 110e9
Uten = 356e6
G = 23e9

#INITIALISE THE FUNCTIONS
blade = Blade_loading(radius, chord_length, taper, skin_thickness, V_flight, rpm, rho, CL, list_x, list_z, LDratio, disc_steps)
blade.lift_distribution()
blade.shear_distribution()
blade.moment_distribution()
blade.profile()
blade.spar_coor()
blade.profile_new()
blade.twist()
blade.center_gravity()
blade.inertia()
blade.centrifugal_force()
blade.area()
blade.bending_stress()
blade.shear_stress()
blade.max_bend()

beam = Simple_Beam()
beam.sections()
beam.lift_distribution(V_flight)
beam.shear_distribution() 
beam.moment_distribution()
beam.profile()
beam.spar_coor(skin_thickness)
beam.profile_new()
beam.twist()
beam.center_gravity(skin_thickness)
beam.inertia(skin_thickness)
beam.area(skin_thickness)
beam.centrifugal_force(skin_thickness)
beam.deflections()
beam.bending_stress()
beam.shear_stress(skin_thickness)
beam.axial()
beam.sigma_total()
beam.von_mises()

#want to loop through the beam and determine a new thickness in each of the beams sections. 
skint_list = [0.04]

for step in range(beam.nseg): 
    first_iteration = True
    skint = 0.04
    skin_thickness = skint_list + list(skint_list[-1]*(np.ones(beam.nseg - step)))
    while first_iteration == True: 
        skin_thickness[step] = skint 
        beam = Simple_Beam()
        beam.sections()
        beam.lift_distribution(V_flight)
        beam.shear_distribution() 
        beam.moment_distribution()
        beam.profile()
        beam.spar_coor(skin_thickness)
        beam.profile_new()
        beam.twist()
        beam.center_gravity(skin_thickness)
        beam.inertia(skin_thickness)
        beam.area(skin_thickness)
        beam.centrifugal_force(skin_thickness)
        beam.deflections()
        beam.bending_stress()
        beam.shear_stress(skin_thickness)
        beam.axial()
        beam.sigma_total()
        beam.von_mises()
        if max(beam.vonmises[step]) >= Uten : 
            #print(max(beam.vonmises[step]),skint)
            skint = skint + 0.001
        if max(beam.vonmises[step]) <  Uten-50e6 and max(beam.vonmises[step]) > 10e6: 
            #print(max(beam.vonmises[step]),skint)
            if skint < 0.002:
                #print(skin_thickness)
                #print(step, max(beam.vonmises[step]), max(beam.stresses[step]), max(beam.tau_list[step]))
                skint_list.append(skint)
                first_iteration = False
            skint = skint - 0.001
        else: 
            skint_list.append(skint)
            #print(skin_thickness)
            #print(step, max(beam.vonmises[step]), max(beam.stresses[step]), max(beam.tau_list[step]))
            first_iteration = False
            
    #beginning the second iterative loop here that can run as many times as it needs 
skint_list = skint_list[1:]
for step in range(beam.nseg):
    iterating = True           
    #print(skint_list)
    while iterating == True: 
        beam = Simple_Beam()
        beam.sections()
        beam.lift_distribution(V_flight)
        beam.shear_distribution() 
        beam.moment_distribution()
        beam.profile()
        beam.spar_coor(skint_list)
        beam.profile_new()
        beam.twist()
        beam.center_gravity(skint_list)
        beam.inertia(skint_list)
        beam.area(skin_thickness)
        beam.centrifugal_force(skint_list)
        beam.deflections()
        beam.bending_stress()
        beam.shear_stress(skint_list)
        beam.axial()
        beam.sigma_total()
        beam.von_mises()
        for themax in beam.vonmises[:step]:
            if max(themax) >= Uten : 
                #print(max(beam.vonmises[step]),skint)
                skint = skint + 0.001
                skint_list[step] = skint 
            if max(themax) <  Uten-50e6 and max(themax) > 10e6: 
                print(max(beam.vonmises[step]),skint)
                skint_list[step] = skint 
                if skint < 0.002:
                    print(skin_thickness)
                    #print(step, max(beam.vonmises[step]), max(beam.stresses[step]), max(beam.tau_list[step]))
                    skint_list[step]= skint
                    iterating = False
                skint = skint - 0.001
                skint_list[step] = skint
            else: 
                skint_list[step] = skint 
                #print(skin_thickness)
                #print(step, max(beam.vonmises[step]), max(beam.stresses[step]), max(beam.tau_list[step]))
                iterating = False
    print(skint_list)

beam = Simple_Beam()
beam.sections()
beam.lift_distribution(V_flight)
beam.shear_distribution() 
beam.moment_distribution()
beam.profile()
beam.spar_coor(skint_list)
beam.profile_new()
beam.twist()
beam.center_gravity(skint_list)
beam.inertia(skint_list)
beam.area(skin_thickness)
beam.centrifugal_force(skint_list)
beam.deflections()
beam.bending_stress()
beam.shear_stress(skint_list)
beam.axial()
beam.sigma_total()
beam.von_mises()


bend = []
for step in range(beam.nseg):
    for x in beam.vonmises[step]:    #change what you want to plot here 
        bend.append(x)
    
bend = np.array(bend)

plot_z = []
for step in range(beam.nseg):                                    #change this to see full blade plot 
    newz = []
    for z in beam.plot_z[step]: 
        z = z + beam.deflist[-1][step]
        newz.append(z)
    plot_z.append(newz)



fig = plt.figure()
ax = fig.add_subplot(111, projection ='3d')    
ax.scatter(beam.plot_x, beam.plot_y, plot_z, c = bend, cmap=plt.jet())
plt.show()