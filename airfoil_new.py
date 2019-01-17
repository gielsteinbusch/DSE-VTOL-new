# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 17:06:32 2019

@author: giel
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 13:01:51 2019

@author: giel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas

#################### IMPORT VALUES

step = 1. / 0.1 #second number is step length (m)
L_HAMRAC = 11.
W_airframe = 1.35
H_airframe = 1.8
H_bottom = 0.3
MTOW = 2006
lf = 4


boom_distance = 0.25
t_sk = 0.002 #m
t_sp = 0.01
d_skids = 2200 #mm
cg = 4170.320885 #mm

#stringer geometry
t_str = 0.0015 #mm
b_str1 = 0.02 #mm
b_str2 = 1.5*b_str1
A_stringer = 3*b_str1*t_str + b_str2*t_str #mm2

## assume alluminium airframe
ve = 0.33
E = 73.1*10**9
sigma_yield = 324*10**6
density = 2.78 #g/cm3


#################### APPLIED LOADS ############################################

## CASE 1: hover
## Forces with their distances to the nose tip
Airframe =       [-200 *9.81 , 4045.]
Control_h =     [-100*9.81 , 8500.]
Control_v =     [-100*9.81 , 11000.]
Avionics =      [-75*9.81 , 500.]
Hoist =         [-50*9.81 , 5150.]
Powertrain =    [-100*9.81, 4988.67]
Engine =        [-280*9.81, 5500.]
Fuel_tank =     [-20*9.81, 5500.]
EMS =           [-130*9.81, 3250.]
Skids =         [-100*9.81, 3250.]

Pilot =         [-85*9.81, 1100.]
Cabin1 =        [-85*9.81, 1700.]
Cabin2 =        [-85*9.81, 4500.]
Fuel =          [-596*9.81 , 5500.]

F_mainrotor =   [MTOW*9.81 , 4988.67]

Forces_case1 = [Airframe , Control_h, Control_v, Avionics, Hoist,
          Powertrain, Engine, Fuel_tank, EMS, Skids, 
          Pilot, Cabin1, Cabin2, Fuel,
          F_mainrotor]

## CASE 2: landing
## Forces with their distances to the nose tip
Vehicle_weight = [-lf*MTOW*9.81,  cg]
Skids1 =        [0.2*-Vehicle_weight[0], cg-0.8*d_skids]
Skids2 =        [0.8*-Vehicle_weight[0], cg+0.2*d_skids]

Forces_case2 = [Vehicle_weight, Skids1, Skids2]

## CASE 3: landing 2.0
Airframe =       [-lf*201.8 *9.81 , 4045.]
Control_h =     [-lf*43.34*9.81 , 8500.]
Control_v =     [-lf*43.34*9.81 , 11000.]
Avionics =      [-lf*136.*9.81 , 500.]
Hoist =         [-lf*50.*9.81 , 5200.]
Powertrain =    [-lf*78.*9.81, 3811.35]
Engine =        [-lf*300.*9.81, 5500.]
Fuel_tank1 =    [-lf*266.5*9.81, 3750.]
Fuel_tank2 =    [-lf*266.5*9.81, 6000.]
EMS =           [-lf*150.*9.81, 3250.]
Skids =         [-lf*100.*9.81, 3250.]

Pilot =         [-lf*85.*9.81, 1100.]
Cabin1 =        [-lf*85.*9.81, 1700.]
Cabin2 =        [-lf*85.*9.81, 4500.]

Skids1 =        [0.2*-Vehicle_weight[0], cg-0.8*d_skids]
Skids2 =        [0.8*-Vehicle_weight[0], cg+0.2*d_skids]

Forces_case3 = [Airframe , Control_h, Control_v, Avionics, Hoist,
          Powertrain, Engine, Fuel_tank1, Fuel_tank2, EMS, Skids, 
          Pilot, Cabin1, Cabin2,
          Skids1, Skids2]


Forces = Forces_case2

for i in Forces: i[1] = i[1]*10**(-3)

momentlist = []
shearlist = []

poslist = []
for pos in np.linspace(0 , L_HAMRAC , int(L_HAMRAC*step)+1):
    shear = 0
    moment = 0
    for i in range(len(Forces)):
        shear += Forces[i][0] * np.heaviside(pos - Forces[i][1], 0.5)
        moment += Forces[i][0] * (pos - Forces[i][1]) * np.heaviside(pos - Forces[i][1], 0.5)
    
    momentlist.append(moment)
    shearlist.append(shear)
    poslist.append(pos)

maxshear = max(max(shearlist),abs(min(shearlist)))
maxmoment = max(max(momentlist),abs(min(momentlist)))


momentlist = []
shearlist = []

poslist = []
for pos in np.linspace(0 , L_HAMRAC , int(L_HAMRAC*step)+1):
    shear = 0
    moment = 0
    for i in range(len(Forces)):
        shear += Forces[i][0] * np.heaviside(pos - Forces[i][1], 0.5)
        moment += Forces[i][0] * (pos - Forces[i][1]) * np.heaviside(pos - Forces[i][1], 0.5)
    
    momentlist.append(moment)
    shearlist.append(shear)
    poslist.append(pos)

maxshear = max(max(shearlist),abs(min(shearlist)))
maxmoment = max(max(momentlist),abs(min(momentlist)))

## creating the profile
a = W_airframe/2
b = H_airframe/2 + 0.2
thetalist = np.linspace(90,0,180)/180*np.pi

x_coorlist = []
y_coorlist = []

r_circle = 0.15*W_airframe
a, b = 2*r_circle, 2*r_circle

x1_coorlist = list(np.linspace(0, 0.5*W_airframe-2*r_circle, 500))
y1_coorlist = list(H_airframe*np.ones(500))
x_coorlist += x1_coorlist
y_coorlist += y1_coorlist

x2_coorlist = []
y2_coorlist = []
for theta in thetalist:
    r = (a*b) / np.sqrt(b**2 * np.cos(theta)**2 + a**2 * np.sin(theta)**2)
    y2_coorlist.append(H_airframe-2*r_circle + r*np.sin(theta))
    x2_coorlist.append(0.5*W_airframe-2*r_circle + r*np.cos(theta))
x_coorlist += x2_coorlist
y_coorlist += y2_coorlist

x3_coorlist = list(0.5*W_airframe*np.ones(500))
y3_coorlist = list(np.linspace(H_airframe-2*r_circle,H_bottom, 500))
x_coorlist += x3_coorlist
y_coorlist += y3_coorlist
    

## adding booms
x1_boomcoor = []
y1_boomcoor = []
xm_boomcoor = []
ym_boomcoor = []
dis = 0
for i in range(len(x_coorlist)-1):
    dis += np.sqrt((x_coorlist[i]-x_coorlist[i+1])**2 + (y_coorlist[i]-y_coorlist[i+1])**2)
    if dis > boom_distance or i == 0:
        x1_boomcoor.append(x_coorlist[i])
        xm_boomcoor.append(-x_coorlist[i])
        y1_boomcoor.append(y_coorlist[i])
        ym_boomcoor.append(y_coorlist[i])
        dis = 0
x1_boomcoor.reverse()
y1_boomcoor.reverse()
x_boomcoor = x1_boomcoor + xm_boomcoor[1:]
y_boomcoor = y1_boomcoor + ym_boomcoor[1:]

x_boomcoor_spar = 2*[0.5*W_airframe] + 2*[1/6*W_airframe] + 2*[-1/6*W_airframe] + 2*[-0.5*W_airframe]
y_boomcoor_spar = 4*[H_bottom, 0] 
x_boomcoor += x_boomcoor_spar
y_boomcoor += y_boomcoor_spar



#calculate centroids
cen_x = 0
Ay = 0
A = 0
for y in y_boomcoor:
    if y > H_bottom:
        A_stringer_skin = A_stringer + t_sk*boom_distance
        Ay += A_stringer_skin*y
        A += A_stringer_skin
    else:
        A_bottom = 2*t_sp*W_airframe + 4*t_sp*H_bottom
        Ay += 0.5*H_bottom * A_bottom
        A += A_bottom
        break
cen_y = Ay/A

# calculating boom areas
boom_area_list = []
d_wspar = W_airframe/3
d_hspar = H_bottom
for i in range(len(x_boomcoor)):
    if i == 0:
        frac1 = (y_boomcoor[-2]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
        A_boom = A_stringer + t_sk*boom_distance/6 * ((2 + frac1) + (2 + frac2))
    elif i == len(x_boomcoor)-8:
        frac1 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac3 = 1
        A_boom = t_sk*boom_distance/6 *(2 + frac1) + t_sp*d_hspar/6 *(2 + frac2) + t_sp*d_wspar/6 *(2 + frac3)
    elif i == len(x_boomcoor)-7 or i == len(x_boomcoor)-1:
        frac1 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = 1
        A_boom = t_sp*d_hspar/6* (2 + frac1) + t_sp*d_wspar/6* (2 + frac2)
    elif i == len(x_boomcoor)-6 or i == len(x_boomcoor)-4:
        frac1 = 1
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac3 = 1
        A_boom = t_sp*d_wspar/6 * ((2 + frac1) + (2 + frac3)) + t_sp*d_hspar/6 * (2 + frac2)
    elif i == len(x_boomcoor)-5 or i == len(x_boomcoor)-3:
        frac1 = 1
        frac2 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac3 = 1   
        A_boom = t_sp*d_wspar/6 * ((2 + frac1) + (2 + frac3)) + t_sp*d_hspar/6 * (2 + frac2)
    elif i == len(x_boomcoor)-2:
        frac1 = (y_boomcoor[0]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac3 = 1 
        A_boom = t_sk*boom_distance/6 * (2 + frac1) + t_sp*d_hspar/6*(2 + frac2) + t_sp*d_wspar/6*(2 + frac3) 
    else:
        frac1 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
        A_boom = A_stringer + t_sk*boom_distance/6 * ((2 + frac1) + (2 + frac2))
    boom_area_list.append(A_boom)

## calculate moment of inertia
Ixx = 0
Iyy = 0 
for i in range(len(x_boomcoor)):
    Ixx += (y_boomcoor[i]-cen_y)**2 * boom_area_list[i]
    Iyy += (x_boomcoor[i]-cen_x)**2 * boom_area_list[i]

## calculate bending moment
sigma_list = []
for i in range(len(x_boomcoor)):
    sigma_list += [maxmoment * (y_boomcoor[i] - cen_y) / Ixx]

## calculate area
A1 = H_bottom*d_wspar
A2 = W_airframe*(H_airframe-r_circle-H_bottom)+ (W_airframe-2*r_circle)*r_circle + 0.5*np.pi*r_circle**2

## calculate shear
## calculate shear stress
q_base_list = []
q_base = 0
q_moment = 0
for i in range(0,len(x_boomcoor)-8):
    q_base += (maxshear/Ixx) * boom_area_list[i] * (y_boomcoor[i] - cen_y)
    print(i, q_base)
    q_moment += -q_base*(x_boomcoor[i+1]-x_boomcoor[i])*(x_boomcoor[i]-cen_x) +\
                q_base*(y_boomcoor[i+1]-y_boomcoor[i])*(y_boomcoor[i]-cen_y)
    q_base_list.append(q_base)

q_base21 = (maxshear/Ixx) * boom_area_list[15] * (y_boomcoor[16] - cen_y)
q_base13 = 
    
    

#q_base_list = []
#q_base = 0
#q_moment = 0
#print('newlist')
#for i in range(len(x_boomcoor)-9,-1,-1):
#    q_base += -(maxshear/Ixx) * boom_area_list[i] * (y_boomcoor[i] - cen_y)
#    print(i, q_base)
#    q_moment += -q_base*(x_boomcoor[i+1]-x_boomcoor[i])*(x_boomcoor[i]-cen_x) +\
#                q_base*(y_boomcoor[i+1]-y_boomcoor[i])*(y_boomcoor[i]-cen_y)
#    q_base_list.append(q_base)
#q_red = -q_moment / (2*A)
#q_tot_list = [x+q_red for x in q_base_list]
#tau_list = [x/t_sk for x in q_tot_list]





plt.axis([-1.5,1.5,-0.5,2])
plt.scatter(x_boomcoor, y_boomcoor)
