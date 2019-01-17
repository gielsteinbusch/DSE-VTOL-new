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
MTOW = 2006
lf = 4

#stringer geometry
t_str = 1.5 #mm
b_str1 = 20 #mm
b_str2 = 1.5*b_str1
A_stringer = 3*b_str1*t_str + b_str2*t_str #mm2

boom_distance = 0.5
t_sk = 0.002 #m
d_skids = 2200 #mm
cg = 4170.320885 #mm

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

a = W_airframe/2
b = H_airframe/2 + 0.2
thetalist = np.linspace(90,270,180)/180*np.pi


## creating the profile, only the left side 
x_coorlist = []
y_coorlist = []
for theta in thetalist:
    r = (a*b) / np.sqrt(b**2 * np.cos(theta)**2 + a**2 * np.sin(theta)**2)
    if r*np.sin(theta) < -H_airframe/2 + 0.2: y_coorlist.append(-H_airframe/2 + 0.2)
    else: y_coorlist.append(r*np.sin(theta))
    x_coorlist.append(r*np.cos(theta))

## adding booms, with equal spacing

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

# calculate centroids (cen_x will obviously be 0)
cen_y = np.average(y_boomcoor)   
cen_x = np.average(x_boomcoor)

# calculating boom areas
boom_area_list = []
A_stringer = A_stringer*10**(-6)  #mm2 -> m2
for i in range(len(x_boomcoor)):
    if i == 0:
        frac1 = (y_boomcoor[-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
    elif i == len(x_boomcoor)-1:
        frac1 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[0]-cen_y) / (y_boomcoor[i]-cen_y)
    else:
        frac1 = (y_boomcoor[i-1]-cen_y) / (y_boomcoor[i]-cen_y)
        frac2 = (y_boomcoor[i+1]-cen_y) / (y_boomcoor[i]-cen_y)
    A_boom = A_stringer + t_sk*boom_distance/6 * ((2 + frac1) + (2 + frac2))
#    print(A_boom, A_stringer,t_sk*boom_distance/6 * ((2 + frac1) + (2 + frac2)))
    boom_area_list.append(A_boom)
    
## calculate moment of inertia
Ixx = 0
Iyy = 0 
for i in range(len(x_boomcoor)):
    Ixx += (y_boomcoor[i]-cen_y)**2 * boom_area_list[i]
    Iyy += (x_boomcoor[i]-cen_x)**2 * boom_area_list[i]
#    print(y_boomcoor[i]-cen_y, boom_area_list[i], Ixx )
## calculate bending stress
sigma_list = []
for i in range(len(x_boomcoor)):
    sigma_list += [maxmoment * (y_boomcoor[i] - cen_y) / Ixx]
#    print(maxmoment, (y_boomcoor[i] - cen_y), Ixx)
    
## area calculation
A = 0
for i in range(len(x_boomcoor)):
    ax, az = x_boomcoor[i], y_boomcoor[i]
    if i == len(x_boomcoor)-1:
        bx, bz = x_boomcoor[0], y_boomcoor[0]
    else: 
        bx, bz = x_boomcoor[i+1], y_boomcoor[i+1]
    u = np.array([ax,az])
    v = np.array([bx,bz])
    A += 0.5*np.cross(v,u)
    
## calculate shear stress
q_base_list = []
q_base = 0
q_moment = 0
for i in range(len(x_boomcoor)):
    q_base += -(maxshear/Ixx) * boom_area_list[i] * (y_boomcoor[i] - cen_y)
    print(boom_area_list[i],(y_boomcoor[i] - cen_y),q_base)
    if i == len(x_boomcoor)-1:
        q_moment += -q_base*(x_boomcoor[0]-x_boomcoor[i])*(x_boomcoor[i]-cen_x) +\
                    q_base*(y_boomcoor[0]-y_boomcoor[i])*(y_boomcoor[i]-cen_y)
    else:
        q_moment += -q_base*(x_boomcoor[i+1]-x_boomcoor[i])*(x_boomcoor[i]-cen_x) +\
                    q_base*(y_boomcoor[i+1]-y_boomcoor[i])*(y_boomcoor[i]-cen_y)
    q_base_list.append(q_base)
q_red = -q_moment / (2*A) 
q_tot_list = [x+q_red for x in q_base_list]
tau_list = [x/t_sk for x in q_tot_list]

## calculate von mises stress
vonmises_list = []
for i in range(len(sigma_list)):
    vonmises_list.append(np.sqrt(sigma_list[i]**2+3*tau_list[i]**2))


maxsigma = max(max(sigma_list),abs(min(sigma_list)))
maxtau = max(max(tau_list),abs(min(tau_list)))
maxvonmises = max(max(vonmises_list),abs(min(vonmises_list)))
##################################### BUCKLING ################################

A_str = 3*b_str1*t_str + b_str2*t_str
A_sk = t_sk*boom_distance

## calculate the critical buckling stress of skin
kc = 4
sigma_cr_sk = np.pi**2 * kc * E * (t_sk/boom_distance)**2 / (12*(1-ve**2))

## calculate the critical buckling stress of stringer
alfa = 0.8
C1 = 0.425 # SSFS for 1 simply supported 1 free for 3 flanges
C2 = 4.0 # SSSS simply supported both sides for 1 flange
n = 0.6 # ????????

sigma_cr_str1 =  sigma_yield* alfa * ( (C1/sigma_yield) * np.pi**2 * E * (t_str / b_str1)**2 / (12*(1-ve**2)))**(1-n)
sigma_cr_str2 =  sigma_yield* alfa * ( (C2/sigma_yield) * np.pi**2 * E * (t_str / b_str2)**2 / (12*(1-ve**2)))**(1-n)
sigma_cr_str = (3*sigma_cr_str1*b_str1*t_str + sigma_cr_str2*b_str2*t_str) / A_str
print(sigma_cr_str1,sigma_cr_str2,sigma_cr_str)


## calculate euler column buckling for stiffeners
cen_y_str = (2*b_str1*t_str*0.5*t_str + b_str2*t_str*0.5*b_str2 + b_str1*t_str*b_str2) / A_str
Ixx_str = (1/12)*b_str1*t_str**3 + (b_str1*t_str)*(cen_y_str-b_str2)**2 +\
          (1/12)*t_str*b_str2**3 + (b_str2*t_str)*(cen_y_str-0.5*b_str2)**2 +\
       2*((1/12)*b_str1*t_str**3 + (b_str1*t_str)*(cen_y_str-0.5*t_str)**2)
rad_gyr = np.sqrt(Ixx_str / A_str )

effective_length_str = rad_gyr*np.sqrt(2*np.pi**2*E / sigma_cr_str)

#sigma_cr = sigma_cr_sk
sigma_cr = (sigma_cr_str*A_str + sigma_cr_sk*A_sk) / (A_str + A_sk)

## calculate mass
L_airframe = 6 
m_sk = density*1000 * boom_distance*len(x_boomcoor)*L_airframe*t_sk
m_str = density*1000 * len(x_boomcoor)*A_stringer*L_airframe
m_frame = density*1000 * boom_distance*len(x_boomcoor)*A_stringer*5
m = m_sk + m_str + m_frame
print(m)



plt.subplot(211)
plt.plot(poslist, momentlist)
plt.plot(poslist, shearlist)
# results
plt.subplot(212)
plt.axis([-1.5,1.5,-1.5,1.5])
#plt.plot(x_coorlist, y_coorlist)
plt.scatter(x_boomcoor, y_boomcoor)
plt.hlines(cen_y,-1,1)

output = pandas.DataFrame([maxsigma, sigma_cr, maxtau, maxvonmises,maxsigma/sigma_cr],
                          ['max sigma:','max critical buckling stress:','max tau', 'max vonmises','sigma ratio'])
                           
                           
print(output)



#def sc(Ixx, Iyy, Ixy, Mx, My, x, y):
#    return ((Ixx*My - Ixy*Mx)*x/(Ixx*Iyy-Ixy**2) + (Iyy*Mx - Ixy*My)*y/(Ixx*Iyy - Ixy**2))
#Ixx = Ixx
#Iyy = Iyy
#Ixy= -0
#x = (x_boomcoor[i]-0)
#y = (y_boomcoor[i]-cen_y)
#Mx = -maxmoment
#My = 0
#print(sc(Ixx, Iyy, Ixy, Mx, My, x, y))