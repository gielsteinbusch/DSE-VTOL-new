# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 12:55:21 2019

@author: archi
"""
import numpy as np 
import matplotlib.pyplot as plt 
from airfoil import list_x, list_z
from sparwars import Blade_loading
from formula_beam import *
from mpl_toolkits.mplot3d import Axes3D

#input values                       (MUST MATCH THE SPARWARS VALUES)
radius = 4.46
taper = 1
chord_length = 1
inc_angle = 5
twist = 10
V_flight = 50
rpm = 41 * 60/(2*np.pi)
rho = 0.5
CL = 1.2
W_aircraft = 2500
LDratio = 15
disc_steps = 10
nseg = disc_steps -1 
skin_thickness = list(0.04*np.ones(nseg))
#material properities           NEED TO BE CHANGED 
den = 1.58e3
E = 110e9
Uten = 356e8
G = 23e9

#call the functions from the other file -----------
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


#set up a class for simple beam theory 
class Simple_Beam:
    def sections(self):
            self.nseg = disc_steps -1 
            self.radius = radius
            self.chord_length = chord_length
            self.taper = taper
            self.disc_steps = disc_steps
            self.V_flight = V_flight
            #making the midpoint of sections into list
            self.length_ds = np.linspace(0.6, radius, self.disc_steps)
            self.sec_points = []
            self.w_segment = self.length_ds[1] - self.length_ds[0]
            for d in self.length_ds[:-1]: 
                midpoint = d + self.w_segment/2
                self.sec_points.append(midpoint)
            #making the midpoint chords a list  
            self.taperchord = np.linspace(chord_length, taper*chord_length , self.disc_steps)
            self.change_chord = (self.taperchord[0] - self.taperchord[1])/2
            self.sec_chord = []
            for c in self.taperchord[:-1]:
                c = c - self.change_chord
                self.sec_chord.append(c)
            #making the midpoint twist a list
            self.twisting = np.deg2rad(-1* np.linspace(inc_angle, inc_angle- twist , self.disc_steps))
            self.sec_twist = []
            for t in self.twisting[:-1]: 
                change_twist = (self.twisting[0]-self.twisting[1])/2
                t = t - change_twist
                self.sec_twist.append(t)
            self.x_coordinates = np.array(list_x)
            self.z_coordinates = np.array(list_z)
            self.LDratio = LDratio
            self.list_z = list_z
            self.list_x = list_x
    #new lift distribution based on the lift at the midpoint of each segment        
    def lift_distribution(self, V_flight):
        self.x_list = []
        self.lift_list = []
        self.totallift = 0
        x = []
        self.Vpoint_list = []
        for segment in range(self.nseg):
            distance_center = self.sec_points[segment] + (1/6)*self.w_segment  
            # we take the velocity at the 2/3 point due to the nature of the lift distribution, however the rest is calculated at average for segment so halfway!!
            V_rotor = 2*np.pi*(rpm/60)*distance_center
            V_total = V_flight + V_rotor
            S = (np.pi * (self.length_ds[segment +1])**2) - sum(x)
            #print(self.length_ds[segment +1])
            x.append(S)
            L = 0.5 * rho * (V_total)**2 * S * CL
            self.lift_list.append(L)
            self.x_list.append(distance_center)
            self.totallift += L
            self.Vpoint_list.append(distance_center)
    #new shear distribution based on the previous lift 
    def shear_distribution(self):
        self.totalmoment = 0 
        self.shear_list = []
        for s in range(self.nseg): 
            shear_mid = ((1/3)*blade.shearforce_list[s] + (2/3)*blade.shearforce_list[s+1])/2
            self.shear_list.append(shear_mid)
            self.totalmoment += self.lift_list[s]*self.x_list[s]
    #new moment on previous shear 
    def moment_distribution(self):
        self.moment_list = []
        for i in range(self.nseg):
           self.moment = (blade.moment_list[i+1] + blade.moment_list[i])/2 
           self.moment_list.append(self.moment)
    #new profiles at the midpoints 
    def profile(self):
        self.profile_x = []
        self.profile_z = []
        self.profile_y = []
        for step in range(self.nseg):
            c = self.sec_chord[step]
            self.profile_x.append(self.x_coordinates*c)  
            self.profile_z.append(self.z_coordinates*c) 
            for i in range(len(self.x_coordinates)):
                self.profile_y.append(self.sec_points[step])
    # adding spar coordinates for new cross-sections 
    def spar_coor(self, skin_thickness):
        self.xspar_list = []
        self.zspar_list = []
        self.len_spar = []
        self.area_spar = []
        self.cen_spar_z = []
        self.loc_spar = []
        for step in range(self.nseg):
            self.t_spar = 2*skin_thickness[step]
            loc_spar_cs = 0.3*self.sec_chord[step]            
            x_spar = list((np.ones(int((len(list_x)- 1)/4)))*(loc_spar_cs))
            z_spar = list(np.linspace(max(self.profile_z[step]), min(self.profile_z[step]),int((len(self.profile_x[step])-1)/4) ))
            len_spar_cs = max(z_spar)-min(z_spar)
            area_spar_cs = len_spar_cs*self.t_spar
            cen_spar_z_cs = len_spar_cs/2
            self.xspar_list.append(x_spar)
            self.zspar_list.append(z_spar)
            self.loc_spar.append(loc_spar_cs)
            self.len_spar.append(len_spar_cs)
            self.area_spar.append(area_spar_cs)
            self.cen_spar_z.append(cen_spar_z_cs)
        self.plot_x = []
        self.plot_y = []
        self.plot_z = []
        for step in range(self.nseg):
            x_all = list(self.profile_x[step][:-1] )+ list(self.xspar_list[step])
            self.plot_x.append(x_all)
            y_all = list(self.Vpoint_list[step] * np.ones(len(self.profile_x[step]) -1 + len(self.xspar_list[step])))
            self.plot_y.append(y_all)
            z_all = list(self.profile_z[step][:-1]) + list(self.zspar_list[step])
            self.plot_z.append(z_all)
#PROFILES OF THE SHEAR STUFF FOR THESE NEW CROSS SECTIONS 
    def profile_new(self):
        self.profile1_x = []
        self.profile2_x = []
        self.profile1_z = []
        self.profile2_z = []
        self.segment1_list = []
        self.segment2_list = []
        for step in range(self.nseg):
            a = float(self.xspar_list[step][0])
            profile1_x_cs = []
            profile2_x_cs = []
            profile1_z_cs = []
            profile2_z_cs = []
            self.count_top1 = 0
            for c in range(len(self.profile_x[step])):
                if self.profile_x[step][c] > a:
                    profile1_x_cs.append(self.profile_x[step][c])
                    profile1_z_cs.append(self.profile_z[step][c])
                else: 
                    if profile1_x_cs[-1] != a:
                        profile2_x_cs.append(profile1_x_cs[-1])
                        profile2_z_cs.append(profile1_z_cs[-1])
                        self.idxtop2 = len(profile1_x_cs) -1
                        for i in range(len(self.xspar_list[step])):
                            profile1_x_cs.append(self.xspar_list[step][i])
                            profile1_z_cs.append(self.zspar_list[step][i])
                        self.idxbottom2 = len(profile1_x_cs) 
                        #print(self.idxtop2, self.idxbottom2)
                    profile2_x_cs.append(self.profile_x[step][c])
                    profile2_z_cs.append(self.profile_z[step][c])
            for i in range(1,len(self.xspar_list[step])+1):
                profile2_x_cs.append(self.xspar_list[step][-i])
                profile2_z_cs.append(self.zspar_list[step][-i])
            self.profile1_x.append(profile1_x_cs)
            self.profile2_x.append(profile2_x_cs)
            self.profile1_z.append(profile1_z_cs)
            self.profile2_z.append(profile2_z_cs)
            segment1_list_cs = []
            segment2_list_cs = []
            for i in range(len(profile1_x_cs)-1):
                segment_length1 = np.sqrt((profile1_x_cs[i+1]-profile1_x_cs[i])**2 + (profile1_z_cs[i+1]-profile1_z_cs[i])**2)
                segment1_list_cs.append(segment_length1)
            for i in range(len(profile2_x_cs)-1):
                segment_length2 = np.sqrt((profile2_x_cs[i+1]-profile2_x_cs[i])**2 + (profile2_z_cs[i+1]-profile2_z_cs[i])**2)
                segment2_list_cs.append(segment_length2)
            self.segment1_list.append(segment1_list_cs)
            self.segment2_list.append(segment2_list_cs)
#twist for the cross sectional points 
    def twist(self):
        self.profile3_x = []
        self.profile4_x = []
        self.profile3_z = []
        self.profile4_z = []
        for step in range(self.nseg):
            self.profile3_x_cs = []
            self.profile4_x_cs = []
            self.profile3_z_cs = []
            self.profile4_z_cs = []
            twiz = self.sec_twist[step]
            for i in range(len(self.profile1_x[step])):
                self.profile3_x_cs.append(self.profile1_x[step][i]*np.cos(twiz) - self.profile1_z[step][i]*np.sin(twiz)) 
                self.profile3_z_cs.append(self.profile1_z[step][i]*np.cos(twiz) + self.profile1_x[step][i]*np.sin(twiz))
            for i in range(len(self.profile2_x[step])):
                self.profile4_x_cs.append(self.profile2_x[step][i]*np.cos(twiz) - self.profile2_z[step][i]*np.sin(twiz)) 
                self.profile4_z_cs.append(self.profile2_z[step][i]*np.cos(twiz) + self.profile2_x[step][i]*np.sin(twiz))
            self.profile3_x.append(self.profile3_x_cs)
            self.profile3_z.append(self.profile3_z_cs)
            self.profile4_x.append(self.profile4_x_cs)
            self.profile4_z.append(self.profile4_z_cs)
            #new centroid calculations for the cross sectionss !!! ------------------------------------
    def center_gravity(self, skin_thickness):
        self.centroids = []
        self.cen_x_list = []
        self.cen_z_list = []
        self.segment_list = []
        for step in range(self.nseg):
            cen_x_list_cs = []
            cen_z_list_cs = []
            segment_list_cs = []
            for i in range(len(self.x_coordinates)-1):
                segment_length = np.sqrt((self.profile_x[step][i+1]-self.profile_x[step][i])**2 + (self.profile_z[step][i+1]-self.profile_z[step][i])**2)
                cen_x = self.profile_x[step][i]+(self.profile_x[step][i+1]-self.profile_x[step][i])/2
                cen_z = self.profile_z[step][i]+(self.profile_z[step][i+1]-self.profile_z[step][i])/2
                segment_list_cs.append(segment_length)
                cen_x_list_cs.append(cen_x) 
                cen_z_list_cs.append(cen_z)
            centroid_x = (np.sum(skin_thickness[step]*np.array(segment_list_cs)*np.array(cen_x_list_cs)) + (self.area_spar[step]*self.loc_spar[step]))   / (np.sum(skin_thickness[step]*np.array(segment_list_cs)) + self.area_spar[step])  
            centroid_z = (np.sum(skin_thickness[step]*np.array(segment_list_cs)*np.array(cen_z_list_cs)) + (self.area_spar[step]*self.cen_spar_z[step])) / (np.sum(skin_thickness[step]*np.array(segment_list_cs)) + self.area_spar[step])
            self.centroids.append([centroid_x,centroid_z])  
            self.cen_x_list.append(cen_x_list_cs)
            self.cen_z_list.append(cen_z_list_cs)
            self.segment_list.append(segment_list_cs)
    #NOW REDOING MOMENTS OF INERTIA FOR THE CROSS SECTIONSSS 
    def inertia(self, skin_thickness):
        self.ix_list = []
        self.iz_list = []
        self.ixz_list = []
        for step in range(self.nseg):
            ix = 0
            iz = 0
            ixz = 0
            twiz = self.twisting[step]
            for i in range(len(self.x_coordinates)-1):
                iz += (self.segment_list[step][i]*skin_thickness[step])*(self.cen_x_list[step][i]-self.centroids[step][0])**2
                ix += (self.segment_list[step][i]*skin_thickness[step])*(self.cen_z_list[step][i]-self.centroids[step][1])**2
                ixz += (self.segment_list[step][i]*skin_thickness[step])*(self.cen_x_list[step][i]-self.centroids[step][0])*(self.cen_z_list[step][i]-self.centroids[step][1])
            ix_spar = (1/12) * self.t_spar * (self.len_spar[step])**3 + self.area_spar[step] * (self.cen_spar_z[step] - self.centroids[step][1])**2
            iz_spar = (1/12) * (self.t_spar)**3 * self.len_spar[step] + self.area_spar[step] * (self.loc_spar[step] - self.centroids[step][0])**2
            ixz_spar = self.area_spar[step] * (self.cen_spar_z[step] - self.centroids[step][1]) * (self.loc_spar[step] - self.centroids[step][0])
            if twiz != 0:
                ix_spar = (ix_spar + iz_spar)/2 + ((ix_spar - iz_spar)/2)*np.cos(2*twiz) - ixz_spar*np.sin(2*twiz)
                iz_spar = (ix_spar + iz_spar)/2 - ((ix_spar - iz_spar)/2)*np.cos(2*twiz) + ixz_spar*np.sin(2*twiz)
                ixz_spar = ((ix_spar - iz_spar)/2)*np.sin(2*twiz) + ixz_spar*np.cos(2*twiz) 
            self.ix_list.append(ix + ix_spar)
            self.iz_list.append(iz + iz_spar)
            self.ixz_list.append(ixz + ixz_spar)
    #doing the internal areas for each section!! 
    def area(self, skin_thickness):
        self.area_list = []
        for step in range(self.nseg):
            A3 = 0
            for i in range(len(self.profile3_x[step])-1):
                ax,az = self.profile3_x[step][i], self.profile3_z[step][i]
                bx,bz = self.profile3_x[step][i+1], self.profile3_z[step][i+1]
                u = np.array([ax-self.xspar_list[step][0],az-self.zspar_list[step][0]])
                v = np.array([bx-self.xspar_list[step][0],bz-self.zspar_list[step][0]])
                parallelogram = np.cross(u,v)
                dA = 0.5*parallelogram
                A3 += dA 
            A4 = 0
            for i in range(len(self.profile4_x[step])-1):
                ax,az = self.profile4_x[step][i], self.profile4_z[step][i]
                bx,bz = self.profile4_x[step][i+1], self.profile4_z[step][i+1]
                u = np.array([ax,az])
                v = np.array([bx,bz])
                parallelogram = np.cross(u,v)
                dA = 0.5*parallelogram
                A4 += dA 
            self.area_list.append([A3+A4,A3,A4])
    #centrifugal force
    def centrifugal_force(self, skin_thickness):
        self.CSAlist = []
        self.siglist = []
        self.mlist = []
        self.ylist= []
        mass = 0 
        self.centrifugal = []
        self.wrs = (rpm*2*np.pi)/60
        self.w_segment = self.length_ds[1]-self.length_ds[0]
        for i in range(self.nseg):
            Aspar = self.area_spar[i]
            Askin = sum(self.segment_list[i])*skin_thickness[i]
            CSA = Aspar + Askin
            self.CSAlist.append(CSA)
            m = CSA*den*self.w_segment
            mass +=m
            #print(m)
            ypoint = self.Vpoint_list[i]
            # we also take the centrifugal force at that pesky 2/3 point as well as it is quadratically distributed
            self.ylist.append(ypoint)
            centri = den*CSA*((ypoint*self.w_segment)/2)*self.wrs**2   
            self.centrifugal.append(centri)
            sigi = den*((ypoint*self.w_segment)/2)*self.wrs**2        
            #sigi = Ni/CSA
            self.siglist.append(sigi)
    
    def deflections(self):
        self.niter = 15
        iteration = 0
        self.deflist = [np.zeros(self.nseg)] 
        while iteration < self.niter:
            deflections = []
            self.moments = []
#            change1 = list(self.deflist[iteration])[0]
#            changedefl = [change1]
#            for i in range(len(self.deflist[0])-1): 
#                change = list(self.deflist[iteration])[i+1] - list(self.deflist[iteration])[i]
#                changedefl.append(change)                   
            #print(changedefl)
            for step in range(self.nseg):
                    if iteration == 0:                
                        self.Ma = self.totalmoment
                    else:
                        self.Ma = sum(np.array(self.lift_list)*(np.array(self.Vpoint_list )-blade.length_ds[0])) - Mp
                        #print(Ma)
                    Ra = self.totallift
                    ix = self.ix_list[step]
                    z = self.Vpoint_list[step] - blade.length_ds[0]
                    W = self.lift_list
                    P = self.centrifugal
                    Pd = np.sum(P) * self.deflist[iteration][step]
                    #print(np.sum(P), self.deflist[iteration][step], Pd)
                    Mp = sum(np.array(P)*self.deflist[iteration])
                    #Pdelta = np.array(P)*np.array(changedefl)
                    
                    Pdelta = 0
                    for x in range(self.nseg):
                        if list(self.deflist[iteration])[-1] > 0:
                            Pdelta += P[x]*(list(self.deflist[iteration])[step]-list(self.deflist[iteration])[x])*np.heaviside(list(self.deflist[iteration])[step]-list(self.deflist[iteration])[x],0.5)
                        else: 
                            Pdelta += P[x]*(list(self.deflist[iteration])[step]-list(self.deflist[iteration])[x])*np.heaviside(list(self.deflist[iteration])[x]-list(self.deflist[iteration])[step],0.5)
                    la = z*np.ones(self.nseg) - (np.array(self.Vpoint_list )-blade.length_ds[0])
                    #print(len(la))
                    la = la * np.heaviside(la, 0.5)
                    #print(la)
                    Wla = np.sum(np.array(W)*la)
                    #print(Wla)
                    Wla3 = np.sum(np.array(W)*(la**3))
                    #print(Wla)
                    Pla2 = np.sum(Mp*(la**2))
                    Pla3 = Pdelta*la**2
                    Pdla2 = Pd*la**2
                    #print(Mp)
                    Moment = self.Ma + Wla - Ra*z - Pdelta + Pd #sum(Mp)
                    #print(Pdelta, P[step], Moment, Pla2)
                    self.moments.append(Moment)
                    #print(Moment)
                    vEI = (self.Ma*(z**2)/2) + Wla3/6 - Ra*(z**3)/6 - Pla3/2 + Pdla2/2
                    #print(Pla2)
                    v = np.average(vEI/(E*ix))
                    deflections.append(v)
            #print(Pdelta)
            self.deflist.append(np.array(deflections))
            iteration +=1 
            
                    #print(Moment)
        #print(self.M_it)
    
            
# BEAM THEORY -----------------------------------------------------------------------------------------------
#    def actual_loads(self): 
#        self.Vlist = []
#        self.Mlist = []
#        self.deflections = []
#        MA = 0 
#        RA = 0
#        for step in range(self.nseg): 
#            l = self.w_segment
#            a = (1/3)*l
#            x = l #(this can be changed depending on where you want to determine moment in each discretised section)
#            ix = self.ix_list[step]
#            P = sum(self.centrifugal[step:])
#            W = sum(self.lift_list[step:])
#            print(ix, P,W, P/W)
#            k = K(E, ix, P)
#            c1 = C1(k,l)
#            c2 = C2(k,l)
#            ca3 = Ca3(k,l,a)
#            ca4 = Ca4(k,l,a)
#            theta_A = theta_A1(W, P, k, l, a)
#            yA = yA1(W, P, k, l ,a)
#            f1 = F1(k,x)
#            f2 = F2(k,x)
#            f3 = F3(k,x)
#            f4 = F4(k,x)
#            fa1 = Fa1(k,x,a)
#            fa2 = Fa2(k,x,a)
#            fa3 = Fa3(k,x,a)
#            fa4 = Fa4(k,x,a)
#            LTv = LTV(W, fa1)
#            LTm = LTM(W,k,fa2)
#            LTth = LTtheta(W,P,fa3)
#            LTy = LTY(W,P,k,fa4)
#            V = RA*f1 + MA*k*f2 + theta_A*P*f1 + LTv
#            M = MA*f1 + (RA/k)*f2 + (theta_A*P/k)*f2 + LTm
#            theta = theta_A*f1 + (MA*k/P)*f2 + (RA/P)*f3 + LTth
#            y = yA + (theta_A/k)*f2 + (MA/P)*f3 + (RA/(P*k))*f4 + LTy #(should just be zero at beams end)
#            print(yA,y)##########################################################- now store shit and sum to complete integration 
#            #RA = V
#            #MA = M
#            self.Vlist.append(V)
#            self.Mlist.append(M)
#            self.deflections.append(yA)
        #print(sum(self.deflections))
        #print(self.Mlist)
    def bending_stress(self):
        self.stress_x_list = []
        self.stress_z_list = []
        self.sigma_list = []
        for step in range(self.nseg):
            stress_x_list_cs = []
            stress_z_list_cs = []
            sigma_list_cs = []
            moment_x = self.moments[step]             ### moment due to centrifugal and lift from the function above
            moment_z = self.moment_list[step]/LDratio ### can just leave as pure lift/drag moment
            ix = self.ix_list[step]
            iz = self.iz_list[step]
            ixz = self.ixz_list[step]
            for i in range(len(self.x_coordinates)-1):
                stress_x = ((iz*moment_x - ixz*moment_z) / (ix*iz - ixz**2)) * (self.cen_z_list[step][i] - self.centroids[step][1] ) 
                stress_z = ((ix*moment_z - ixz*moment_x) / (ix*iz - ixz**2)) * (self.cen_x_list[step][i] - self.centroids[step][0] ) 
                sigma = stress_x + stress_z
                stress_x_list_cs.append(stress_x)
                stress_z_list_cs.append(stress_z)
                sigma_list_cs.append(sigma)
            for i in range(len(self.xspar_list[step])):
                stress_x_spar =  ((iz*moment_x - ixz*moment_z) / (ix*iz - ixz**2)) * (self.zspar_list[step][i] - self.centroids[step][1] )
                stress_z_spar = ((ix*moment_z - ixz*moment_x) / (ix*iz - ixz**2)) * (self.xspar_list[step][i] - self.centroids[step][0] ) 
                sigma_spar = stress_x_spar + stress_z_spar
                sigma_list_cs.append(sigma_spar)
            self.stress_x_list.append(stress_x_list_cs)  
            self.stress_z_list.append(stress_z_list_cs)
            self.sigma_list.append(sigma_list_cs)
   
    def shear_stress(self, skin_thickness):
        self.tau_list = []
        self.qb_list3 = []
        self.qb_list4 = []
        self.total_shear_list = []
        self.intM3_list = []
        self.red_shear3 = []
        self.intM4_list = []
        self.red_shear4 = []
        self.sec_flows = []
        self.qflow3 = []
        self.qflow4 = []
        self.shear_flow = []
        for step in range(self.nseg):
            twiz = self.twisting[step]
            ix = self.ix_list[step]
            iz = self.iz_list[step]
            ixz = self.ixz_list[step]
            Sz = self.shear_list[step]
            Sx = Sz / LDratio
            if twiz <  0: 
                Sz = Sz*np.cos(twiz)  - Sx*np.sin(twiz)
                Sx = Sx*np.cos(twiz) + Sz*np.sin(twiz)
            if twiz > 0 : 
                Sz = Sz*np.cos(twiz) - Sx*np.sin(twiz)
                Sx = Sx*np.cos(twiz) + Sz*np.sin(twiz)
            qx_coef = -(Sx*ix - Sz*ixz) / (ix*iz - ixz**2)
            qz_coef = -(Sz*iz - Sx*ixz) / (ix*iz - ixz**2)
            qb3 = 0
            qb3_list_cs = []
            for i in range(self.idxtop2 +1 ):
                qbx0 = qx_coef * skin_thickness[step] * (self.profile3_x[step][i+1] - self.centroids[step][0]) * self.segment1_list[step][i]
                qbz0 = qz_coef * skin_thickness[step] * (self.profile3_z[step][i+1] - self.centroids[step][1]) * self.segment1_list[step][i]
                qb3 += (qbx0 + qbz0)
                qb3_list_cs.append(qb3)
            for i in range(self.idxtop2, self.idxtop2 +1 +len(self.xspar_list[0]) ):
                qbx0 = qx_coef * self.t_spar * (self.profile3_x[step][i+1] - self.centroids[step][0]) * self.segment1_list[step][i]
                qbz0 = qz_coef * self.t_spar * (self.profile3_z[step][i+1] - self.centroids[step][1]) * self.segment1_list[step][i]
                qb3 += (qbx0 + qbz0)
                qb3_list_cs.append(qb3)
            for i in range(self.idxbottom2 , len(self.profile3_x[0])-1):
                qbx0 = qx_coef * skin_thickness[step] * (self.profile3_x[step][i+1] - self.centroids[step][0]) * self.segment1_list[step][i]
                qbz0 = qz_coef * skin_thickness[step] * (self.profile3_z[step][i+1] - self.centroids[step][1]) * self.segment1_list[step][i]
                qb3 += (qbx0 + qbz0)
                qb3_list_cs.append(qb3)
            self.qb_list3.append(qb3_list_cs)
            qb4 = 0
            qb4_list_cs = []
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]) + 1):
                qbx0 = qx_coef * skin_thickness[step] * (self.profile4_x[step][i+1] - self.centroids[step][0]) * self.segment2_list[step][i]
                qbz0 = qz_coef * skin_thickness[step] * (self.profile4_z[step][i+1] - self.centroids[step][1]) * self.segment2_list[step][i]
                qb4 += (qbx0 + qbz0)
                qb4_list_cs.append(qb4)
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]), len(self.profile4_x[0]) -1):
                qbx0 = qx_coef * self.t_spar * (self.profile4_x[step][i+1] - self.centroids[step][0]) * self.segment2_list[step][i]
                qbz0 = qz_coef * self.t_spar * (self.profile4_z[step][i+1] - self.centroids[step][1]) * self.segment2_list[step][i]
                qb4 += (qbx0 + qbz0)
                qb4_list_cs.append(qb4)
            self.qb_list4.append(qb4_list_cs)
            Mint3 = 0 
            for i in range( len(self.profile3_x[0])-1):
                 qbx3 = ((self.profile3_x[step][i+1] - self.profile3_x[step][i]) / self.segment1_list[step][i]) *self.qb_list3[step][i]
                 qbz3 = ((self.profile3_z[step][i+1] - self.profile3_z[step][i]) / self.segment1_list[step][i]) *self.qb_list3[step][i]
                 m3x = qbx3 * self.segment1_list[step][i] *  (self.profile3_x[step][i] - self.profile3_x[step][self.idxbottom2])
                 m3z = qbz3 * self.segment1_list[step][i] * (self.profile3_z[step][i] - self.profile3_z[step][self.idxbottom2])
                 Mint3 += (m3x + m3z)
            self.intM3_list.append(Mint3)
            self.red_shear3.append(Mint3/(-2*self.area_list[step][1]))
            Mint4 = 0 
            for i in range(len(self.profile4_x[0]) -1 ):
                qbx4 = ((self.profile4_x[step][i+1] - self.profile4_x[step][i]) / self.segment2_list[step][i]) *self.qb_list4[step][i]
                qbz4 = ((self.profile4_z[step][i+1] - self.profile4_z[step][i]) / self.segment2_list[step][i]) *self.qb_list4[step][i]
                m4x = qbx4 * self.segment2_list[step][i] *  (self.profile4_x[step][i] - self.profile4_x[step][len(self.profile4_x[0]) - len(self.xspar_list[0])])
                m4z = qbz4 * self.segment2_list[step][i] * (self.profile4_z[step][i] - self.profile4_z[step][len(self.profile4_x[0]) - len(self.xspar_list[0])])
                Mint4 += (m4x + m4z)
            self.intM4_list.append(Mint4)
            self.red_shear4.append(Mint4/(-2*self.area_list[step][2]))
            self.qR3 = []
            # now you have found all that you need for the first equation, only now need to do the bredt-batho stuff..
            ##do for section 3 
            co3 = 1/(2 * self.area_list[step][1] * G)
            co4 = 1/(2 * self.area_list[step][2] * G)
            Bb3 = 0
            Bbq3 = 0
            for i in range(self.idxtop2 +1 ):
                Bb3 += co3 * self.segment1_list[step][i]/skin_thickness[step]
                Bbq3 += co3 * (self.segment1_list[step][i]*self.qb_list3[step][i])/skin_thickness[step]
            for i in range(self.idxtop2, self.idxtop2 +1 +len(self.xspar_list[0]) ):
                Bb3 += co3 * self.segment1_list[step][i]/self.t_spar
                Bbq3 += co3 * (self.segment1_list[step][i]*self.qb_list3[step][i])/self.t_spar
            for i in range(self.idxbottom2 , len(self.profile3_x[0])-1):
                Bb3 += co3 * self.segment1_list[step][i]/skin_thickness[step]
                Bbq3 += co3 * (self.segment1_list[step][i]*self.qb_list3[step][i])/skin_thickness[step]
            Bb4 = 0 
            Bbq4 = 0 
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]) + 1):
                Bb4 += co4 * self.segment2_list[step][i]/skin_thickness[step]
                Bbq4 += co4 * (self.segment2_list[step][i]*self.qb_list4[step][i])/skin_thickness[step]
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]), len(self.profile4_x[0]) -1):
                Bb4 += co4 * self.segment2_list[step][i]/self.t_spar
                Bbq4 += co4 * (self.segment2_list[step][i]*self.qb_list4[step][i])/self.t_spar
            A = np.mat([[2*self.area_list[step][1], 2*self.area_list[step][2]], 
                        [-Bb3 , Bb4]])    
            B = np.mat([[- (Mint3 + Mint4)],
                         [ (Bbq3 - Bbq4) ]])
            qvec = np.linalg.solve(A,B)
            self.sec_flows.append(qvec)
            #print(Bbq3, Bbq4)
            #need to now add the redundant shear of each section to the base shear of each section 
            Q3 = []
            for i in range(len(self.profile3_x[0])):
                Q3.append( self.qb_list3[step][i] + float(self.sec_flows[step][0]))
            Q4 = []
            for i in range(len(self.profile4_x[0])): 
                Q4.append( self.qb_list4[step][i] + float(self.sec_flows[step][1]))
            self.qflow3.append(Q3)
            self.qflow4.append(Q4)
            
            # Now need to put all the shear flows into one list in the right order and with the two spar sides summed together properly 
            qpart1 = []
            taupart1 = []
            for i in range(self.idxtop2 ): 
                qpart1.append(Q3[i])
                taupart1.append(Q3[i]/skin_thickness[step]) 
            qpart2 = []
            taupart2 = []
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]) ):
                qpart2.append(Q4[i])
                taupart2.append(Q4[i]/skin_thickness[step])
            qpart3 = [] 
            taupart3 = []
            for i in range(self.idxbottom2 , len(self.profile3_x[0])-1):
                qpart3.append(Q3[i])
                taupart3.append(Q3[i]/skin_thickness[step])
            qspar_sec4 = []
            qspar_sec3 = []
            for i in range(self.idxtop2, self.idxtop2 +1 +len(self.xspar_list[0]) -1 ):
                qspar_sec3.append(Q3[i])
            for i in range(len(self.profile4_x[0]) - len(self.xspar_list[0]) -1, len(self.profile4_x[0]) -1):
                qspar_sec4.append(Q4[i])
            #print(qspar_sec3,qspar_sec4)
            #now need to combine the two spars 
            qpart4 = []
            taupart4 = []
            for i in range(len(qspar_sec4)): 
                qpart4.append(qspar_sec3[i] + qspar_sec4[i])
                taupart4.append((qspar_sec3[i] + qspar_sec4[i])/self.t_spar)
            self.flow_list = [] + qpart1 + qpart2 + qpart3 + qpart4
            self.shear_flow.append(self.flow_list)
            self.tau = [] + taupart1 + taupart2 + taupart3 + taupart4 
            self.tau_list.append(self.tau)
        #print(Q3, qpart3, qvec)
        
        
        
        
    def axial(self):
        self.ax_stress = []
        P = []
        #print(self.centrifugal)
        for step in range(self.nseg):
            p = sum(self.centrifugal[step:])
            P.append(p)
        #print(P)
        for step in range(self.nseg):
            sig_norm = (P[step])/(self.CSAlist[step])
            self.ax_stress.append(sig_norm)
        #print(self.ax_stress)
        
    def sigma_total(self): 
        self.stresses = []
        for step in range(self.nseg):
            ax = self.ax_stress[step]
            sigmas = []
            for i in range(len(self.sigma_list[0])): 
                sigma = self.sigma_list[step][i] + ax
                sigmas.append(sigma)
            self.stresses.append(sigmas)
    
    def von_mises(self):
        self.vonmises = []
        for step in range(self.nseg):
            von_mises_cs = []
            for i in range(len(self.sigma_list[0])):
                von_mises_segment = np.sqrt((self.stresses[step][i])**2 + 3*(self.tau_list[step][i])**2)
                von_mises_cs.append(von_mises_segment)
            self.vonmises.append(von_mises_cs)

    ## we can add and alter the stress calculations after the beam theory stuff is finalised ########################################
            
    
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

#for i in range(beam.niter):
    #plt.plot(beam.Vpoint_list, beam.deflist[i])
#
#for i in range(beam.nseg):
#    plt.plot(beam.profile3_x[i], beam.profile3_z[i])
#    plt.plot(beam.profile4_x[i], beam.profile4_z[i]) 
# 
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
#    
