# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 16:56:31 2018

@author: dprat
"""
import re
import os
import numpy as np
from geomdl import BSpline
# from geomdl.visualization import VisMPL
# from matplotlib import cm
from numpy import linalg as la

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import rhino

class Properties:
    def __init__(self):
        pass

class GDF:
    def __init__(self,file_name):
        self.file_name = file_name
        properties = Properties
        self.properties = properties
        self.read_gdf()
        self.scaled = False
        # self.props()
  
    def read_gdf(self):
        fileID = open(self.file_name)
        arqGDF=fileID.readlines()
        fileID.close()   
        # standard for reading numbers
        padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"        
        # Line 1 -> Head
        # Line 2 -> ULEN e g
        [ulen, g] = [float(i) for i in re.findall(padrao, arqGDF[1])]
        #Line 3 -> Ix e Iy
        [Isx, Isy] = [int(i) for i in re.findall(padrao, arqGDF[2])]
        #Line 4 -> Npatch e IGDEF
        [Npatch, IGDEF] = [int(i) for i in re.findall(padrao, arqGDF[3])]
        # reading numbers
        x = []
        
        for j in range(4,len(arqGDF)):
            k = [float(i) for i in re.findall(padrao,arqGDF[j])]
            for i in k:
                x.append(i)
        
        # Variables initialization
        NUG = [None]*Npatch
        NVG = [None]*Npatch
        KUG = [None]*Npatch
        KVG = [None]*Npatch
        NUA = [None]*Npatch
        NVA = [None]*Npatch
        Mu = [None]*Npatch
        Mv = [None]*Npatch
        NB = [None]*Npatch
        VKNTUG = [None]*Npatch
        VKNTVG = [None]*Npatch
        XCOEF_1 = [None]*Npatch
        XCOEF_2 = [None]*Npatch
        XCOEF_3 = [None]*Npatch
        
        line = 0
        for ii in range(Npatch):
            NUG[ii] = int(x[line])
            line+=1
            NVG[ii] = int(x[line])
            line+=1
            KUG[ii] = int(x[line])
            line+=1
            KVG[ii] = int(x[line])
            line+=1
            
            NUA[ii] = NUG[ii] + 2*KUG[ii] - 1
            NVA[ii] = NVG[ii] + 2*KVG[ii] - 1
            Mu[ii] = NUG[ii] + KUG[ii] - 1
            Mv[ii] = NVG[ii] + KVG[ii] - 1
            NB[ii] = Mu[ii]*Mv[ii]
            
            VKNTUG[ii] = []
            for j in range(int(NUA[ii])):
                VKNTUG[ii].append(x[line])
                line+=1
            
            VKNTVG[ii] = []
            for j in range(int(NVA[ii])):
                VKNTVG[ii].append(x[line])
                line+=1
                
            XCOEF_1[ii] = []
            XCOEF_2[ii] = []
            XCOEF_3[ii] = []
            for j in range(int(NB[ii])):
                XCOEF_1[ii].append(x[line])
                line+=1
                XCOEF_2[ii].append(x[line])
                line+=1
                XCOEF_3[ii].append(x[line])
                line+=1
            
#        [L,B,D,K] = evalExtremes(Isx,Isy,XCOEF_1,XCOEF_2,XCOEF_3)
#        print([L,B,D,K])
        self.ulen = ulen
        self.g = g
        self.Isx = Isx
        self.Isy = Isy
        self.Npatch = Npatch
        self.IGDEF = IGDEF
        self.NUG = NUG
        self.NVG = NVG
        self.KUG = KUG
        self.KVG = KVG
        self.NUA = NUA
        self.NVA = NVA
        self.Mu = Mu
        self.Mv = Mv
        self.NB = NB
        self.VKNTUG = VKNTUG
        self.VKNTVG = VKNTVG
        self.XCOEF_1 = XCOEF_1
        self.XCOEF_2 = XCOEF_2
        self.XCOEF_3 = XCOEF_3
        
    def gdf2crtlpoints(self):
    
        ctrlpoints = [None]*self.Npatch
        
        for ii in range(self.Npatch):
            x = np.array(self.XCOEF_1[ii])
            x = x.reshape(self.Mv[ii],self.Mu[ii])
            x = x.transpose()
            x = x.reshape(self.Mu[ii]*self.Mv[ii],)
            y = np.array(self.XCOEF_2[ii])
            y = y.reshape(self.Mv[ii],self.Mu[ii])
            y = y.transpose()
            y = y.reshape(self.Mu[ii]*self.Mv[ii],)
            z = np.array(self.XCOEF_3[ii])
            z = z.reshape(self.Mv[ii],self.Mu[ii])
            z = z.transpose()
            z = z.reshape(self.Mu[ii]*self.Mv[ii],)
            
            ctrlpoints[ii] = []
            for jj in range(len(x)):
                ctrlpoints[ii].append([x[jj],y[jj],z[jj]])
            
        return ctrlpoints
    
    def props(self, T=[], KG=[]):
        print(' ')
        print('***************************')
        print('  Properties of ' + self.file_name) 
        print(' ')
        # Variables Initialization     
        faces = [None]*self.Npatch
        centroid = [None]*self.Npatch
        vertices = [None]*self.Npatch
        normal = [None]*self.Npatch
        area = [None]*self.Npatch
        area_total = [None]*self.Npatch
        surf = [None]*self.Npatch
        pontos = self.gdf2crtlpoints()
        Volx = [None]*self.Npatch
        Voly = [None]*self.Npatch
        Volz = [None]*self.Npatch
        Vx = [None]*self.Npatch
        Vy = [None]*self.Npatch
        Vz = [None]*self.Npatch
        LCF = [None]*self.Npatch
        BM = [None]*self.Npatch
        BMl = [None]*self.Npatch
        Awl = [None]*self.Npatch
        
        self.irr=[]
        [L,B,D,KGe] = self.evalExtremes()

        if T == []:
            T = D

        rho = 1.025

        # Evaluate vertices and faces
        new_faces = []
        for ii in range(self.Npatch):
            # determines which patches have irr=1
            if abs(np.max(self.XCOEF_3[ii])-np.min(self.XCOEF_3[ii]))<0.2:
                self.irr.append(1)
            else:
                self.irr.append(0)
            
            surf[ii] = BSpline.Surface()
            surf[ii].degree_u = self.KUG[ii] - 1
            surf[ii].degree_v = self.KVG[ii] - 1        
            surf[ii].set_ctrlpts(pontos[ii],int(self.Mu[ii]),int(self.Mv[ii]))
            surf[ii].knotvector_u = self.VKNTUG[ii]
            surf[ii].knotvector_v = self.VKNTVG[ii] 
    
            # Set evaluation delta
            #surf.delta = 0.025
            surf[ii].delta_u = .02 #1/(2*self.Mu[ii]) #0.025
            surf[ii].delta_v = .02 #1/(2*self.Mv[ii]) #0.025#
            
            surf[ii].evaluate()
            surf[ii].tessellate() # The tesselation function provides triangles faces in which the evaluations will be done
            # if ii ==1:
            #     # Set visualization component
            #     surf[ii].vis = VisMPL.VisSurface(ctrlpts=False, legend=False)

            #     # Plot the surface
            #     surf[ii].render(colormap=cm.terrain)
            
            vertices[ii] = np.array([list(surf[ii].tessellator.vertices[i].data) for i in range(len(surf[ii].tessellator.vertices))])
            faces[ii] = np.array([list(surf[ii].tessellator.faces[i].data) for i in range(len(surf[ii].tessellator.faces))])

            vertices[ii][:,2] = vertices[ii][:,2] + np.abs(KGe) # Move keel to z=0

            z_max = np.max(vertices[ii][:,2])

            z_min = np.min(vertices[ii][:,2])

            z_test = vertices[ii][:,2]
            
            if self.irr[ii] == 0:
                ## Find faces with all vertices above waterline
                id_above = z_test[faces[ii]] > T
                id_all_above = np.sum(id_above,1) == 3
                id_2_above = np.sum(id_above,1) == 2
                id_1_above = np.sum(id_above,1) == 1
                id_all_under = np.sum(id_above,1) == 0
                # Cut panels above water line

                # create new panels
                #Panels with all vtx under water line
                vtx_0 = vertices[ii][faces[ii][id_all_under],:]

                #Panels with 1 vtx above water line
                vtx_1 = vertices[ii][faces[ii][id_1_above],:]
                new_faces1 = []
                for vtx_aux in vtx_1:
                    test_aux = vtx_aux[:,2] > T
                    x_aux = vtx_aux[:,0]
                    y_aux = vtx_aux[:,1]
                    z_aux = vtx_aux[:,2]
                    # Case 1
                    if np.array_equal(test_aux, np.array([True, False, False])):
                        
                        p1_x = np.interp(T, [z_aux[0], z_aux[1]], [x_aux[0], x_aux[1]])
                        p1_y = np.interp(T, [z_aux[0], z_aux[1]], [y_aux[0], y_aux[1]])

                        p2_x = np.interp(T, [z_aux[0], z_aux[2]], [x_aux[0], x_aux[2]])
                        p2_y = np.interp(T, [z_aux[0], z_aux[2]], [y_aux[0], y_aux[2]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([p1, vtx_aux[1,:], vtx_aux[2,:]])
                        fc2 = np.array([p1, vtx_aux[1,:], p2])
                        
                        new_faces1.append(fc1)
                        new_faces1.append(fc2)
                    # Case 2
                    elif np.array_equal(test_aux,np.array([False, True, False])):

                        p1_x = np.interp(T, [z_aux[0], z_aux[1]], [x_aux[0], x_aux[1]])
                        p1_y = np.interp(T, [z_aux[0], z_aux[1]], [y_aux[0], y_aux[1]])

                        p2_x = np.interp(T, [z_aux[1], z_aux[2]], [x_aux[1], x_aux[2]])
                        p2_y = np.interp(T, [z_aux[1], z_aux[2]], [y_aux[1], y_aux[2]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([vtx_aux[0,:], p1, vtx_aux[2,:]])
                        fc2 = np.array([p1, p2, vtx_aux[2,:]])

                        new_faces1.append(fc1)
                        new_faces1.append(fc2)
                    # Case 3
                    else:

                        p1_x = np.interp(T, [z_aux[0], z_aux[2]], [x_aux[0], x_aux[2]])
                        p1_y = np.interp(T, [z_aux[0], z_aux[2]], [y_aux[0], y_aux[2]])

                        p2_x = np.interp(T, [z_aux[1], z_aux[2]], [x_aux[1], x_aux[2]])
                        p2_y = np.interp(T, [z_aux[1], z_aux[2]], [y_aux[1], y_aux[2]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([vtx_aux[0,:], vtx_aux[1,:], p2])
                        fc2 = np.array([vtx_aux[0,:], p2, p1])
                        
                        new_faces1.append(fc1)
                        new_faces1.append(fc2)

                # ############################
                # # debug (must be commented)    
                # fig = plt.figure()
                # ax = fig.add_subplot(111, projection='3d') 
                # nfx = np.reshape(np.array(new_faces1),(len(new_faces1)*3,3))
                # nfx1 = np.reshape(np.array(vtx_1),(len(vtx_1)*3,3))
                # ax.scatter(nfx[:,0], nfx[:,1], nfx[:,2], zdir='z', s=20)
                # ax.scatter(nfx1[:,0], nfx1[:,1], nfx1[:,2], zdir='z', s=20)
                # for nfxx in new_faces1:
                #     ct = [np.mean(nfxx[:,0]), np.mean(nfxx[:,1]), np.mean(nfxx[:,2])]
                #     ux = [nfxx[1,0] - nfxx[0,0], nfxx[1,1] - nfxx[0,1], nfxx[1,2] - nfxx[0,2]]
                #     vx = [nfxx[2,0] - nfxx[0,0], nfxx[2,1] - nfxx[0,1], nfxx[2,2] - nfxx[0,2]]
                #     nx = np.cross(ux,vx)/2
                #     ax.plot([ct[0], ct[0]+nx[0]], [ct[1], ct[1]+nx[1]], [ct[2],ct[2]+nx[2]], color='black')
                #     ax.plot(nfxx[[0,1,2,0],0], nfxx[[0,1,2,0],1], nfxx[[0,1,2,0],2], color='blue')

                # for nfxx in vtx_1:
                #     ct = [np.mean(nfxx[:,0]), np.mean(nfxx[:,1]), np.mean(nfxx[:,2])]
                #     ux = [nfxx[1,0] - nfxx[0,0], nfxx[1,1] - nfxx[0,1], nfxx[1,2] - nfxx[0,2]]
                #     vx = [nfxx[2,0] - nfxx[0,0], nfxx[2,1] - nfxx[0,1], nfxx[2,2] - nfxx[0,2]]
                #     nx = np.cross(ux,vx)/2
                #     ax.plot([ct[0], ct[0]+nx[0]], [ct[1], ct[1]+nx[1]], [ct[2],ct[2]+nx[2]],color='red')
                #     ax.plot(nfxx[[0,1,2,0],0], nfxx[[0,1,2,0],1], nfxx[[0,1,2,0],2], color='green')

                # # plt.show()
                # ############################

                #Panels with 2 vtx above water line
                vtx_2 = vertices[ii][faces[ii][id_2_above],:]
                new_faces2 = []
                for vtx_aux in vtx_2:
                    test_aux = vtx_aux[:,2] > T
                    x_aux = vtx_aux[:,0]
                    y_aux = vtx_aux[:,1]
                    z_aux = vtx_aux[:,2]
                    # Case 1
                    if np.array_equal(test_aux, np.array([True, True, False])):
                        
                        p1_x = np.interp(T, [z_aux[0], z_aux[2]], [x_aux[0], x_aux[2]])
                        p1_y = np.interp(T, [z_aux[0], z_aux[2]], [y_aux[0], y_aux[2]])

                        p2_x = np.interp(T, [z_aux[1], z_aux[2]], [x_aux[1], x_aux[2]])
                        p2_y = np.interp(T, [z_aux[1], z_aux[2]], [y_aux[1], y_aux[2]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([p1, p2, vtx_aux[2,:]])
                        
                        new_faces2.append(fc1)
                    # Case 2
                    elif np.array_equal(test_aux,np.array([False, True, True])):

                        p1_x = np.interp(T, [z_aux[0], z_aux[1]], [x_aux[0], x_aux[1]])
                        p1_y = np.interp(T, [z_aux[0], z_aux[1]], [y_aux[0], y_aux[1]])

                        p2_x = np.interp(T, [z_aux[0], z_aux[2]], [x_aux[0], x_aux[2]])
                        p2_y = np.interp(T, [z_aux[0], z_aux[2]], [y_aux[0], y_aux[2]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([vtx_aux[0,:], p1, p2])

                        new_faces2.append(fc1)
                    # Case 3
                    else:

                        p1_x = np.interp(T, [z_aux[1], z_aux[2]], [x_aux[1], x_aux[2]])
                        p1_y = np.interp(T, [z_aux[1], z_aux[2]], [y_aux[1], y_aux[2]])

                        p2_x = np.interp(T, [z_aux[1], z_aux[0]], [x_aux[1], x_aux[0]])
                        p2_y = np.interp(T, [z_aux[1], z_aux[0]], [y_aux[1], y_aux[0]])

                        p1 = np.array([p1_x, p1_y, T])
                        p2 = np.array([p2_x, p2_y, T])

                        fc1 = np.array([p2, vtx_aux[1,:], p1])
                        
                        new_faces2.append(fc1)

                # ############################
                # # debug (must be commented)    
                # fig = plt.figure()
                # ax2 = fig.add_subplot(111, projection='3d') 
                # nfx = np.reshape(np.array(new_faces2),(len(new_faces2)*3,3))
                # nfx1 = np.reshape(np.array(vtx_1),(len(vtx_1)*3,3))
                # ax2.scatter(nfx[:,0], nfx[:,1], nfx[:,2], zdir='z', s=20)
                # ax2.scatter(nfx1[:,0], nfx1[:,1], nfx1[:,2], zdir='z', s=20)
                # for nfxx in new_faces2:
                #     ct = [np.mean(nfxx[:,0]), np.mean(nfxx[:,1]), np.mean(nfxx[:,2])]
                #     ux = [nfxx[1,0] - nfxx[0,0], nfxx[1,1] - nfxx[0,1], nfxx[1,2] - nfxx[0,2]]
                #     vx = [nfxx[2,0] - nfxx[0,0], nfxx[2,1] - nfxx[0,1], nfxx[2,2] - nfxx[0,2]]
                #     nx = np.cross(ux,vx)/2
                #     ax2.plot([ct[0], ct[0]+nx[0]], [ct[1], ct[1]+nx[1]], [ct[2],ct[2]+nx[2]], color='black')
                #     ax2.plot(nfxx[[0,1,2,0],0], nfxx[[0,1,2,0],1], nfxx[[0,1,2,0],2], color='blue')

                # for nfxx in vtx_1:
                #     ct = [np.mean(nfxx[:,0]), np.mean(nfxx[:,1]), np.mean(nfxx[:,2])]
                #     ux = [nfxx[1,0] - nfxx[0,0], nfxx[1,1] - nfxx[0,1], nfxx[1,2] - nfxx[0,2]]
                #     vx = [nfxx[2,0] - nfxx[0,0], nfxx[2,1] - nfxx[0,1], nfxx[2,2] - nfxx[0,2]]
                #     nx = np.cross(ux,vx)/2
                #     ax2.plot([ct[0], ct[0]+nx[0]], [ct[1], ct[1]+nx[1]], [ct[2],ct[2]+nx[2]],color='red')
                #     ax2.plot(nfxx[[0,1,2,0],0], nfxx[[0,1,2,0],1], nfxx[[0,1,2,0],2], color='green')

                # plt.show()
                # ############################

                if new_faces1 == []:
                    new_faces1 = np.array(new_faces1).reshape((0,3,3))

                if new_faces2 == []:
                    new_faces2 = np.array(new_faces2).reshape((0,3,3))

                new_faces_aux = np.append(vtx_0, new_faces1, axis=0)
                new_faces_aux = np.append(new_faces_aux, new_faces2, axis=0)
                new_faces.append(new_faces_aux)
            else:
                new_faces.append([])

        # ############################
        # # debug (must be commented)    
        # fig = plt.figure()
        # ax2 = fig.add_subplot(111, projection='3d') 
        # nfx = np.reshape(np.array(new_faces),(len(new_faces)*3,3))
        
        # ax2.scatter(nfx[:,0], nfx[:,1], nfx[:,2], zdir='z', s=20)
        
        # for nfxx in new_faces:
        #     ct = [np.mean(nfxx[:,0]), np.mean(nfxx[:,1]), np.mean(nfxx[:,2])]
        #     ux = [nfxx[1,0] - nfxx[0,0], nfxx[1,1] - nfxx[0,1], nfxx[1,2] - nfxx[0,2]]
        #     vx = [nfxx[2,0] - nfxx[0,0], nfxx[2,1] - nfxx[0,1], nfxx[2,2] - nfxx[0,2]]
        #     nx = np.cross(ux,vx)/2
        #     ax2.plot([ct[0], ct[0]+nx[0]], [ct[1], ct[1]+nx[1]], [ct[2],ct[2]+nx[2]], color='black')
        #     ax2.plot(nfxx[[0,1,2,0],0], nfxx[[0,1,2,0],1], nfxx[[0,1,2,0],2], color='blue')

        # plt.show()
        # ############################

        for ii in range(len(new_faces)):
            
            normal[ii] = []
            area[ii] = []
            centroid[ii] = []

            for fc in new_faces[ii]:
                # Centroide Evaluation
                centroid[ii].append([np.mean(fc[:,0]), np.mean(fc[:,1]), np.mean(fc[:,2]) - T])
                
                # Normal and area evaluation
                u = [fc[1,0] - fc[0,0], fc[1,1] - fc[0,1], fc[1,2] - fc[0,2]]
                v = [fc[2,0] - fc[0,0], fc[2,1] - fc[0,1], fc[2,2] - fc[0,2]]
                normal_aux = - np.cross(u,v)/2 # the tesselation generates normals point to inside vessel, so it is necessary the negative signal
                area[ii].append(la.norm(normal_aux))
                normal[ii].append(normal_aux)
            
            area[ii] = np.array(area[ii])
            area_total[ii] = np.sum(area[ii])
            centroid[ii] = np.array(centroid[ii])
            normal[ii] = np.array(normal[ii])
            prod = 2*np.sum([self.Isx, self.Isy])
            if normal[ii].size == 0:
                Volx[ii] = 0.0
                Voly[ii] = 0.0
                Volz[ii] = 0.0
                Vx[ii] = 0.0
                Vy[ii] = 0.0
                Vz[ii] = 0.0
                LCF[ii] = 0.0
                BM[ii]  = 0.0
                BMl[ii] = 0.0
                Awl[ii] = 0.0
            else:
                Volx[ii] = np.sum(normal[ii][:,0] * centroid[ii][:,0]) * prod
                Voly[ii] = np.sum(normal[ii][:,1] * centroid[ii][:,1]) * prod
                Volz[ii] = np.sum(normal[ii][:,2] * centroid[ii][:,2]) * prod
                Vx[ii] = np.sum(normal[ii][:,0] * centroid[ii][:,0]**2) / np.sum(normal[ii][:,0] * centroid[ii][:,0]) / 2
                Vy[ii] = np.sum(normal[ii][:,1] * centroid[ii][:,1]**2) / np.sum(normal[ii][:,0] * centroid[ii][:,0]) / 2
                Vz[ii] = np.sum(normal[ii][:,2] * centroid[ii][:,2]**2) / np.sum(normal[ii][:,0] * centroid[ii][:,0]) / 2
                LCF[ii] = np.sum(normal[ii][:,2] * centroid[ii][:,0]) / np.sum(normal[ii][:,2])
                BM[ii]  = -np.sum(normal[ii][:,2] * centroid[ii][:,1]**2) / np.sum(normal[ii][:,2] * centroid[ii][:,2])
                BMl[ii] = -np.sum(normal[ii][:,2] * centroid[ii][:,0]**2) / np.sum(normal[ii][:,2] * centroid[ii][:,2])
                Awl[ii] = -np.sum(normal[ii][:,2])
            
        flag_irr =  abs(np.array(self.irr)-1)   
        
        Volx = np.sum(flag_irr * Volx)
        Voly = np.sum(flag_irr * Voly)
        Volz = np.sum(flag_irr * Volz)
        
        if self.Isx == 1:
            Vx = np.float64(0.0)
        else:
            Vx = np.sum(flag_irr * Vx)
            
        if self.Isy == 1:
            Vy = np.float64(0.0)
        else:
            Vy = np.sum(flag_irr * Vy)

        if ((self.Isx == 1) and (self.Isy == 1)):
            awl_coef = 4
        elif ((self.Isx == 1) or (self.Isy == 1)):
            awl_coef = 2
        else:
            awl_coef = 1
        
        Vz = np.sum(flag_irr * Vz)
        
        CB = np.array([Vx, Vy, Vz])
        LCF = np.sum(flag_irr * LCF)
        BM = np.sum(flag_irr * BM)
        BMl = np.sum(flag_irr * BMl)
        Awl = np.sum(flag_irr * Awl) * awl_coef
        Volume = np.mean([Volx, Voly, Volz])
        KB = T + Vz

        if KG == []:
            KG = KGe

        GM  = KB + BM  - KG
        GMl = KB + BMl - KG
        
        ZPOT = KG - T

        # Mass matrix
        Rxx = 0.35 * B
        Ryy = 0.25 * L
        coef = 0.05 # five percent of critical damping

        Mass = np.zeros((6,6))
        Mass[0,0] = Volume * rho
        Mass[1,1] = Volume * rho
        Mass[2,2] = Volume * rho
        Mass[3,3] = Volume * rho * (Rxx ** 2)
        Mass[4,4] = Volume * rho * (Ryy ** 2)
        Mass[5,5] = Volume * rho * (Rxx ** 2) + Volume * rho * (Ryy ** 2)

        wn3 = np.sqrt(rho * self.g * Awl / (Mass[2,2] * 1.3))
        wn4 = np.sqrt(GM * self.g * Mass[1,1] / (Mass[3,3] * 1.5))
        wn5 = np.sqrt(GMl * self.g * Mass[1,1] / (Mass[4,4] * 2))

        Be = np.zeros((6,6))
        Be[2,2] = coef * 2 * wn3 * Mass[2,2] * 1.3
        Be[3,3] = coef * 2 * wn4 * Mass[3,3] * 1.5
        Be[4,4] = coef * 2 * wn5 * Mass[4,4] * 2

        # prints
        print('water density = ' + "{:.3f}".format(rho) + ' t/m³')
        print('L = ' + "{:.2f}".format(L) + ' m')
        print('B = ' + "{:.2f}".format(B) + ' m')
        print('D = ' + "{:.2f}".format(D) + ' m')
        print('T = ' + "{:.2f}".format(T) + ' m')
        print('KG = ' + "{:.2f}".format(KG) + ' m')
        print('CB = [' + "{:.2f}".format(CB[0]) + ', ' + "{:.2f}".format(CB[1]) + ', ' + "{:.2f}".format(CB[2]) + '] m')
        print('KB = ' + "{:.2f}".format(KB) + ' m')
        print('BM = ' + "{:.2f}".format(BM) + ' m')
        print('BMl = ' + "{:.2f}".format(BMl) + ' m')
        print('GM = ' + "{:.2f}".format(GM) + ' m')
        print('GMl = ' + "{:.2f}".format(GMl) + ' m')
        print('ZPOT = ' + "{:f}".format(ZPOT) + ' m')
        print('Awl = ' + "{:.2f}".format(Awl) + ' m²')
        print('Vol = ' + "{:.2f}".format(Volume) + ' [' + "{:.2f}".format(Volx) + ', ' + "{:.2f}".format(Voly) + ', ' + "{:.2f}".format(Volz) + '] m³') 
        print(' ')
        print('  End of gdf_class.props()') 
        print('***************************')
        print(' ')
        
        # Exporting properties
        self.properties.L = L
        self.properties.B = B
        self.properties.D = D
        self.properties.T = T
        self.properties.KG = KG
        self.properties.rho = rho
        self.properties.area = area
        self.properties.area_total = area_total
        self.properties.normal = normal
        self.properties.centroid = centroid
        self.properties.Volx = Volx
        self.properties.Voly = Voly
        self.properties.Volz = Volz
        self.properties.Volume = Volume
        self.properties.Vx = Vx
        self.properties.Vy = Vy
        self.properties.Vz = Vz
        self.properties.CB = CB
        self.properties.LCF = LCF
        self.properties.KB = KB
        self.properties.Awl = Awl
        self.properties.BM = BM
        self.properties.BMl = BMl
        self.properties.GM = GM
        self.properties.GMl = GMl
        self.properties.ZPOT = ZPOT
        self.properties.Mass = Mass
        self.properties.Be = Be
        pass

    def gdf_props_write(self, prop_name, gdf_out=[]):
        if gdf_out==[]:
            gdf_out = self.file_name

        prop = open(prop_name,'w')
        # prints
        prop.write('Properties of model: ' + gdf_out + '\n\n')
        prop.write('water density = ' + "{:.3f}".format(self.properties.rho) + ' t/m³' + '\n')
        prop.write('L = ' + "{:.2f}".format(self.properties.L) + ' m\n')
        prop.write('B = ' + "{:.2f}".format(self.properties.B) + ' m\n')
        prop.write('D = ' + "{:.2f}".format(self.properties.D) + ' m\n')
        prop.write('T = ' + "{:.2f}".format(self.properties.T) + ' m\n')
        prop.write('KG = ' + "{:.2f}".format(self.properties.KG) + ' m\n')
        prop.write('CB = [' + "{:.2f}".format(self.properties.CB[0]) + ', ' + "{:.2f}".format(self.properties.CB[1]) + ', ' + "{:.2f}".format(self.properties.CB[2]) + ']' + ' m\n')
        prop.write('KB = ' + "{:.2f}".format(self.properties.KB) + ' m\n')
        prop.write('BM = ' + "{:.2f}".format(self.properties.BM) + ' m\n')
        prop.write('BMl = ' + "{:.2f}".format(self.properties.BMl) + ' m\n')
        prop.write('GM = ' + "{:.2f}".format(self.properties.GM) + ' m\n')
        prop.write('GMl = ' + "{:.2f}".format(self.properties.GMl) + ' m\n')
        prop.write('ZPOT = ' + "{:f}".format(self.properties.ZPOT) + ' m\n')
        prop.write('Awl = ' + "{:.2f}".format(self.properties.Awl) + ' m²\n')
        prop.write('Vol = ' + "{:.2f}".format(self.properties.Volume) + ' [' + "{:.2f}".format(self.properties.Volx) + ', ' + "{:.2f}".format(self.properties.Voly) + ', ' + "{:.2f}".format(self.properties.Volz) + ']' + ' m³\n') 
        prop.write('\nFORCE FILE:\n')
        prop.write('1\n')
        for m1 in self.properties.Mass:
            for m2 in m1:
                prop.write("{:.5e}".format(m2))
                prop.write(" ")
            prop.write('\n')
        prop.write('1\n')
        for m1 in self.properties.Be:
            for m2 in m1:
                prop.write("{:.5e}".format(m2))
                prop.write(" ")
            prop.write('\n')
        prop.write('1\n')
        for m1 in np.zeros((6,6)):
            for m2 in m1:
                prop.write("{:.5e}".format(m2))
                prop.write(" ")
            prop.write('\n')    

        pass
    
    def gdf_write(self, gdf_out, path_out=[], LWL_rhino=False):
        if path_out != []:
            # os.mkdir(path_out,exist_ok=True)
            os.makedirs(path_out, exist_ok=True)
            gdf_name = path_out + '\\' + gdf_out
        else:
            gdf_name = gdf_out
        
        gdf = open(gdf_name,'w')

        gdf.write('SCALED GDF FILE\n')
        gdf.write(str(self.ulen) + ' ' + str(self.g) + ' ULEN GRAV\n')
        gdf.write(str(self.Isx) + ' ' + str(self.Isy) + ' ISX  ISY\n')
        gdf.write(str(np.sum(self.irr)) + ' ' + str(self.IGDEF) + ' NPATCH IGDEF\n')
        for ii in range(self.Npatch):
            if self.irr[ii] == 0:
                gdf.write(str(self.NUG[ii]) + ' ' +  str(self.NVG[ii]) + '\n')
                gdf.write(str(self.KUG[ii]) + ' ' +  str(self.KVG[ii]) + '\n')
                for lt in self.VKNTUG[ii]:
                    gdf.write(' ' + str(lt) + '\n')

                for lt in self.VKNTVG[ii]:
                    gdf.write(' ' + str(lt) + '\n')
                
                for lt in range(len(self.XCOEF_1[ii])):
                    gdf.write(' ' + str(self.XCOEF_1[ii][lt] + self.properties.CB[0]) +' ' + str(self.XCOEF_2[ii][lt]) + ' ' + str(self.XCOEF_3[ii][lt]) + '\n')
                gdf.write('\n')
        gdf.close()

        self.gdf_props_write(gdf_name[:-4] + '_props.txt', gdf_out)

        if LWL_rhino == True:
            rhino.LWL_ho(gdf_out, self.properties.ZPOT, True, path_gdf=path_out)
        pass
    
    def gdf_scale(self, Lf, Bf, Df, T=[], KG=[]):
        # Final Dimensions
        #    Lf = 250
        #    Bf = 45
        #    Df = 15
        
        [L,B,D,K] = self.evalExtremes()
        
        if KG == []:
            KG = K

        if hasattr(self,'DimOriginal') == False:
            self.DimOriginal = [L,B,D,K]
        
        # Exceptions for flat surfaces
        if L==0:
            L=1
            
        if B==0:
            B=1
            
        if D==0:
            D=1
            
        #Scale Factors
        f_L = Lf/L
        f_B = Bf/B
        f_D = Df/D
        
        XCOEF_1_n = [np.array(self.XCOEF_1[jj])*f_L for jj in range(self.Npatch)]
        XCOEF_2_n = [np.array(self.XCOEF_2[jj])*f_B for jj in range(self.Npatch)]
        XCOEF_3_n = [(np.array(self.XCOEF_3[jj])+K)*f_D-KG for jj in range(self.Npatch)]
                
        self.XCOEF_1 = XCOEF_1_n
        self.XCOEF_2 = XCOEF_2_n
        self.XCOEF_3 = XCOEF_3_n
        
        [L1,B1,D1,K1] = self.evalExtremes()
        
        self.props(T)
        pass
    
    def evalExtremes(self):
        x_max = max([max(self.XCOEF_1[i]) for i in range(len(self.XCOEF_1))])
        x_min = min([min(self.XCOEF_1[i]) for i in range(len(self.XCOEF_1))])
        
        y_max = max([max(self.XCOEF_2[i]) for i in range(len(self.XCOEF_2))])
        y_min = min([min(self.XCOEF_2[i]) for i in range(len(self.XCOEF_2))])
        
        z_max = max([max(self.XCOEF_3[i]) for i in range(len(self.XCOEF_3))])
        z_min = min([min(self.XCOEF_3[i]) for i in range(len(self.XCOEF_3))])
        
        L = abs(x_max-x_min)
        B = abs(y_max-y_min)
        D = abs(z_max-z_min)
        
        if self.Isx==1:
            L*=2
        
        if self.Isy==1:
            B*=2
            
        KGe = abs(z_min)
        
#        print('Extremos: L = ' + str(L) + '\n          B = ' + str(B) + '\n          D = ' + str(D) + '\n          KGe = ' + str(KGe))
        
        return [L,B,D,KGe]

        



# # Example
# # loading gdf file
# painel = GDF('ship2.gdf')

# # scale gdf file
# painel.gdf_scale(300, 40, 20, T=10, KG=14)

# # exporting gdf
# painel.gdf_write('ship_scaled.gdf')