# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:46:44 2019

Function to read DAT file from WAMIT output

@author: dprat
"""
import matplotlib.pyplot as plt
import numpy as np
import re
import os

files=os.listdir()
padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"

for file in files:
    if file.endswith('pan.dat'):
        dat_name = file
        print('## DAT file Found: ' + file)
              
dat_arq_handle = open(dat_name,'r')
arq_dat = dat_arq_handle.readlines()
dat_arq_handle.close()

cont=0
pos_zone=[]
zone_type=[]
for line in arq_dat:
    if 'ZONE' in line or 'zone' in line:
        pos_zone.append(cont)
        if 'F=FEPOINT' in line:
            zone_type.append('FEPOINT')
        elif 'F=POINT' in line:
            zone_type.append('POINT')
    cont += 1

I = []
J = []
for ii in pos_zone:
    # Read I, J
    pos_aux = re.search('I=',arq_dat[ii])
    pos_I = pos_aux.start()
    pos_aux = re.search('J=',arq_dat[ii])
    pos_J = pos_aux.start()
    pos_aux = re.search('F=',arq_dat[ii])
    pos_F = pos_aux.start()
    I_aux = re.findall(padrao,arq_dat[ii][pos_I:pos_J])
    J_aux = re.findall(padrao,arq_dat[ii][pos_J:pos_F])
    I.append(int(I_aux[0]))
    J.append(int(J_aux[0]))
    points=[]
    faces=[]
    if 'F=FEPOINT' in arq_dat[ii]:
        points_txt = arq_dat[ii+1:ii+1+I[-1]]
        faces_txt = arq_dat[ii+1+I[-1]:ii+2+I[-1]+J[-1]]
        for pts in points_txt:
            points.append([float(i) for i in re.findall(padrao,pts)])
        points = np.array(points)
        for fcs in faces_txt:
            faces.append([int(i) for i in re.findall(padrao,fcs)])
        faces = np.array(faces)
    if 'F=POINT' in arq_dat[ii]:
        points = arq_dat[ii+1:ii+1+I[-1]+J[-1]]
        
    
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
import scipy as sp

ax = a3.Axes3D(pl.figure())

x_max,y_max,z_max = np.max(points,axis=0)
x_min,y_min,z_min = np.min(points,axis=0)


for i in faces:
    vtx = [[points[i[0]-1,0],points[i[0]-1,1],points[i[0]-1,2]],
           [points[i[1]-1,0],points[i[1]-1,1],points[i[1]-1,2]],
           [points[i[2]-1,0],points[i[2]-1,1],points[i[2]-1,2]],
           [points[i[3]-1,0],points[i[3]-1,1],points[i[3]-1,2]]]
    tri = a3.art3d.Poly3DCollection([vtx])
#    tri.set_color(colors.rgb2hex([0,0,1])
    tri.set_edgecolor('k')
    ax.add_collection3d(tri)

max_range = np.array([x_max-x_min, y_max-y_min, z_max-z_min]).max() / 2.0

mid_x = (x_max+x_min) * 0.5
mid_y = (y_max+y_min) * 0.5
mid_z = (z_max+z_min) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

pl.show()