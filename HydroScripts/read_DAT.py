# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 17:46:44 2019

Function to read DAT file from WAMIT output

@author: dprat
"""
import numpy as np
import re
import os

def low_order():
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
    points_list=[]
    faces_list=[]
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
            points_txt = arq_dat[ii+1:ii+1+I[-1]+J[-1]]
            for pts in points_txt:
                points.append([float(i) for i in re.findall(padrao,pts)])
        points = np.array(points)
        points_list.append(points)    
        faces_list.append(faces)
    return [points_list,faces_list]        