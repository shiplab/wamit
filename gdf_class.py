# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 16:56:31 2018

@author: dprat
"""
import re
import numpy as np
from geomdl import BSpline
#from geomdl.visualization import VisMPL as vis
from numpy import linalg as la

class Propriedades:
    def __init__(self):
        pass

class GDF:
    def __init__(self,nome_arquivo):
        self.nome_arquivo = nome_arquivo
        propriedades = Propriedades
        self.propriedades = propriedades
        self.le_gdf()
        self.props()
  
    def le_gdf(self):
        fileID = open(self.nome_arquivo)
        arqGDF=fileID.readlines()
        fileID.close()   
        #padrao para ler números
        padrao = "[+-]?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?"        
        #Linha 1 -> Cabeçalho
        #Linha 2 -> ULEN e g
        [ulen,g] = [float(i) for i in re.findall(padrao,arqGDF[1])]
        #Linha 3 -> Ix e Iy
        [Isx,Isy] = [int(i) for i in re.findall(padrao,arqGDF[2])]
        #Linha 4 -> Npatch e IGDEF
        [Npatch,IGDEF] = [int(i) for i in re.findall(padrao,arqGDF[3])]
        # Leitura dos numeros
        x = []
        
        for j in range(4,len(arqGDF)):
            k = [float(i) for i in re.findall(padrao,arqGDF[j])]
            for i in k:
                x.append(i)
        
        #Inicialização de variáveis
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
            
            NUA[ii] = NUG[ii] + 2*KUG[ii] - 1;
            NVA[ii] = NVG[ii] + 2*KVG[ii] - 1;
            Mu[ii] = NUG[ii] + KUG[ii] - 1;
            Mv[ii] = NVG[ii] + KVG[ii] - 1;
            NB[ii] = Mu[ii]*Mv[ii];
            
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
            
#        [L,B,D,K] = avaliaExtremos(Isx,Isy,XCOEF_1,XCOEF_2,XCOEF_3)
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
    
    def props(self):       
        faces = [None]*self.Npatch
        centroide = [None]*self.Npatch
        vertices = [None]*self.Npatch
        normais = [None]*self.Npatch
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
        
        self.irr=[]
        [L,B,D,KG] = self.avaliaExtremos()
        print('L = ' + "{:.2f}".format(L))
        print('B = ' + "{:.2f}".format(B))
        print('D = ' + "{:.2f}".format(D))
        print('KG = ' + "{:.2f}".format(KG))
        
        #definição do calado manual (deve ser implementado como chamada de função)
        
        T = 10
        
        for ii in range(self.Npatch):
            # determina quais quais patchs são irr=1
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
            surf[ii].delta_u = 1/(2*self.Mu[ii]) #0.025
            surf[ii].delta_v = 1/(2*self.Mv[ii]) #0.025#
            
            surf[ii].evaluate()
            surf[ii].tessellate()
            
            vertices[ii] = [list(surf[ii].tessellator.vertices[i].data) for i in range(len(surf[ii].tessellator.vertices))]
            faces[ii] = [list(surf[ii].tessellator.triangles[i].vertex_ids_zero) for i in range(len(surf[ii].tessellator.triangles))]
            
            ##### IMPLEMENTAR FUNÇÃO QUE REDEFINE VERTICES E FACES COM BASE NO CALADO
            # move quilha pra Z=0
            
            ###########################################################################
            
            normais[ii] = []
            area[ii] = []
            centroide[ii] = []

            for fc in faces[ii]:
                # Calculo do centroide
                xc = np.mean([vertices[ii][fc[0]][0], vertices[ii][fc[1]][0], vertices[ii][fc[2]][0]])
                yc = np.mean([vertices[ii][fc[0]][1], vertices[ii][fc[1]][1], vertices[ii][fc[2]][1]])
                zc = np.mean([vertices[ii][fc[0]][2], vertices[ii][fc[1]][2], vertices[ii][fc[2]][2]])+(KG-D)              
                centroide[ii].append([xc, yc, zc])
                
                # Calculo da normal e da area
                u = np.array([vertices[ii][fc[1]][0] - vertices[ii][fc[0]][0], vertices[ii][fc[1]][1] - vertices[ii][fc[0]][1], vertices[ii][fc[1]][2] - vertices[ii][fc[0]][2]])
                v = np.array([vertices[ii][fc[2]][0] - vertices[ii][fc[0]][0], vertices[ii][fc[2]][1] - vertices[ii][fc[0]][1], vertices[ii][fc[2]][2] - vertices[ii][fc[0]][2]])
                
                normal_aux = np.cross(u,v)/2
                area[ii].append(la.norm(normal_aux))
                normais[ii].append(normal_aux)
#                if la.norm(normal_aux) == 0:
#                    normais[ii].append(np.array([0, 0, 0]))
#                else:
#                    normais[ii].append(normal_aux/la.norm(normal_aux))
            
            area[ii] = np.array(area[ii])
            area_total[ii] = np.sum(area[ii])
            centroide[ii] = np.array(centroide[ii])
            normais[ii] = np.array(normais[ii])
            prod = 2*np.sum([self.Isx, self.Isy])
            Volx[ii] = np.sum(normais[ii][:,0]*centroide[ii][:,0])*prod
            Voly[ii] = np.sum(normais[ii][:,1]*centroide[ii][:,1])*prod
            Volz[ii] = np.sum(normais[ii][:,2]*centroide[ii][:,2])*prod
            Vx[ii] = np.sum(normais[ii][:,0]*centroide[ii][:,0]**2)/np.sum(normais[ii][:,0]*centroide[ii][:,0])/2
            Vy[ii] = np.sum(normais[ii][:,1]*centroide[ii][:,1]**2)/np.sum(normais[ii][:,0]*centroide[ii][:,0])/2
            Vz[ii] = np.sum(normais[ii][:,2]*centroide[ii][:,2]**2)/np.sum(normais[ii][:,0]*centroide[ii][:,0])/2
            LCF[ii] = np.sum(normais[ii][:,2]*centroide[ii][:,0])/np.sum(normais[ii][:,2]);
            BM[ii]  = -np.sum(normais[ii][:,2]*centroide[ii][:,1]**2)/np.sum(normais[ii][:,2]*centroide[ii][:,2]);
            BMl[ii] = -np.sum(normais[ii][:,2]*centroide[ii][:,0]**2)/np.sum(normais[ii][:,2]*centroide[ii][:,2]);
            
        flag_irr =  abs(np.array(self.irr)-1)   
        
        Volx = np.sum(flag_irr*Volx)
        Voly = np.sum(flag_irr*Voly)
        Volz = np.sum(flag_irr*Volz)
        
        if self.Isx == 1:
            Vx = np.float64(0.0)
        else:
            Vx = np.sum(flag_irr*Vx)
            
        if self.Isy == 1:
            Vy = np.float64(0.0)
        else:
            Vy = np.sum(flag_irr*Vy)
        
        Vz = np.sum(flag_irr*Vz) - (KG-D)
        
        CB = np.array([Vx,Vy,Vz])
        print('CB = [' + "{:.2f}".format(CB[0]) + ',' + "{:.2f}".format(CB[1]) + ',' + "{:.2f}".format(CB[2]) + ']')
        
        LCF = np.sum(flag_irr*LCF)
        BM = np.sum(flag_irr*BM)
        print('BM = ' + "{:.2f}".format(BM))
        BMl = np.sum(flag_irr*BMl)
        print('BMl = ' + "{:.2f}".format(BMl))
        
        Volume = np.mean([Volx,Voly,Volz])
        print('Vol = ' + "{:.2f}".format(Volume) + ' [' + "{:.2f}".format(Volx) + ',' + "{:.2f}".format(Voly) + ',' + "{:.2f}".format(Volz) + ']') 
        
        KB = KG + CB[2]
        print('KB = ' + "{:.2f}".format(KB))

        GM  = KB + BM  - KG        
        print('GM = ' + "{:.2f}".format(GM))
        
        GMl = KB + BMl - KG
        print('GMl = ' + "{:.2f}".format(GMl))
        
        self.propriedades.vertices = vertices
        self.propriedades.faces = faces
        self.propriedades.area = area
        self.propriedades.area_total = area_total
        self.propriedades.normais = normais
        self.propriedades.centroide = centroide
        self.propriedades.Volx = Volx
        self.propriedades.Voly = Voly
        self.propriedades.Volz = Volz
        self.propriedades.Volume = Volume
        self.propriedades.Vx = Vx
        self.propriedades.Vy = Vy
        self.propriedades.Vz = Vz
        self.propriedades.CB = CB
        self.propriedades.LCF = LCF
        self.propriedades.BM = BM
        self.propriedades.BMl = BMl
        self.propriedades.GM = GM
        self.propriedades.GMl = GMl
        pass
    
    def escreveGDF(self,gdf_saida):
#        [ulen,g,Isx,Isy,Npatch,IGDEF,NUG,NVG,KUG,KVG,NUA,NVA,Mu,Mv,NB,VKNTUG,VKNTVG,XCOEF_1,XCOEF_2,XCOEF_3] = infoGDF
        gdf = open(gdf_saida,'w')
        gdf.write('Arquivo gdf escalado\n')
        gdf.write(str(self.ulen) + ' ' + str(self.g) + ' ULEN GRAV\n')
        gdf.write(str(self.Isx) + ' ' + str(self.Isy) + ' ISX  ISY\n')
        gdf.write(str(self.Npatch) + ' ' + str(self.IGDEF) + ' NPATCH IGDEF\n')
        for ii in range(self.Npatch):
            gdf.write(str(self.NUG[ii]) + ' ' +  str(self.NVG[ii]) + '\n')
            gdf.write(str(self.KUG[ii]) + ' ' +  str(self.KVG[ii]) + '\n')
            for lt in self.VKNTUG[ii]:
                gdf.write(' ' + str(lt) + '\n')

            for lt in self.VKNTVG[ii]:
                gdf.write(' ' + str(lt) + '\n')
            
            for lt in range(len(self.XCOEF_1[ii])):
                gdf.write(' ' + str(self.XCOEF_1[ii][lt]) +' ' + str(self.XCOEF_2[ii][lt]) + ' ' + str(self.XCOEF_3[ii][lt]) + '\n')
            gdf.write('\n')
        gdf.close()
        pass
    
    def escalagdf(self,Lf,Bf,Df):
        #    [ulen,g,Isx,Isy,Npatch,IGDEF,NUG,NVG,KUG,KVG,NUA,NVA,Mu,Mv,NB,VKNTUG,VKNTVG,XCOEF_1,XCOEF_2,XCOEF_3] = infoGDF
        #    gdf_saida = 'scaled.gdf'
        # Dimensoes finais
        #    Lf = 250
        #    Bf = 45
        #    Df = 15
        
        [L,B,D,K] = self.avaliaExtremos()
        
        if hasattr(self,'DimOriginal') == False:
            self.DimOriginal = [L,B,D,K]
        
        # Exceções para superfícies planas    
        if L==0:
            L=1
            
        if B==0:
            B=1
            
        if D==0:
            D=1
            
        #Fatores de escala    
        f_L = Lf/L
        f_B = Bf/B
        f_D = Df/D
        
        XCOEF_1_n = [np.array(self.XCOEF_1[jj])*f_L for jj in range(self.Npatch)]
        XCOEF_2_n = [np.array(self.XCOEF_2[jj])*f_B for jj in range(self.Npatch)]
        XCOEF_3_n = [(np.array(self.XCOEF_3[jj])-K)*f_D for jj in range(self.Npatch)]
                
        self.XCOEF_1 = XCOEF_1_n
        self.XCOEF_2 = XCOEF_2_n
        self.XCOEF_3 = XCOEF_3_n
        
        [L1,B1,D1,K1] = self.avaliaExtremos()
        
        self.props()
        pass
    
    def avaliaExtremos(self):
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
            
        KG = abs(z_min)
        
#        print('Extremos: L = ' + str(L) + '\n          B = ' + str(B) + '\n          D = ' + str(D) + '\n          KG = ' + str(KG))
        
        return [L,B,D,KG]

#painel = GDF('aliv2.gdf')
#print(painel.propriedades.area_total)

#import matplotlib.pyplot as plt
#
#plt.figure()
#id_plot = 0
#vertices = painel.propriedades.vertices[id_plot]
#for fc in painel.propriedades.faces[id_plot]:
#    x = [vertices[fc[0]][0], vertices[fc[1]][0], vertices[fc[2]][0], vertices[fc[0]][0]]
#    y = [vertices[fc[0]][1], vertices[fc[1]][1], vertices[fc[2]][1], vertices[fc[0]][1]]
#    plt.plot(x,y,'-b')
#
#for ct in painel.propriedades.centroide[id_plot]:
#    plt.plot(ct[0],ct[1],'x')
#    
#plt.show()