# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:05:46 2019

@author: dprat
"""
import os
import analysis_wamit
import matplotlib.pyplot as plt

#pastas = ['Calado_10_1\Amazon_L228B400T101D15',
#         'Calado_10_1\BWAustria_L226B366T101D15',
#        'Calado_10_1\DanCisne_L207B322T101D15',
#        'Calado_10_1\Saltram_L225B366T101D15',
#        'Calado_10_1\Steena_L182B400T101D15',
#        'Calado_10_1\Tanker50k_L211B323T101D15',
#        'Calado_8_7\Amazon_L228B400T087D15',
#        'Calado_8_7\BWAustria_L226B366T087D15',
#        'Calado_8_7\DanCisne_L207B322T087D15',
#        'Calado_8_7\Saltram_L225B366T087D15',
#        'Calado_8_7\Steena_L182B400T087D15',
#        'Calado_8_7\Tanker50k_L211B323T087D15',
#        'Calado_Ballasted\Amazon_L228B400T071D15',
#        'Calado_Ballasted\BWAustria_L226B366T069D15',
#        'Calado_Ballasted\DanCisne_L207B322T078D15',
#        'Calado_Ballasted\Saltram_L225B366T069D15',
#        'Calado_Ballasted\Steena_L182B400T076D15',
#        'Calado_Ballasted\Tanker50k_L211B323T075D15']

pastas = ['calado_11\Amazon_L228B400T115D15',
'calado_11\BWAustria_L226B366T112D15',
'calado_11\DanCisne_L207B322T117D15',
'calado_11\Saltram_L225B366T112D15',
'calado_11\Steena_L182B400T115D15',
'calado_11\Tanker50k_L211B323T117D15',
'calado_11_8\Amazon_L228B400T119D15',
'calado_11_8\BWAustria_L226B366T118D15',
'calado_11_8\DanCisne_L207B322T119D15',
'calado_11_8\Saltram_L225B366T119D15',
'calado_11_8\Tanker50k_L211B323T119D15',
'calado_13\Amazon_L228B400T120D15',
'calado_13\DanCisne_L207B322T130D15',
'calado_13\Saltram_L225B366T120D15',
'calado_13\Tanker50k_L211B323T126D15',
'calado_9_4\Amazon_L228B400T090D15',
'calado_9_4\BWAustria_L226B366T090D15',
'calado_9_4\DanCisne_L207B322T090D15',
'calado_9_4\Saltram_L225B366T090D15',
'calado_9_4\Steena_L182B400T090D15',
'calado_9_4\Tanker50k_L211B323T090D15']

for pst in pastas:
    os.chdir(pst)
#    m = analysis_wamit.drift_forces_momentum(1)
    r = analysis_wamit.raos(1)
    os.chdir('..')
    nome = pst.split('\\')
    plt.savefig(nome[1])
    os.chdir('..')
    plt.close('all')
    
    
    