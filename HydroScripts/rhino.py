import shutil
import subprocess
import re
import os

def LWL_ho(arq_gdf, ZPOT, close_rhino=True, path_gdf=[]):
    # LWL_ho function to create panel IRR=1 using Rhino 5
    # installed on standard folde "C:\Program Files\Rhinoceros 5.0 (64-bit)\System\Rhino.exe"
    #
    # Important: This works only with high order model with Y symmetry and unique ship surface
    #
    # Inputs:
    # arq_gdf = 'ship.gdf' Name of gdf file
    # ZPOT: WAMIT parameter to define the water line. Generally: ZPOT = KG - T
    # close_rhino: Flag to close (=True) or not (=False) Rhino programn after the confection of surface
    # local: path to the gdf file, if it is empty local = pwd

    if path_gdf != []:
        os.chdir(path_gdf)

    os.getcwd()
    
    #backup of original file
    shutil.copyfile(arq_gdf, arq_gdf + '_original')

    #removes the symmetry line to avoid prompts
    fileID = open(arq_gdf,'r')
    arqGDF = fileID.readlines()
    fileID.close()
    arqGDF[2] = '0 0\n'

    fileID = open(arq_gdf,'w')
    arqGDF = fileID.writelines(arqGDF)
    fileID.close()

    strComando = '_SetActiveViewport Front _SelAll _Section -500,' + str(-ZPOT) + ' 500,' + str(-ZPOT) + ' Enter ' +\
    '_FitCrv _DeleteInpu=Yes _Degree 3 _AngleTolerance 0.1 0.001 Enter _Divide 1 _SelNone '+\
    '_SelPt _CurveThroughPt _Degree 3 _CurveType Interpolated _Knots Uniform Enter _SelNone ' +\
    '_SelCrv _-Loft _Simplify None Enter _SelNone ' +\
    '_SelSrf _-Export ""' + os.getcwd() + '\\' + arq_gdf + '"" Y Enter'\

    if close_rhino==True:
        strComando = '"C:\\Program Files\\Rhinoceros 5.0 (64-bit)\\System\\Rhino.exe" /nosplash /notemplate /runscript="' +\
        strComando + ' _-RunScript (Rhino.DocumentModified False) -_exit" "' + arq_gdf + '"'

    else:
        strComando = '"C:\\Program Files\\Rhinoceros 5.0 (64-bit)\\System\\Rhino.exe" /nosplash /notemplate /runscript="' + strComando + '" "' + arq_gdf + '"'
    
    # print(strComando)
    
    subprocess.run(strComando)

    if path_gdf != []:
        os.chdir('..')


#debug
# LWL('ship_granel_001.gdf', -0.8, True)