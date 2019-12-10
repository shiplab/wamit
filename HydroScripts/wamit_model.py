import gdf_class
import numpy as np

def scaled(gdf_file, Lf, Bf, Df, T, KG, path_out='WAMIT'):
    ship = gdf_class.GDF(gdf_file)
    ship.gdf_scale(Lf, Bf, Df, T, KG)
    ship.gdf_write('ship.gdf', path_out, LWL_rhino=True)
    write_config_wam(path_out)
    write_fnames_wam(path_out)
    write_force_frc(path_out)
    write_ship_frc(path_out, ship.properties.Mass, ship.properties.Be)
    write_poten_pot(path_out, ship.properties.ZPOT)


def write_config_wam(path_out):
    ## Config.wam
    config = open(path_out + '/config.wam','w')
    config.write('ISOR=0\n')
    config.write('IRR(1)=1\n')
    config.write('ISOLVE=1\n')
    config.write('IQUAD=0\n')
    config.write('ILOG=1\n')
    config.write('IDIAG=0\n')
    config.write('IPERIN=2\n')
    config.write('IPEROUT=1\n')
    config.write('MONITR=0\n')
    config.write('NUMHDR=0\n')
    config.write('MAXSCR=12000\n')
    config.write('MAXMIT=12\n')
    config.write('MAXITT=40\n')
    config.write('IALTPOT=2\n')
    config.write('IALTFRC=3\n')
    config.write('IALTFRCN=2\n')
    config.write('IPOTEN=1\n')
    config.write('IFORCE=1\n')
    config.write('IPLTDAT=1\n')
    config.write('ILOWHI=1\n')
    config.write('KSPLIN=3\n')
    config.write('IQUADO=3\n')
    config.write('IQUADI=4\n')
    config.write('VMAXOPT9=-1\n')
    config.write('PANEL_SIZE=3\n')
    config.write('IGENMDS(1)=0\n')
    config.write('NEWMDS(1)=0\n')
    config.write('ITRIMWL=1\n')
    config.write('NCPU=20\n')
    config.write('RAMGBMAX=120\n')
    config.write('USERID_PATH=c:\\wamitv7061\\\n')

    config.close()

def write_fnames_wam(path_out):
    ## FNAMES.WAM

    fnames = open(path_out + '/fnames.wam','w')

    fnames.write('poten.pot\n')
    fnames.write('config.wam\n')
    fnames.write('force.frc\n')
    fnames.close()

def write_force_frc(path_out):
    ## FORCE.FRC
    force = open(path_out + '/force.frc','w')
    force.write('FORCE.FRC\n')
    force.write('1 1 1 1 0 0 0 1 0\n')
    force.write('1.025\n')
    force.write('ship.frc\n')
    force.write('0\n')
    force.write('0\n')
    force.write('0\n')
    force.write('0\n')
    force.close()

def write_ship_frc(path_out, M, Be):
    ## SHIP.FRC
    shipfrc = open(path_out + '/ship.frc','w')
    shipfrc.write('SHIP.FRC\n')
    shipfrc.write('1 1 1 1 0 0 0 1 0\n')
    shipfrc.write('1.025\n')
    shipfrc.write('0 0 0\n')
    shipfrc.write('1\n')
    for m1 in M:
        for m2 in m1:
            shipfrc.write("{:.5e}".format(m2))
            shipfrc.write(" ")
        shipfrc.write('\n')
    shipfrc.write('1\n')
    for m1 in Be:
        for m2 in m1:
            shipfrc.write("{:.5e}".format(m2))
            shipfrc.write(" ")
        shipfrc.write('\n')
    shipfrc.write('1\n')
    for m1 in np.zeros((6,6)):
        for m2 in m1:
            shipfrc.write("{:.5e}".format(m2))
            shipfrc.write(" ")
        shipfrc.write('\n')    
    shipfrc.write('0\n')
    shipfrc.write('0\n')
    shipfrc.close()

def write_poten_pot(path_out, ZPOT):
    ## POTEN.POT
    poten = open(path_out + '/poten.pot','w')
    poten.write('POTEN.POT\n')
    poten.write('-1\n')
    poten.write('1 1\n')
    poten.write('60\n')
    poten.write('0.01 0.05 0.1 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5 0.525 0.55 0.575 0.6 0.625 0.65 0.675 0.7 0.725 0.75 0.775 0.8 0.825 0.85 0.875 0.9 0.925 0.95 0.975 1 1.025 1.05 1.075 1.1 1.125 1.15 1.175 1.2 1.225 1.25 1.275 1.3 1.325 1.35 1.375 1.4 1.425 1.45 1.475 1.5 1.525 1.55\n')
    poten.write('25\n')
    poten.write('0 15 30 45 60 75 90 105 120 135 150 165 180 195 210 225 240 255 270 285 300 315 330 345 360\n')
    poten.write('1\n')
    poten.write('ship.gdf\n')
    poten.write('0 0 ' + str(ZPOT) + ' 0\n')
    poten.write('1 1 1 1 1 1\n')
    poten.close()

#debug
#scaled('ship.gdf', 200, 36, 15, 11, 10.2, path_out='WAMIT')