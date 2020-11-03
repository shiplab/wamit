import numpy as np

## Linear mooring restoration matrix according the paper by
## Natarajan et al. (1995). Analysis of moorings of berthed ship. Marine Structures, 9 (5), 481-499.
## Chernjawski, M. (1980). Mooring of Surface Vessels To Piers. Marine Technology, 17(1), 1–7.

# Fenders

kif = np.array([5257.0, 5257.0, 5257.0, 5257.0])
# kif = np.array([0.0, 0.0, 0.0, 0.0])
# Chocks positions (relative to the ship's CoG)
Cxf = np.array([-55.5, -52.0, 52.0, 55.5])
Cyf = np.array([ 23.2,  23.2, 23.2, 23.2])
Czf = np.array([  0.0,   0.0,  0.0,  0.0])
# Berth positions (relative to the ship's CoG)
Bxf = np.array([-55.5, -52.0, 52.0, 55.5])
Byf = np.array([ 25.2,  25.2, 25.2, 25.2])
Bzf = np.array([  0.0,   0.0,  0.0,  0.0])

# line vector
x0f = Bxf - Cxf
y0f = Byf - Cyf
z0f = Bzf - Czf

l0f = np.sqrt(x0f**2 + y0f**2 + z0f**2)

Axf = x0f / l0f
Ayf = y0f / l0f
Azf = z0f / l0f


# Lines
EA = np.array([475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97, 475.97])
EA_tail = np.array([57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29, 57.29])
# Chocks positions (relative to the ship's CoG)
Cx = np.array([-145.6, -145.6, -145.6, -139.7, -138.2, -134.0, -131.0, -109.3, -107.4, 92.4, 94.7, 122.3, 124.4, 133.0, 135.3, 146.4, 147.3, 147.3])
Cy = np.array([   2.9,    4.5,    8.1,   16.3,   16.9,   19.1,   19.9,   23.2,   23.2, 23.0, 22.8,  17.4,  16.9,  13.0,  12.0,   3.0,   2.0,  -2.0])
Cz = np.array([  13.5 ,  13.5,   13.5,   13.5,   13.5,   13.5,   13.5,   13.5,   13.5, 13.5, 13.5,  13.5,  13.5,  13.5,  13.5,  13.5,  13.5,  13.5])
# Distance from chock to bitt
Le = np.array([4.5, 4.5, 4.5, 27.0, 27.0, 7.8, 7.8, 7.8, 7.8, 42.7, 42.7, 26.0, 26.0, 8.9, 8.9, 9.5, 10.2, 10.2])
# Tail length
Ltail = np.array([11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0])
# Berth positions (relative to the ship's CoG)
Bx = np.array([-191.0, -191.0,-191.0, -147.0, -147.0, -102.0, -102.0, -54.0, -54.0, 54.0, 54.0, 147.0, 147.0, 152.0, 152.0, 191.0, 191.0, 191.0])
By = np.array([  63.2,   63.2,  63.2,   63.2,   63.2,   63.2,   63.2,  27.2,  27.2, 27.2, 27.2,  63.2,  63.2,  63.2,  63.2,  63.2,  63.2,  63.2])
Bz = np.array([   1.1,    1.1,   1.1,    1.1,    1.1,    1.1,    1.1,   1.1,   1.1,  1.1,  1.1,   1.1,   1.1,   1.1,   1.1,   1.1,   1.1,   1.1])

# line vector
x0 = Bx - Cx
y0 = By - Cy
z0 = Bz - Cz

l0 = np.sqrt(x0**2 + y0**2 + z0**2)

Ax = x0 / l0
Ay = y0 / l0
Az = z0 / l0

Inc_down = np.arctan(Az)*180/np.pi

L_total = l0 + Le

ki_wire = EA / (L_total-Ltail)
ki_tail = EA_tail / Ltail

ki = 1/(1/ki_wire + 1/ki_tail)


Cx = np.hstack((Cx, Cxf))
Cy = np.hstack((Cy, Cyf))
Cz = np.hstack((Cz, Czf))
Ax = np.hstack((Ax, Axf))
Ay = np.hstack((Ay, Ayf))
Az = np.hstack((Az, Azf))
ki = np.hstack((ki, kif))


Ci = []
Ki = []
Ai = []
Bi = []
Qi = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
for cx,cy,cz,ax,ay,az,kki in zip(Cx,Cy,Cz,Ax,Ay,Az,ki):
    CCi = np.array([[0, -cz, cy], [cz, 0, -cx], [-cy, cx, 0]])
    ki_11 = ax * ax
    ki_12 = ax * ay
    ki_13 = ax * az
    ki_21 = ay * ax
    ki_22 = ay * ay
    ki_23 = ay * az
    ki_31 = az * ax
    ki_32 = az * ay
    ki_33 = az * az
    KKi = kki*np.array([[ki_11, ki_12, ki_13], [ki_21, ki_22, ki_23], [ki_31, ki_32, ki_33]])
    AAi = np.dot(-KKi, CCi)
    BBi = np.dot(CCi, AAi)
    QQi = np.hstack((np.vstack((KKi, AAi.transpose())), np.vstack((AAi, BBi))))
    Ki.append(KKi)
    Ai.append(AAi)
    Bi.append(BBi)
    Qi += QQi


L_total = np.round(L_total, decimals=1)
Inc_down = np.round(Inc_down, decimals=0)

mr = open('C:\\Users\\dprat\\Documents\\PYTHON\\mooring_tests\\stiff.txt','w')

for qi in Qi:
    mr.write('    '.join('{:.5E}'.format(x) for x in qi)+'\n')
            
mr.write('\n')
mr.close()

x = 0
