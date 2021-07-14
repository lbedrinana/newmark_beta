# 
#=========================================================          
#     newmark.py
#=========================================================
# Coded by Luis Bedriñana,
#          Universidad de Ingenieria y Tecnologia - UTEC    
#          Jun, 2020
# Version: 1.1
# 
# Script to analyze the dynamic linear response of 1-DOF subjected to GM record by the
# Newmark beta method
#
# INPUT:
#   per: period of the analyzed system
#   damp: Equivalent damping ratio of the analyzed system
#   GM_acc: Ground motion record read form a text file 
#
# OUTPUT:
#   Out_acce: Time history of the acceleration response
#   Out_disp: Time history of the displacement response 
#   peak_acc: Peak acceleration
#   peak_disp: peak displacement
#
# Imports 
import numpy as np
import matplotlib.pyplot as plt
import os
#
# FUNCTIONS
# ---------------------------------------------------------------
# Reading GM record from file
# be careful with the format of the input file
def readGM (fname):
    #fname should include the path
    fhand = open(fname,'r')
    GMlist=list()
    PGA = None
    print('Reading from input file...', fhand)
    ii = 0
    for line in fhand:
        line = line.rstrip()
        if ii == 0:
            npoints = line
        elif ii == 1:
            dtime = line
        else :
            GMlist.append(line)
            peaktemp = float(line)
            if PGA is None or abs(peaktemp) > PGA:
                PGA = abs(peaktemp)
        ii = ii + 1     
    return npoints, dtime, PGA, GMlist 
#
# ------------------------------------------------------------------
#
# ------------------------------ Newmark Method ---------------
# This formulation only works for linear systems
def newbet(acc_s, damp, ww, dtime, gamma, beta, dd, vv, aa):
    #
    cc = 2.*damp*ww
    # constants
    #
    bdv1 = gamma/(beta*dtime)
    bdv2 = 1.-(gamma/beta)
    bdv3 = 1.- (0.5*gamma/beta)
    bda1 = 1./(beta*dtime*dtime)
    bda2 = 1./(beta*dtime)
    bda3 = (0.5/beta)-1.
    p1 = bda1+(bdv1)*cc
    p2 = bda2 + ((gamma/beta) - 1.)*cc
    p3 = bda3 + ((0.5*gamma/beta) - 1.)*dtime*cc
    #equivalent linear stiffness
    kk = ww*ww+bdv1*cc+bda1
    ppn = -acc_s+p1*dd+p2*vv+p3*aa
    # New displacement 
    ddn=ppn/kk
    # New velocity
    vn=bdv1*(ddn-dd)+bdv2*vv+dtime*bdv3*aa
    # New in acceleration
    an=bda1*(ddn-dd)-bda2*vv-bda3*aa
    return ddn, vn, an
#
# ========================== MAIN ============================
# Setting initial values for input record
# Input motion path
fpath = '.\\Input_GM\\'
# GM text file name
fname = 'centro_NS.txt'
fname_t = fpath + fname
outRecord = readGM(fname_t)
# Input data for the system
per = 1.0
damp = 0.02
#
npoints = outRecord[0]
dtime = float(outRecord[1])
PGA = float(outRecord[2])
GMlist = outRecord[3]
GM_acc = np.array(GMlist, dtype='float32')
# Checking the input data
print('Number of points:', npoints)
print('Time interval of data:', dtime)
print (GM_acc)
print('Peak Ground acceleration:',PGA)
#
# Preparing the data for calculations
gamma = 0.5
beta = 0.25
pi = 3.14159265359
# initial response
dd = 0.
vv = 0.
aa = -GM_acc[0]
#
ww = 2.*pi/per
print(ww)
#
Xd = list()
Xat = list()
#
Xd_peak = None
Xat_peak = None
# looping in time
for xwave in GM_acc:
    acc_s = xwave
    # Calling the function
    #outNewbet = newbet(acc_s, damp, ww, dtime, gamma, beta, dd, vv, aa)
    dd, vv, aa = newbet(acc_s, damp, ww, dtime, gamma, beta, dd, vv, aa)
    # tuplet output ddn, vn, an, ant
    # updating values
    #dd = float(outNewbet[0])
    #vv = float(outNewbet[1])
    #aa = float(outNewbet[2])
    #
    # total acceleration
    aat = aa + acc_s
    #
    # saving output 
    Xd.append(dd)
    #Xa.append(aa)
    Xat.append(aat)
    if Xd_peak is None or abs(dd) > Xd_peak:
        Xd_peak = abs(dd)
    #
    if Xat_peak is None or abs(aat) > Xat_peak:
        Xat_peak = abs(aat)
#
# Peak values
print ('Peak displacement: ', Xd_peak)
print ('Peak total acceleration: ', Xat_peak)
# Output and figures 
#
# Path for output
fpath2 = '.\\Output_Newmark\\'
if not os.path.exists(fpath2):
    os.makedirs(fpath2)
#
time = np.arange(0.0, dtime*int(npoints), dtime)
# 
plt.title("Displacement response", fontsize = 20)
plt.ylabel("Displacement", fontsize=14)
plt.xlabel("Time", fontsize=14)
plt.plot(time,Xd)
plt.savefig(fpath2 + "TopDisp.png")
#
fig, axarr=plt.subplots(2,1,sharex=True)
axarr[0].plot(time,Xat)
axarr[0].set_title('Acceleration Response')
axarr[1].plot(time,GM_acc)
axarr[1].set_title('Ground Motion')
#
plt.savefig(fpath2 +"Acce.png")
# Writing csv output file
fname2 = fpath2 + "outwave.csv"
with open(fname2,'w') as fhand2:
    for i in range(len(Xat)):
        out_string = ""
        out_string += str(time[i])
        out_string += "," + str(Xat[i])
        out_string += "," + str(Xd[i])
        out_string += "\n"
        fhand2.write(out_string)   
fhand2.close()
plt.show()
#
