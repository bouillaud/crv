import math
import numpy as np
import sys
import os
import uproot
import matplotlib.pyplot as plt
from scipy import stats, optimize
import pylandau

# from funs import *

plt.rcParams["font.size"]=10
plt.rcParams["axes.labelsize"]=12

### get data ###

inFileName = "MC7_CTHana0.root"
outFileName = "{}.txt".format(inFileName)
file = uproot.open(inFileName)

tree = file['tree_CTH']
tree_arr = tree.arrays()

x = tree_arr["x"]
y = tree_arr["y"]
z = tree_arr["z"]
px = tree_arr["px"]
py = tree_arr["py"]
pz = tree_arr["pz"]
t = tree_arr["t"]

ix = tree_arr["ix"]
iy = tree_arr["iy"]
iz = tree_arr["iz"]
ipx = tree_arr["ipx"]
ipy = tree_arr["ipy"]
ipz = tree_arr["ipz"]
it = tree_arr["it"]

ox = tree_arr["ox"]
oy = tree_arr["oy"]
oz = tree_arr["oz"]
opx = tree_arr["opx"]
opy = tree_arr["opy"]
opz = tree_arr["opz"]
ot = tree_arr["ot"]

pid = tree_arr["Pid"]
trackid = tree_arr["TrackId"]
nodeid = tree_arr["nodeId"]
dE = tree_arr["dE"]
Len = tree_arr["len"]
dEcorr = tree_arr["dEcorr"]
opid = tree_arr["oPid"]

iprocess = tree_arr["iProcess"]
# oprocess = tree_arr["oProcess"]

omin = -10000
ox = ox[ox>omin]
oy = oy[oy>omin]
oz = oz[oz>omin]
opx = opx[opx>omin]
opy = opy[opy>omin]
opz = opz[opz>omin]
ot = ot[ot>omin]

####

Nhits = len(x)
Nprimary = len(ix)
Nparents = len(ox)

ptypes = {}
for idd in pid:
	if idd not in ptypes.keys():
		ptypes[idd] = 1
	else:
		ptypes[idd] += 1
ptypes = dict(sorted(ptypes.items()))
print(ptypes)


print(Nhits, Nprimary, Nparents)

rho = (y**2 + z**2)**0.5
irho = (iy**2 + iz**2)**0.5
orho = (oy**2 + oz**2)**0.5

# part = pid[1]

# plt.scatter(y[pid==part], z[pid==part], s=8+np.log10(dE[pid==part]), c=x[pid==part])
# plt.show()

# plt.scatter(iy[pid==part], iz[pid==part], c=ix[pid==part])
# plt.show()

fig1, ax1 = plt.subplots()
sc1 = ax1.scatter(z, y, s=8+np.log10(dE), c=x)
ax1.set_ylabel('y [mm]')
ax1.set_xlabel('z [mm]')
fig1.colorbar(sc1, label='x [mm]')
ax1.set_title("CTH hits")

fig3, ax3 = plt.subplots()
sc3 = ax3.scatter(rho, x, s=8+np.log10(dE))
ax3.set_ylabel('rho [mm]')
ax3.set_xlabel('x [mm]')
ax3.set_xlim(6000, 9500)
ax3.set_ylim(2000, 10000)
ax3.set_title("CTH hits")
# fig3.colorbar(sc3, label='x [mm]')

fig2, ax2 = plt.subplots()
sc2 = ax2.scatter(iz, iy, s=8+np.log10(dE), c=ix)
ax2.set_ylabel('y [mm]')
ax2.set_xlabel('z [mm]')
fig2.colorbar(sc2, label='x [mm]')
ax2.set_title("Primary particle origins")

fig5, ax5 = plt.subplots()
sc5 = ax5.scatter(irho, ix, s=8+np.log10(dE))
ax5.set_ylabel('rho [mm]')
ax5.set_xlabel('x [mm]')
ax5.set_title("Primary particle origins")
ax5.set_xlim(6000, 9500)
ax5.set_ylim(2000, 10000)
# fig3.colorbar(sc3, label='x [mm]')

fig6, ax6 = plt.subplots()
sc6 = ax6.scatter(oz, oy, c=ox, s=8)
ax6.set_ylabel('y [mm]')
ax6.set_xlabel('z [mm]')
ax6.set_title("Parent particle origins")
fig6.colorbar(sc6, label='x [mm]')

fig4, ax4 = plt.subplots()
sc4 = ax4.scatter(orho, ox, s=8)
ax4.set_ylabel('rho [mm]')
ax4.set_xlabel('x [mm]')
ax4.set_title("Parent particle origins")
ax4.set_xlim(6000, 9500)
ax4.set_ylim(2000, 10000)
# fig4.colorbar(sc4, label='x [mm]')ls




plt.show()


