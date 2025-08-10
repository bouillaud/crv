import math
import numpy as np
import sys
import os
import uproot
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import stats, optimize
import pylandau
from pathlib import Path

# from funs import *

plt.rcParams["font.size"]=10
plt.rcParams["axes.labelsize"]=12
plt.rcParams["figure.dpi"]=200
plt.rcParams['figure.figsize'] = [6, 4]

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

## define functions

def binx(bins):
    x = []
    for i in range(len(bins)-1):
        x.append( (bins[i]+bins[i+1])/2 )
    return np.array(x)

### get data ###

inFileName = sys.argv[1]

outFileStem = "{}".format(Path(inFileName).stem)

file = uproot.open(inFileName)

tree = file['tree_CTH']
tree_arr = tree.arrays()

x = np.array(tree_arr["x"])
y = np.array(tree_arr["y"])
z = np.array(tree_arr["z"])
px = np.array(tree_arr["px"])
py = np.array(tree_arr["py"])
pz = np.array(tree_arr["pz"])
t = np.array(tree_arr["t"])

seg = np.array(tree_arr["seg"])
hod = np.array(tree_arr["hod"])
cry = np.array(tree_arr["cry"])

ix = np.array(tree_arr["ix"])
iy = np.array(tree_arr["iy"])
iz = np.array(tree_arr["iz"])
ipx = np.array(tree_arr["ipx"])
ipy = np.array(tree_arr["ipy"])
ipz = np.array(tree_arr["ipz"])
it = np.array(tree_arr["it"])

ox = np.array(tree_arr["ox"])
oy = np.array(tree_arr["oy"])
oz = np.array(tree_arr["oz"])
opx = np.array(tree_arr["opx"])
opy = np.array(tree_arr["opy"])
opz = np.array(tree_arr["opz"])
ot = np.array(tree_arr["ot"])

pid = np.array(tree_arr["iPid"])
trackid = np.array(tree_arr["TrackId"])
nodeid = np.array(tree_arr["nodeId"])
dE = np.array(tree_arr["dE"])
Len = np.array(tree_arr["len"])
dEcorr = np.array(tree_arr["dEcorr"])
opid = np.array(tree_arr["oPid"])

iprocess = np.array(tree_arr["iProcess"])
# oprocess = tree_arr["oProcess"]

omin = -10000
ox = ox[ox>omin]
oy = oy[oy>omin]
oz = oz[oz>omin]
opx = opx[opx>omin]
opy = opy[opy>omin]
opz = opz[opz>omin]
ot = ot[ot>omin]

x0, y0, z0 = 6670, 0, 7650  ## CTH center / origin of CTH coord system

# cartesian coordinates in CTH system
xc = x-x0
yc = y-y0
zc = z-z0

# polar coordinates in CTH system
rho = ((y-y0)**2 + (z-z0)**2)**0.5
# irho = (iy**2 + iz**2)**0.5
# orho = (oy**2 + oz**2)**0.5
phi = np.arctan2(z-z0, y-y0)
# iphi = np.arctan2(iy, iz)
# ophi = np.arctan2(oy, oz)

p2 = np.sqrt(px**2 + py**2 + pz**2)
ip2 = np.sqrt(ipx**2 + ipy**2 + ipz**2)
op2 = np.sqrt(opx**2 + opy**2 + opz**2)

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
print("ptypes")
print(ptypes)
print("Nhits, Nprimary, Nparents:")
print(Nhits, Nprimary, Nparents)

optypes = {}
for idd in opid:
	if idd not in optypes.keys():
		optypes[idd] = 1
	else:
		optypes[idd] += 1
optypes = dict(sorted(optypes.items()))
print("optypes")
print(optypes)

npot = 4e7 ## number of protons on target

### constants
protonRate = 2.5e12 ## number of protons per second
bunchDt = 1316e-9 ## average bunch duration
npotPerBunch = 1.63e7 ## numbers of pot per bunch

simTime0 = npot/protonRate ## hits per sec method (old)
simTime =  (npot/npotPerBunch) * bunchDt  ## hits per bunch method

# print("{:.2e}".format(simTime0))
print("Sim time = {:.2e} s".format(simTime))

totalHitRate = Nhits/simTime
hitRatePerCounter = totalHitRate/256

print("Total hit rate = {:.2e} Hz".format(totalHitRate))
print("Avg hit rate per counter = {:.2e} Hz".format(hitRatePerCounter))

upOuter = seg[(hod==0) & (cry==0)]
upInner = seg[(hod==0) & (cry==1)]
doOuter = seg[(hod==1) & (cry==0)]
doInner = seg[(hod==1) & (cry==1)]

## plot hit rate 

# hTOT = plt.hist(seg, histtype='step', bins=64, label="all counters", weights=(1/simTime)*np.ones_like(seg))
# plt.xticks(np.linspace(0, 64, 7), np.linspace(0, 360, 7).astype(int))
# plt.xlabel(r"$\varphi$ $(^\circ)$")
# plt.ylabel("Hit rate per counter (Hz)")
# plt.legend(ncol=2)
# # plt.ylim(0, 3e8)

hUO = plt.hist(upOuter, histtype='step', bins=64, label="Upstream Outer", weights=(1/simTime)*np.ones_like(upOuter))
hUI = plt.hist(upInner, histtype='step', bins=64, label="Upstream Inner", weights=(1/simTime)*np.ones_like(upInner))
hDO = plt.hist(doOuter, histtype='step', bins=64, label="Downstream Outer", weights=(1/simTime)*np.ones_like(doOuter))
hDI = plt.hist(doInner, histtype='step', bins=64, label="Downstream Inner", weights=(1/simTime)*np.ones_like(doInner))
plt.xticks(np.linspace(0, 64, 7), np.linspace(0, 360, 7).astype(int))
plt.xlabel(r"$\varphi$ $(^\circ)$")
plt.ylabel("Hit rate per counter (Hz)")
plt.legend(ncol=2)
plt.ylim(0, 3.5e8)

hrpcUO = np.mean(hUO[0])
hrpcUI = np.mean(hUI[0])
hrpcDO = np.mean(hDO[0])
hrpcDI = np.mean(hDI[0])
# hrpcTOT = np.mean(hTOT[0])/4
print("hrpc upOuter = {:.2e} Hz".format(hrpcUO))
print("hrpc upInner = {:.2e} Hz".format(hrpcUI))
print("hrpc doOuter = {:.2e} Hz".format(hrpcDO))
print("hrpc doInner = {:.2e} Hz".format(hrpcDI))
# print("hrpc total = {:.2e} Hz".format(hrpcTOT))

# plt.show()
# plt.savefig("../plots/counterHitRates_{}.pdf".format(outFileStem), bbox_inches='tight')
