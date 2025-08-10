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

from matplotlib.patches import Rectangle
import shapely
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.colors as mcolors
from shapely.geometry import Polygon, Point

plt.rcParams["font.size"]=12
plt.rcParams["axes.labelsize"]=12
plt.rcParams["figure.dpi"]=200
plt.rcParams['figure.figsize'] = [6, 5.5]

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')


class ModuleGeoGrid:

        def __init__(self, xmin=-2150, xmax=2150, ymin=-2150, ymax=2150, ngrid=100):
            self.Xmax = xmax
            self.Xmin = xmin
            self.Ymax = ymax
            self.Ymin = ymin
            self.ngrid = ngrid ## grid step (default 1 cm2)
            self.makeGrid()

        def makeGrid(self):
            ngrid = self.ngrid
            dl = (self.Xmax - self.Xmin) / ngrid
            xBounds = {}
            yBounds = {}
            k=0 ## unique dic key
            for i in range(ngrid):
                for j in range(ngrid):
                    xBounds[k] = np.array([self.Xmin + i*dl, self.Xmin + (i+1)*dl]) ## key: strip index, entry: strip bounds (xmin, xmax)
                    yBounds[k] = np.array([self.Ymin + j*dl, self.Ymin + (j+1)*dl]) ## key: strip index, entry: strip bounds (ymin, ymax)
                    k += 1
            # print("done")
            self.xBounds = xBounds
            self.yBounds = yBounds
            return self.xBounds, self.yBounds

class ModuleGeoRect1:

        def __init__(self, nstrips=10, xmin=-2120, xmax=2120, ymin=-2120, ymax=2120, L=3000, l=1240):
            self.nstrips = nstrips
            self.Xmax = xmax
            self.Xmin = xmin
            self.Ymax = ymax
            self.Ymin = ymin
            self.L = L
            self.l = l
            self.D = self.Xmax-self.Xmin
            print("Rectangular module dimensions: L = {}, l = {}, D = {}".format(self.L, self.l, self.D))
            self.l = (self.Xmax-self.Xmin) - self.L
            self.pitch = self.L/self.nstrips
            self.makeModules()
            self.makeModuleSet()

        def makeModule(self):
            moduleXmax = self.Xmin + self.L
            moduleXmin = self.Xmin
            moduleYmax = self.Ymax
            moduleYmin = self.Ymax - self.l
            return np.array([moduleXmin, moduleXmax]), np.array([moduleYmin, moduleYmax])

        def makeModules(self):
            ### module 1
            self.xBounds1, self.yBounds1 = self.makeModule()
            ### module 2 
            xBounds2 = self.xBounds1*np.cos(-np.pi/2) + self.yBounds1*np.sin(-np.pi/2)
            yBounds2 = -self.xBounds1*np.sin(-np.pi/2) + self.yBounds1*np.cos(-np.pi/2)
            self.xBounds2 = np.array([min(xBounds2), max(xBounds2)])
            self.yBounds2 = np.array([min(yBounds2), max(yBounds2)])
            ### module 3
            self.xBounds3 = self.xBounds1 + self.l
            self.yBounds3 = self.yBounds1 - self.L
            ### module 4
            self.xBounds4 = self.xBounds2 + self.L
            self.yBounds4 = self.yBounds2 + self.l

        def makeStrips(self):
            dL = self.L/self.nstrips
            xStripsBounds = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(self.nstrips):
                xStripsBounds[i] = np.array([self.Xmin + i*dL, self.Xmin + (i+1)*dL])
                yStripsBounds[i] = np.array([self.Ymax - self.l, self.Ymax])
            # print("done")
            return xStripsBounds, yStripsBounds
        
        def makeModuleSet(self):
            ### module 1
            self.xStripsBounds1, self.yStripsBounds1 = self.makeStrips()
            ### module 2
            xStripsBounds2 = { i : self.xStripsBounds1[i]*np.cos(-np.pi/2) + self.yStripsBounds1[i]*np.sin(-np.pi/2) for i in range(self.nstrips) }
            yStripsBounds2 = { i : -self.xStripsBounds1[i]*np.sin(-np.pi/2) + self.yStripsBounds1[i]*np.cos(-np.pi/2) for i in range(self.nstrips) }
            self.xStripsBounds2 = { i : np.array([min(xStripsBounds2[i]), max(xStripsBounds2[i])]) for i in range(self.nstrips) }
            self.yStripsBounds2 = { i : np.array([min(yStripsBounds2[i]), max(yStripsBounds2[i])]) for i in range(self.nstrips) }
            # ### module 3
            self.xStripsBounds3 = { i : self.xStripsBounds1[i] + self.l for i in range(self.nstrips) }
            self.yStripsBounds3 = { i : self.yStripsBounds1[i] - self.L for i in range(self.nstrips) }
            # ### module 4
            self.xStripsBounds4 = { i : self.xStripsBounds2[i] + self.L for i in range(self.nstrips) }
            self.yStripsBounds4 = { i : self.yStripsBounds2[i] + self.l for i in range(self.nstrips) }
            # print("done")


class ModuleGeoRect2:

        def __init__(self, nstrips=10, xmin=-2120, xmax=2120, ymin=-2120, ymax=2120, L=3000, l=1240):
            self.nstrips = nstrips
            self.Xmax = xmax
            self.Xmin = xmin
            self.Ymax = ymax
            self.Ymin = ymin
            self.D = self.Xmax-self.Xmin
            self.L = L
            self.l = l
            print("Rectangular module dimensions: L = {}, l = {}, D = {}".format(self.L, self.l, self.D))
            # self.l = (self.Xmax-self.Xmin) - self.L
            self.pitch = self.l/self.nstrips   ## l for rect2, L for rect1
            self.makeModules()
            self.makeModuleSet()

        def makeModule(self):
            moduleXmax = self.Xmin + self.L
            moduleXmin = self.Xmin
            moduleYmax = self.Ymax
            moduleYmin = self.Ymax - self.l
            return np.array([moduleXmin, moduleXmax]), np.array([moduleYmin, moduleYmax])

        def makeModules(self):
            ### module 1
            self.xBounds1, self.yBounds1 = self.makeModule()
            ### module 2 
            xBounds2 = self.xBounds1*np.cos(-np.pi/2) + self.yBounds1*np.sin(-np.pi/2)
            yBounds2 = -self.xBounds1*np.sin(-np.pi/2) + self.yBounds1*np.cos(-np.pi/2)
            self.xBounds2 = np.array([min(xBounds2), max(xBounds2)])
            self.yBounds2 = np.array([min(yBounds2), max(yBounds2)])
            ### module 3
            self.xBounds3 = self.xBounds1 + self.l
            self.yBounds3 = self.yBounds1 - self.L
            ### module 4
            self.xBounds4 = self.xBounds2 + self.L
            self.yBounds4 = self.yBounds2 + self.l

        def makeStrips(self):
            dl = self.l/self.nstrips
            xStripsBounds = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(self.nstrips):
                xStripsBounds[i] = np.array([self.Xmin, self.Xmin + self.L])
                yStripsBounds[i] = np.array([self.Ymax - self.l + i*dl, self.Ymax - self.l + (i+1)*dl])
            # print("done")
            return xStripsBounds, yStripsBounds
        
        def makeModuleSet(self):
            ### module 1
            self.xStripsBounds1, self.yStripsBounds1 = self.makeStrips()
            ### module 2
            xStripsBounds2 = { i : self.xStripsBounds1[i]*np.cos(-np.pi/2) + self.yStripsBounds1[i]*np.sin(-np.pi/2) for i in range(self.nstrips) }
            yStripsBounds2 = { i : -self.xStripsBounds1[i]*np.sin(-np.pi/2) + self.yStripsBounds1[i]*np.cos(-np.pi/2) for i in range(self.nstrips) }
            self.xStripsBounds2 = { i : np.array([min(xStripsBounds2[i]), max(xStripsBounds2[i])]) for i in range(self.nstrips) }
            self.yStripsBounds2 = { i : np.array([min(yStripsBounds2[i]), max(yStripsBounds2[i])]) for i in range(self.nstrips) }
            # ### module 3
            self.xStripsBounds3 = { i : self.xStripsBounds1[i] + self.l for i in range(self.nstrips) }
            self.yStripsBounds3 = { i : self.yStripsBounds1[i] - self.L for i in range(self.nstrips) }
            # ### module 4
            self.xStripsBounds4 = { i : self.xStripsBounds2[i] + self.L for i in range(self.nstrips) }
            self.yStripsBounds4 = { i : self.yStripsBounds2[i] + self.l for i in range(self.nstrips) }
            # print("done")


### TENTATIVE
class ModuleGeoTrapezoidal:
    def __init__(self, nstrips=10, xmin=-2120, xmax=2120, ymin=-2120, ymax=2120, H=1240):
        self.nstrips = nstrips
        self.Xmax = xmax
        self.Xmin = xmin
        self.Ymax = ymax
        self.Ymin = ymin
        self.L = (self.Xmax-self.Xmin) ## outer radius trapezoid length
        self.H = H ## trapezoid height
        self.l = self.L - 2*self.H ## inner radius trapezoid length
        print("Trapezoidal module dimensions: L = {}, l = {}, H = {}".format(self.L, self.l, self.H))
        self.outPitch = self.L/self.nstrips ## outer pitch of trapezoid module
        self.inPitch = self.l/self.nstrips ## inner pitch of trapezoid module
        self.makeStrips()

    def makeStrips(self):

        dL = self.outPitch 
        dl = self.inPitch
        nstrips = self.nstrips

        self.module1strips = {}  # Dict of polygons
        self.module2strips = {}  # Dict of polygons
        self.module3strips = {}  # Dict of polygons
        self.module4strips = {}  # Dict of polygons

        ### module 1
        for i in range(nstrips):
            y0 = self.Ymax - self.H
            y1 = self.Ymax
            x0 = self.Xmin + self.H + i*dl
            x1 = self.Xmin + i*dL
            x2 = self.Xmin + self.H + (i+1)*dl
            x3 = self.Xmin + (i+1)*dL

            # Trapezoid corners (bottom wider than top)
            poly = Polygon([
                (x1, y1),  # top left
                (x3, y1),  # top right
                (x2, y0),  # bottom right
                (x0, y0)   # bottom left
            ])
            self.module1strips[i] = poly

        ### module 2
        for i in range(nstrips):
            x0 = self.Xmin + self.H
            x1 = self.Xmin 
            y0 = self.Ymin + self.H + i*dl
            y1 = self.Ymin + i*dL
            y2 = self.Ymin + self.H + (i+1)*dl
            y3 = self.Ymin + (i+1)*dL

            # Trapezoid corners (bottom wider than top)
            poly = Polygon([
                (x1, y3),  # top left
                (x0, y2),  # top right
                (x0, y0),  # bottom right
                (x1, y1)   # bottom left
            ])
            self.module2strips[i+nstrips] = poly

        ### module 3
        for i in range(nstrips):
            y0 = self.Ymin + self.H
            y1 = self.Ymin
            x0 = self.Xmin + self.H + i*dl
            x1 = self.Xmin + i*dL
            x2 = self.Xmin + self.H + (i+1)*dl
            x3 = self.Xmin + (i+1)*dL

            # Trapezoid corners (bottom wider than top)
            poly = Polygon([
                (x0, y0),  # top left
                (x2, y0),  # top right
                (x3, y1),  # bottom right
                (x1, y1)   # bottom left
            ])
            self.module3strips[i+2*nstrips] = poly

        ### module 4
        for i in range(nstrips):
            x0 = self.Xmax - self.H
            x1 = self.Xmax
            y0 = self.Ymin + self.H + i*dl
            y1 = self.Ymin + i*dL
            y2 = self.Ymin + self.H + (i+1)*dl
            y3 = self.Ymin + (i+1)*dL

            # Trapezoid corners (bottom wider than top)
            poly = Polygon([
                (x0, y2),  # top left
                (x1, y3),  # top right
                (x1, y1),  # bottom right
                (x0, y0)   # bottom left
            ])
            self.module4strips[i+3*nstrips] = poly
        
        self.allStrips = self.module1strips | self.module2strips | self.module3strips | self.module4strips



class CRVHitMap:

    def __init__(self, inFilePath, nPOT=99995000, mode='front'):
        self.nPOT = nPOT
        self.outFileName = "{}.txt".format(Path(inFilePath).stem)
        self.file = uproot.open(inFilePath)
        self.fileKeys = self.file.keys()
        self.treeFront = self.file['treeFrontCRV']
        self.treeBack = self.file['treeBackCRV']
        self.treeFront_arr = self.treeFront.arrays()
        self.treeBack_arr = self.treeBack.arrays()
        self.getHitPos(mode)
        # self.getGeometry()
        
    def getHitPos(self, mode='front'):
        if mode=='front':
            treeFront_arr = self.treeFront_arr
            self.IsFirstHit = treeFront_arr.IsFirstHit
            self.RPCweights = treeFront_arr.RpcSensitivity[self.IsFirstHit==True]
            self.HitT = np.array(treeFront_arr.PositionT)[self.IsFirstHit==True]
            self.HitX = np.array(treeFront_arr.PositionX)[self.IsFirstHit==True]
            self.HitY = np.array(treeFront_arr.PositionY)[self.IsFirstHit==True]
            self.HitZ = np.array(treeFront_arr.PositionZ)[self.IsFirstHit==True]
            self.nhits = len(self.HitT)

    def getGeometry(self, nstrips=10, L=3000, l=1240, xmin=-2120, xmax=2120, ymin=-2120, ymax=2120, geo_type="rectangular strips 1"):
        self.geo_type = geo_type
        if geo_type=='grid':
            self.geometry = ModuleGeoGrid(ngrid=nstrips, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        elif geo_type=='rectangular strips 1':
            self.geometry = ModuleGeoRect1(nstrips=nstrips, L=L, l=l, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        elif geo_type=='rectangular strips 2':
            self.geometry = ModuleGeoRect2(nstrips=nstrips, L=L, l=l, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        elif geo_type=='trapezoidal strips':
            self.geometry = ModuleGeoTrapezoidal(nstrips=nstrips, H=L, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

    def getHitsPerChannel(self):
        self.POTperSec = 2.5e12 
        self.hitRateFactor = self.POTperSec / self.nPOT
        if self.geo_type=='rectangular strips 1' or self.geo_type=='rectangular strips 2':
            nstrips = self.geometry.nstrips
            self.module1Hits = {i : 0 for i in range(nstrips)}
            self.module2Hits = {i : 0 for i in range(nstrips)}
            self.module3Hits = {i : 0 for i in range(nstrips)}
            self.module4Hits = {i : 0 for i in range(nstrips)}
            for ihit in range(self.nhits):
                # iterate over hits
                xi = self.HitY[ihit]
                yi = self.HitZ[ihit]
                weight = self.RPCweights[ihit]
                # check if hit is in module
                if self.geometry.xBounds1[0]<=xi<self.geometry.xBounds1[1] and self.geometry.yBounds1[0]<=yi<self.geometry.yBounds1[1]:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds1[i][0]<=xi<self.geometry.xStripsBounds1[i][1] and self.geometry.yStripsBounds1[i][0]<=yi<self.geometry.yStripsBounds1[i][1]:
                            self.module1Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module
                elif self.geometry.xBounds2[0]<=xi<self.geometry.xBounds2[1] and self.geometry.yBounds2[0]<=yi<self.geometry.yBounds2[1]:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds2[i][0]<=xi<self.geometry.xStripsBounds2[i][1] and self.geometry.yStripsBounds2[i][0]<=yi<self.geometry.yStripsBounds2[i][1]:
                            self.module2Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module
                if self.geometry.xBounds3[0]<=xi<self.geometry.xBounds3[1] and self.geometry.yBounds3[0]<=yi<self.geometry.yBounds3[1]:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds3[i][0]<=xi<self.geometry.xStripsBounds3[i][1] and self.geometry.yStripsBounds3[i][0]<=yi<self.geometry.yStripsBounds3[i][1]:
                            self.module3Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module
                if self.geometry.xBounds4[0]<=xi<self.geometry.xBounds4[1] and self.geometry.yBounds4[0]<=yi<self.geometry.yBounds4[1]:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds4[i][0]<=xi<self.geometry.xStripsBounds4[i][1] and self.geometry.yStripsBounds4[i][0]<=yi<self.geometry.yStripsBounds4[i][1]:
                            self.module4Hits[i] += 1*weight*self.hitRateFactor
                            break

        elif self.geo_type=='trapezoidal strips':
            self.stripHits = {i: 0 for i in self.geometry.allStrips}
            for ihit in range(self.nhits):
                xi = self.HitY[ihit]  # Assuming same convention
                yi = self.HitZ[ihit]
                point = Point(xi, yi)
                weight = self.RPCweights[ihit]
                for i, polygon in self.geometry.allStrips.items():
                    if polygon.contains(point):
                        self.stripHits[i] += 1 * weight * self.hitRateFactor
                        break

        elif self.geo_type=='grid':
            ngrid2 = self.geometry.ngrid**2
            gridHits = {i : 0 for i in range(ngrid2)}
            xBounds = self.geometry.xBounds
            yBounds = self.geometry.yBounds
            for ihit in range(self.nhits):
                # iterate over hits
                xi = self.HitY[ihit]
                yi = self.HitZ[ihit]
                weight = self.RPCweights[ihit]
                # check if hit is in grid
                for i in range(ngrid2):
                    if xBounds[i][0]<=xi<xBounds[i][1] and yBounds[i][0]<=yi<yBounds[i][1]:
                        gridHits[i] += 1*weight*self.hitRateFactor
                        break
            self.gridHits = gridHits




def plot_strip_hitmap(crvmap, modules=[1, 2, 3, 4], figsize=(8, 6), plotHits=True):
    fig, ax = plt.subplots(figsize=figsize)

    # Collect all hit values for normalization
    all_hits = []
    for imod in modules:
        hits = getattr(crvmap, f"module{imod}Hits")
        all_hits.extend(hits.values())

    vmax = max(all_hits) if all_hits else 1e-6
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cmap = plt.cm.YlOrRd
    print("max_hit_rate = {:.2f} kHz/strip".format(vmax/1000))

    for imod in modules:
        hits = getattr(crvmap, f"module{imod}Hits")
        xbounds_dict = getattr(crvmap.geometry, f"xStripsBounds{imod}")
        ybounds_dict = getattr(crvmap.geometry, f"yStripsBounds{imod}")
        
        for i in range(len(hits)):
            if i not in xbounds_dict or i not in ybounds_dict:
                continue

            x0, x1 = xbounds_dict[i].tolist()
            y0, y1 = ybounds_dict[i].tolist()
            val = hits[i]
            color = cmap(norm(val))

            rect = Rectangle((x0, y0), x1 - x0, y1 - y0,
                             facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(rect)

    # Set axis limits explicitly
    all_x = []
    all_y = []
    for imod in modules:
        xbounds_dict = getattr(crvmap.geometry, f"xStripsBounds{imod}")
        ybounds_dict = getattr(crvmap.geometry, f"yStripsBounds{imod}")
        for i in range(len(xbounds_dict)):
            x0, x1 = xbounds_dict[i].tolist()
            y0, y1 = ybounds_dict[i].tolist()
            all_x.extend([x0, x1])
            all_y.extend([y0, y1])

    ax.set_xlim(min(all_x) - 100, max(all_x) + 100)
    ax.set_ylim(min(all_y) - 100, max(all_y) + 100)

    ax.set_xlabel("Y (mm)", size=14)
    ax.set_ylabel("Z (mm)", size=14)
    ax.set_aspect('equal')
    ax.set_title("Front CRV RPCs: rectangular strips")
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="Weighted Hit Rate (kHz/strip)")
    cb.ax.yaxis.label.set_size(14)
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'{x/1000:.0f}'))

    if plotHits==True:
        ax.scatter(crvmap.HitY, crvmap.HitZ, s=0.3, alpha=0.1)

    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()


### TENTATIVE
def plot_trapezoidal_strip_hitmap(crvmap, figsize=(8, 6), cmap_name='YlOrRd', plotHits=True):
    fig, ax = plt.subplots(figsize=figsize)

    # Get trapezoid shapes and hit values
    strips = crvmap.geometry.allStrips         # {index: shapely Polygon}
    hits = crvmap.stripHits                 # {index: hit count}

    # Normalize colors
    all_values = list(hits.values())
    vmax = max(all_values) if all_values else 1e-6
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cmap = plt.cm.get_cmap(cmap_name)

    for i, poly in strips.items():
        val = hits.get(i, 0)
        color = cmap(norm(val))

        coords = np.array(poly.exterior.coords)
        patch = MplPolygon(coords, closed=True, facecolor=color,
                           edgecolor='black', linewidth=0.5)
        ax.add_patch(patch)

    # Adjust limits to fit all polygons
    all_x = [x for poly in strips.values() for x, _ in poly.exterior.coords]
    all_y = [y for poly in strips.values() for _, y in poly.exterior.coords]
    ax.set_xlim(min(all_x) - 100, max(all_x) + 100)
    ax.set_ylim(min(all_y) - 100, max(all_y) + 100)

    ax.set_aspect('equal')
    ax.set_xlabel("Y (mm)", size=14)
    ax.set_ylabel("Z (mm)", size=14)
    ax.set_title("Front CRV RPCs: trapezoidal strips", size=14)
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label("Weighted Hit Rate (kHz/strip)", fontsize=14)
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'{x/1000:.0f}'))

    if plotHits==True:
        ax.scatter(crvmap.HitY, crvmap.HitZ, s=0.3, alpha=0.1)

    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_grid_hitmap(crvmap, figsize=(8, 6), plotHits=True):
    fig, ax = plt.subplots(figsize=figsize)

    # Collect all hit values for normalization
    all_hits = []
    hits = crvmap.gridHits
    all_hits.extend(hits.values())

    vmax = max(all_hits) if all_hits else 1e-6
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cmap = plt.cm.YlOrRd

    xbounds_dict = crvmap.geometry.xBounds
    ybounds_dict = crvmap.geometry.yBounds
    
    for i in range(len(hits)):
        if i not in xbounds_dict or i not in ybounds_dict:
            continue

        x0, x1 = xbounds_dict[i].tolist()
        y0, y1 = ybounds_dict[i].tolist()
        val = hits[i]
        color = cmap(norm(val))

        rect = Rectangle((x0, y0), x1 - x0, y1 - y0,
                            facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect)

    # Set axis limits explicitly
    all_x = []
    all_y = []
    for i in range(len(xbounds_dict)):
        x0, x1 = xbounds_dict[i].tolist()
        y0, y1 = ybounds_dict[i].tolist()
        all_x.extend([x0, x1])
        all_y.extend([y0, y1])

    ax.set_xlim(min(all_x) - 100, max(all_x) + 100)
    ax.set_ylim(min(all_y) - 100, max(all_y) + 100)

    ax.set_xlabel("Y (mm)", size=14)
    ax.set_ylabel("Z (mm)", size=14)
    ax.set_aspect('equal')
    ax.set_title("RPC-weighted front CRV hit distribution")
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=r"Weighted Hit Rate (kHz/dm$^2$)")
    cb.ax.yaxis.label.set_size(14)
    cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: f'{x/1000:.0f}'))

    if plotHits==True:
        ax.scatter(crvmap.HitY, crvmap.HitZ, s=0.3, alpha=0.1)

    # plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.show()