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

# matplotlib.rcParams['mathtext.fontset'] = 'custom'
# matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
# matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
# matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')
# matplotlib.rcParams['mathtext.fontset'] = 'stix'
# matplotlib.rcParams['font.family'] = 'STIXGeneral'
# matplotlib.pyplot.title(r'ABC123 vs $\mathrm{ABC123}^{123}$')

def getT2T3fromT1(R, L1, h1, gamma, h3):
    l = 2*R/(1+np.sqrt(2))
    d1 = np.sqrt(((L1-l)/2)**2 + h1**2)
    alpha1 = np.arctan((L1-l)/(2*h1)) + np.pi/2
    alpha2 = np.pi + gamma - alpha1
    e2 = d1*np.sin(alpha2 - np.pi/2)
    h2 = d1*np.cos(alpha2 - np.pi/2)
    L2 = (R + h3 - L1/2)/np.cos(gamma)
    f2 = L2 - e2 - l
    beta2 = np.arctan(f2/h2) + np.pi/2
    print("(alpha1, alpha2, beta2) = ({:.2f}, {:2f}, {:2f}) pi".format(alpha1/np.pi, alpha2/np.pi, beta2/np.pi))
    d2 = np.sqrt(f2**2 + h2**2)
    # d22 = h2/np.cos(beta2 - np.pi/2)
    L3 = 2*np.sqrt(d2**2 - h3**2) + l
    return l, d1, h2, d2, L2, L3

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


class ModuleGeoOctoTrapezoidal:
    def __init__(self, nstrips, xcenter, ycenter, R, L1, h1, h3, h5, gamma):
        self.nstrips = nstrips ## number of strips per module
        self.xcenter, self.ycenter = xcenter, ycenter ## center of inner octogon
        self.getParameters(R, L1, h1, h3, h5, gamma)
        # print("Module short base = {}".format(self.l))
        # print("Module long bases = {}".format(self.Ls))
        # print("Module heights = {}".format(self.Hs))
        # print("Module angles = {} pi".format(self.alphas/np.pi))
        # self.outPitches = self.Ls/self.nstrips ## outer pitches of trapezoid modules
        # self.inPitch = self.l/self.nstrips ## inner pitches of trapezoid modules
        # print("Module long pitches = {}".format(self.outPitches))
        # print("Module short pitch = {}".format(self.inPitch))
        self.makeStrips()

    def getParameters(self, R, L1, h1, h3, h5, gamma):

        self.R, self.L1, self.h1, self.h3, self.h5, self.gamma = R, L1, h1, h3, h5, gamma
        self.Da = 2*R + h1 + h5 ## total horizontal length
        self.Db = 2*(R + h3) ## total vertical length

        l, d1, h2, d2, L2, L3 = getT2T3fromT1(R, L1, h1, gamma, h3)
        print("l, d1, h2, d2, L2, L3 = ")
        print(l, d1, h2, d2, L2, L3)
        self.l = l ## trapeze short base
        self.d1, self.h2, self.d2, self.L2, self.L3 = d1, h2, d2, L2, L3

        l, d3, h4, d4, L4, L5 = getT2T3fromT1(R, L3, h3, gamma, h5)
        print("l, d3, h4, d4, L4, L5 = ")
        print(l, d3, h4, d4, L4, L5)
        self.d3, self.h4, self.d4, self.L4, self.L5 = d3, h4, d4, L4, L5

        self.ht1 = R + h1 - L3/2
        self.ht3 = R + h3 - L1/2
        self.htm1 = R + h5 - L3/2
        self.ht5 = R + h3 - L5/2
        
    def makeStrips(self):
        nstrips = self.nstrips
        xc, yc = self.xcenter, self.ycenter
        L1, L2, L3, L4, L5 = self.L1, self.L2, self.L3, self.L4, self.L5
        l = self.l
        dl = l/nstrips
        R = self.R
        d1, d2, d3, d4 = self.d1, self.d2, self.d3, self.d4
        h1, h2, h3, h4, h5 = self.h1, self.h2, self.h3, self.h4, self.h5
        ht1, htm1, ht3, ht5 = self.ht1, self.htm1, self.ht3, self.ht5

        self.module1strips = {}  # Dict of polygons
        self.module2strips = {}  
        self.module3strips = {}  
        self.module4strips = {}  
        self.module5strips = {} 
        self.module6strips = {}  
        self.module7strips = {}  
        self.module8strips = {}  
        self.module9strips = {}  
        self.module10strips = {}  
        self.module11strips = {}  
        self.module12strips = {}  

        ### module 1
        dX = L1/nstrips
        # dY = 0
        dx = dl
        # dy = 0
        for i in range(nstrips):
            x0 = xc - L1/2 + i*dX
            y0 = yc + R + h1
            x1 = xc - L1/2 + (i+1)*dX
            y1 = y0
            x2 = xc - l/2 + i*dx
            y2 = yc + R
            x3 = xc - l/2 + (i+1)*dx
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module1strips[i] = poly

        ### module 2
        alpha = -np.pi/4 ## module angle
        # alpha = -self.alphas[0] ## angle between module 1 and module 2
        dL2 = L2/nstrips
        dX = dL2*np.cos(alpha)
        dY = dL2*np.sin(alpha)
        dx = dl*np.cos(alpha)
        dy = dl*np.sin(alpha)
        for i in range(nstrips):
            x0 = xc + L1/2 + i*dX
            y0 = yc + R + h1 + i*dY
            x1 = xc + L1/2 + (i+1)*dX
            y1 = yc + R + h1 + (i+1)*dY
            x2 = xc + l/2 + i*dx
            y2 = yc + R + i*dy
            x3 = xc + l/2 + (i+1)*dx
            y3 = yc + R + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module2strips[i+nstrips] = poly

        ### module 3
        # alpha = -np.pi/4 ## module angle
        # dX = dLs[1]*np.cos(alpha)
        dY = -L3/nstrips
        # dx = dl*np.cos(alpha)
        dy = -dl
        for i in range(nstrips):
            x0 = xc + R + h3
            y0 = yc + L3/2 + i*dY
            x1 = x0
            y1 = yc + L3/2 + (i+1)*dY
            x2 = xc + R
            y2 = yc + l/2 + i*dy
            x3 = x2
            y3 = yc + l/2 + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module3strips[i+2*nstrips] = poly

        ### module 4
        alpha = np.pi/4 ## module angle
        dL4 = L4/nstrips
        dX = dL4*np.cos(alpha)
        dY = dL4*np.sin(alpha)
        dx = dl*np.cos(alpha)
        dy = dl*np.sin(alpha)
        for i in range(nstrips):
            x0 = xc + L5/2 + i*dX
            y0 = yc - R - h5 + i*dY
            x1 = xc + L5/2 + (i+1)*dX
            y1 = yc - R - h5 + (i+1)*dY
            x2 = xc + l/2 + i*dx
            y2 = yc - R + i*dy
            x3 = xc + l/2 + (i+1)*dx
            y3 = yc - R + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module4strips[i+3*nstrips] = poly
        
        ### module 5
        dX = L5/nstrips
        # dY = 0
        dx = dl
        # dy = 0
        for i in range(nstrips):
            x0 = xc - L5/2 + i*dX
            y0 = yc - R - h5
            x1 = xc - L5/2 + (i+1)*dX
            y1 = y0
            x2 = xc - l/2 + i*dx
            y2 = yc - R
            x3 = xc - l/2 + (i+1)*dx
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module5strips[i+4*nstrips] = poly
        
        ### module 6
        alpha = 3*np.pi/4 ## module angle
        dL4 = L4/nstrips
        dX = dL4*np.cos(alpha)
        dY = dL4*np.sin(alpha)
        dx = dl*np.cos(alpha)
        dy = dl*np.sin(alpha)
        for i in range(nstrips):
            x0 = xc - L5/2 + i*dX
            y0 = yc - R - h5 + i*dY
            x1 = xc - L5/2 + (i+1)*dX
            y1 = yc - R - h5 + (i+1)*dY
            x2 = xc - l/2 + i*dx
            y2 = yc - R + i*dy
            x3 = xc - l/2 + (i+1)*dx
            y3 = yc - R + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module6strips[i+5*nstrips] = poly

        ### module 7
        # alpha = -np.pi/4 ## module angle
        # dX = dLs[1]*np.cos(alpha)
        dY = -L3/nstrips
        # dx = dl*np.cos(alpha)
        dy = -dl
        for i in range(nstrips):
            x0 = xc - R - h3
            y0 = yc + L3/2 + i*dY
            x1 = x0
            y1 = yc + L3/2 + (i+1)*dY
            x2 = xc - R
            y2 = yc + l/2 + i*dy
            x3 = x2
            y3 = yc + l/2 + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module7strips[i+6*nstrips] = poly

        ### module 8
        alpha = -3*np.pi/4 ## module angle
        dL2 = L2/nstrips
        dX = dL2*np.cos(alpha)
        dY = dL2*np.sin(alpha)
        dx = dl*np.cos(alpha)
        dy = dl*np.sin(alpha)
        for i in range(nstrips):
            x0 = xc - L1/2 + i*dX
            y0 = yc + R + h1 + i*dY
            x1 = xc - L1/2 + (i+1)*dX
            y1 = yc + R + h1 + (i+1)*dY
            x2 = xc - l/2 + i*dx
            y2 = yc + R + i*dy
            x3 = xc - l/2 + (i+1)*dx
            y3 = yc + R + (i+1)*dy

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)   # bottom left
            ])
            self.module8strips[i+7*nstrips] = poly

        ### module 9 (triangular)
        # alpha = np.pi/4 ## module angle
        Lx = ht3
        Ly = ht1
        dLx = Lx/self.nstrips
        dLy = Ly/self.nstrips
        for i in range(nstrips):
            x0 = xc - R - h3 + i*dLx
            y0 = yc + L3/2 + i*dLy
            x1 = xc - R - h3 + (i+1)*dLx
            y1 = yc + L3/2 + (i+1)*dLy
            x2 = x0
            y2 = yc + R + h1
            x3 = x1
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)
            ])
            self.module9strips[i+8*nstrips] = poly

        ### module 10 (triangular)
        # alpha = np.pi/4 ## module angle
        Lx = ht3
        Ly = ht1
        dLx = Lx/self.nstrips
        dLy = Ly/self.nstrips
        for i in range(nstrips):
            x0 = xc + R + h3 - i*dLx
            y0 = yc + L3/2 + i*dLy
            x1 = xc + R + h3 - (i+1)*dLx
            y1 = yc + L3/2 + (i+1)*dLy
            x2 = x0
            y2 = yc + R + h1
            x3 = x1
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)
            ])
            self.module10strips[i+9*nstrips] = poly

        ### module 11 (triangular)
        # alpha = np.pi/4 ## module angle
        Lx = ht5
        Ly = htm1
        dLx = Lx/self.nstrips
        dLy = Ly/self.nstrips
        for i in range(nstrips):
            x0 = xc - R - h3 + i*dLx
            y0 = yc - L3/2 - i*dLy
            x1 = xc - R - h3 + (i+1)*dLx
            y1 = yc - L3/2 - (i+1)*dLy
            x2 = x0
            y2 = yc - R - h5
            x3 = x1
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)
            ])
            self.module11strips[i+10*nstrips] = poly

        ### module 12 (triangular)
        # alpha = np.pi/4 ## module angle
        Lx = ht5
        Ly = htm1
        dLx = Lx/self.nstrips
        dLy = Ly/self.nstrips
        for i in range(nstrips):
            x0 = xc + R + h3 - i*dLx
            y0 = yc - L3/2 - i*dLy
            x1 = xc + R + h3 - (i+1)*dLx
            y1 = yc - L3/2 - (i+1)*dLy
            x2 = x0
            y2 = yc - R - h5
            x3 = x1
            y3 = y2

            poly = Polygon([
                (x0, y0),  # top left
                (x1, y1),  # top right
                (x3, y3),  # bottom right
                (x2, y2)
            ])
            self.module12strips[i+11*nstrips] = poly

        octoDicts = [self.module1strips, self.module2strips, self.module3strips, self.module4strips,
                 self.module5strips, self.module6strips, self.module7strips, self.module8strips,
        ]
        trianDicts = [self.module9strips, self.module10strips, self.module11strips, self.module12strips]
        allDicts = octoDicts.copy()
        allDicts.extend(trianDicts)
        
        self.allStrips = {}
        self.octoStrips = {}
        self.trianStrips = {}
        for d in allDicts:
            self.allStrips |= d
        for d in octoDicts:
            self.octoStrips |= d
        for d in trianDicts:
            self.trianStrips |= d



class ModuleGeoRect2b:

        def __init__(self, nstrips=10, xmin=-2120, xmax=2120, ymin=-2120, ymax=2120, Ls=[3000, 3000, 3000, 2800], ls=[1040, 1240, 1240, 1240]):
            self.nstrips = nstrips
            self.Xmax = xmax
            self.Xmin = xmin
            self.Ymax = ymax
            self.Ymin = ymin
            self.D = self.Xmax-self.Xmin
            self.Ls = np.array(Ls)
            self.ls = np.array(ls)
            # print("Rectangular module dimensions: L = {}, l = {}, D = {}".format(self.L, self.l, self.D))
            # self.pitch = self.l/self.nstrips   ## l for rect2, L for rect1
            self.makeStrips()

        def makeStrips(self):
            nstrips = self.nstrips
            dls = self.ls/nstrips
            halfLs = self.Ls/2
            yoff = (self.Ymax-self.Ymin) - self.Ls[1] - self.ls[0]
            self.yoff = yoff
            print("y_off = {}".format(yoff))
            ### module 1
            xStripsBounds1 = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds1 = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(nstrips):
                xStripsBounds1[i] = np.array([self.Xmin, self.Xmin + halfLs[0]])
                yStripsBounds1[i] = np.array([self.Ymax - self.ls[0] + i*dls[0], self.Ymax - self.ls[0] + (i+1)*dls[0]])
                xStripsBounds1[i+nstrips] = np.array([self.Xmin + halfLs[0], self.Xmin + 2*halfLs[0]])
                yStripsBounds1[i+nstrips] = yStripsBounds1[i]
            self.xStripsBounds1, self.yStripsBounds1 = xStripsBounds1, yStripsBounds1
            ### module 2
            xStripsBounds2 = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds2 = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(nstrips):
                xStripsBounds2[i] = np.array([self.Xmin + i*dls[1], self.Xmin + (i+1)*dls[1]])
                yStripsBounds2[i] = np.array([self.Ymin + yoff, self.Ymin + yoff + halfLs[1]])
                xStripsBounds2[i+nstrips] = xStripsBounds2[i]
                yStripsBounds2[i+nstrips] = np.array([self.Ymin + yoff + halfLs[1], self.Ymin + yoff + 2*halfLs[1]])
            self.xStripsBounds2, self.yStripsBounds2 = xStripsBounds2, yStripsBounds2
            ### module 3
            xStripsBounds3 = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds3 = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(nstrips):
                xStripsBounds3[i] = np.array([self.Xmax - 2*halfLs[2], self.Xmax - halfLs[2]])
                yStripsBounds3[i] = np.array([self.Ymin + yoff + i*dls[2], self.Ymin + yoff + (i+1)*dls[2]])
                xStripsBounds3[i+nstrips] = np.array([self.Xmax - halfLs[2], self.Xmax])
                yStripsBounds3[i+nstrips] = yStripsBounds3[i]
            self.xStripsBounds3, self.yStripsBounds3 = xStripsBounds3, yStripsBounds3
            ### module 4
            xStripsBounds4 = {} ## key: strip index, entry: strip bounds (xmin, xmax)
            yStripsBounds4 = {} ## key: strip index, entry: strip bounds (ymin, ymax)
            for i in range(nstrips):
                xStripsBounds4[i] = np.array([self.Xmax - self.ls[3] + i*dls[3], self.Xmax - self.ls[3] + (i+1)*dls[3]])
                yStripsBounds4[i] = np.array([self.Ymax - 2*halfLs[3] , self.Ymax - halfLs[3]])
                xStripsBounds4[i+nstrips] = xStripsBounds4[i]
                yStripsBounds4[i+nstrips] = np.array([self.Ymax - halfLs[3], self.Ymax])
            self.xStripsBounds4, self.yStripsBounds4 = xStripsBounds4, yStripsBounds4

            self.xAllStrips = self.xStripsBounds1 | self.xStripsBounds2 | self.xStripsBounds3 | self.xStripsBounds4
            self.yAllStrips = self.yStripsBounds1 | self.yStripsBounds2 | self.yStripsBounds3 | self.yStripsBounds4


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
        elif geo_type=='rectangular strips 2b':
            self.geometry = ModuleGeoRect2b(nstrips=nstrips, Ls=L, ls=l, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    
    # def getNewGeometry(self, nstrips=10, xcenter=0, ycenter=200, Ls=[1800, 1600, 1600], l=730, Hs=[1000, 1100, 1100, 1100]):
    #     self.geometry = ModuleGeoOctoTrapezoidal(nstrips, xcenter, ycenter, Ls, l, Hs)
    #     self.geo_type = "octogonal trapezoidal"
    # def getNewGeometry(self, nstrips=10, xcenter=0, ycenter=200, Da=4300, Db=4300, l=730, L1=1800, H1=1000, gamma=np.pi/4):
    #     self.geometry = ModuleGeoOctoTrapezoidal(nstrips, xcenter, ycenter, Da, Db, l, L1, H1, gamma)
    #     self.geo_type = "octogonal trapezoidal"
    def getNewGeometry(self, nstrips=10, center=[0, 200], R=900, L1=2000, h1=1000, h3=1200, h5=1500, gamma=np.pi/4):
        self.geometry = ModuleGeoOctoTrapezoidal(nstrips, center[0], center[1], R, L1, h1, h3, h5, gamma)
        self.geo_type = "octogonal trapezoidal"

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
                elif self.geometry.xBounds3[0]<=xi<self.geometry.xBounds3[1] and self.geometry.yBounds3[0]<=yi<self.geometry.yBounds3[1]:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds3[i][0]<=xi<self.geometry.xStripsBounds3[i][1] and self.geometry.yStripsBounds3[i][0]<=yi<self.geometry.yStripsBounds3[i][1]:
                            self.module3Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module
                elif self.geometry.xBounds4[0]<=xi<self.geometry.xBounds4[1] and self.geometry.yBounds4[0]<=yi<self.geometry.yBounds4[1]:
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

        elif self.geo_type=='octogonal trapezoidal':
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
        
        elif self.geo_type=='rectangular strips 2b':
            geo = self.geometry
            yoff = geo.yoff
            nstrips = 2*geo.nstrips ## n strips per module includes 2 PCBs
            self.module1Hits = {i : 0 for i in range(nstrips)}
            self.module2Hits = {i : 0 for i in range(nstrips)}
            self.module3Hits = {i : 0 for i in range(nstrips)}
            self.module4Hits = {i : 0 for i in range(nstrips)}
            for ihit in range(self.nhits):
                # iterate over hits
                xi = self.HitY[ihit]
                yi = self.HitZ[ihit]
                weight = self.RPCweights[ihit]
                # check if hit is in module 1
                if geo.Xmin<=xi<(geo.Xmin+geo.Ls[0]) and (geo.Ymax-geo.ls[0])<=yi<geo.Ymax:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds1[i][0]<=xi<self.geometry.xStripsBounds1[i][1] and self.geometry.yStripsBounds1[i][0]<=yi<self.geometry.yStripsBounds1[i][1]:
                            self.module1Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module 2
                elif geo.Xmin<=xi<(geo.Xmin+geo.ls[1]) and (geo.Ymin+yoff)<=yi<(geo.Ymin+geo.Ls[1]+yoff):
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds2[i][0]<=xi<self.geometry.xStripsBounds2[i][1] and self.geometry.yStripsBounds2[i][0]<=yi<self.geometry.yStripsBounds2[i][1]:
                            self.module2Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module 3
                elif (geo.Xmax-geo.Ls[2])<=xi<geo.Xmax and (geo.Ymin+yoff)<=yi<(geo.Ymin+geo.ls[2]+yoff):
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds3[i][0]<=xi<self.geometry.xStripsBounds3[i][1] and self.geometry.yStripsBounds3[i][0]<=yi<self.geometry.yStripsBounds3[i][1]:
                            self.module3Hits[i] += 1*weight*self.hitRateFactor
                            break
                # check if hit is in module 4
                elif (geo.Xmax-geo.ls[3])<=xi<geo.Xmax and (geo.Ymax-geo.Ls[3])<=yi<geo.Ymax:
                    for i in range(nstrips):
                        # check if hit is in strip
                        if self.geometry.xStripsBounds4[i][0]<=xi<self.geometry.xStripsBounds4[i][1] and self.geometry.yStripsBounds4[i][0]<=yi<self.geometry.yStripsBounds4[i][1]:
                            self.module4Hits[i] += 1*weight*self.hitRateFactor
                            break




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