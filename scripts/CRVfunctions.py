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

def xBins(bins):
    return np.array([ (bins[i]+bins[i+1])/2 for i in range(len(bins)-1) ])

def getPDFfromDist(t_arr, nbins):
    pdf_t, bins_t = np.histogram(t_arr, nbins, density=True)
    xbins_t = xBins(bins_t)
    return pdf_t, xbins_t

def invCDFGen(cdf, xbins, nsamples):
    # draw y-value from uniform distribution
    u_arr = np.random.uniform(0, 1, nsamples)
    # find cdf point index of y-value
    indices = np.searchsorted(cdf, u_arr)
    # get x-value (x_i = F-1(y_i))
    x_dist = xbins[indices]
    return x_dist

def generateFromDist(t_arr, nbins, nsamples):
    pdf_t, xbins_t = getPDFfromDist(t_arr, nbins)
    norm = np.sum(pdf_t)
    print("pdf norm = {:.2e}".format(norm))
    cdf_t = np.cumsum(pdf_t) / norm
    t_dist = invCDFGen(cdf_t, xbins_t, nsamples)
    return t_dist

def getPDFfromDist_2D(x_data, y_data, nbins):
    hist2d, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=nbins, density=True)
    xbins, ybins = xBins(x_edges), xBins(y_edges)
    return hist2d, xbins, ybins

def invCDFGen_2D(cdf, tbins, sbins, nsamples, twoDshape):
    # draw y-value from uniform distribution
    u_arr = np.random.uniform(0, 1, nsamples)
    # find cdf point index of y-value
    indices = np.searchsorted(cdf, u_arr)
    # convert flat indices to 2D bin indices
    t_idx, s_idx = np.unravel_index(indices, twoDshape) 
    # get x-value (x_i = F-1(y_i))
    t_dist = tbins[t_idx]
    s_dist = sbins[s_idx]
    return t_dist, s_dist

def generateFromDist_2D(x_data, y_data, nbins, nsamples):
    hist2d, xbins, ybins = getPDFfromDist_2D(x_data, y_data, nbins)
    pdf = hist2d.flatten()
    norm = np.sum(pdf)
    print("pdf norm = {:.2e}".format(norm))
    cdf = np.cumsum(pdf) / norm
    x_dist, y_dist = invCDFGen_2D(cdf, xbins, ybins, nsamples, hist2d.shape)
    return x_dist, y_dist

def getPDFfromDist_3D(t_data, x_data, y_data, nbinst, nbinsxy):
    hist3d, edges = np.histogramdd([t_data, x_data, y_data], np.array([nbinst, nbinsxy, nbinsxy]), density=True)
    t_edges, x_edges, y_edges = edges
    tbins, xbins, ybins = xBins(t_edges), xBins(x_edges), xBins(y_edges)
    return hist3d, tbins, xbins, ybins

def invCDFGen_3D(cdf, tbins, sbins, rbins, nsamples, twoDshape):
    # draw y-value from uniform distribution
    u_arr = np.random.uniform(0, 1, nsamples)
    # find cdf point index of y-value
    indices = np.searchsorted(cdf, u_arr)
    # convert flat indices to 2D bin indices
    t_idx, s_idx, r_idx = np.unravel_index(indices, twoDshape) 
    # get x-value (x_i = F-1(y_i))
    t_dist = tbins[t_idx]
    s_dist = sbins[s_idx]
    r_dist = rbins[r_idx]
    return t_dist, s_dist, r_dist

def generateFromDist_3D(t_data, x_data, y_data, nbinst, nbinsxy, nsamples):
    hist3d, tbins, xbins, ybins = getPDFfromDist_3D(t_data, x_data, y_data, nbinst, nbinsxy)
    pdf = hist3d.flatten()
    norm = np.sum(pdf)
    print("pdf norm = {:.2e}".format(norm))
    cdf = np.cumsum(pdf) / norm
    t_dist, x_dist, y_dist = invCDFGen_3D(cdf, tbins, xbins, ybins, nsamples, hist3d.shape)
    return t_dist, x_dist, y_dist