# coding: utf-8

#!/usr/bin/env python

__author__ = "Renan Alves de Oliveira, Thiago dos Santos Pereira"
__license__ = "CC BY-NC 4.0"
__version__ = "1.0"
__maintainer__ = "Renan Alves de Oliveira"
__email__ = "fisica.renan@gmail.com"

##################
# Basic packages #
##################

import healpy as hp
import numpy as np
from numba import jit

################################################################
# Auxiliary functions - For generate norms, angles and sort it #
################################################################

def listanorms(nmax):
    """
    Function that gives all possible norms and angles. Perhaps it might contain duplicates.

    Args:
        nmax (int): Maximum possible norm.

    Returns:
        Array containing [sqrt(norm), nx, ny, nz, theta, and phi]. 
    """
    lista = [] #Create an empty list.
    nr = np.arange(0, 3*nmax + 1) #Range for nx and ny with include nx = ny = 0.
    nzr = np.arange(1, 3*nmax + 1) #Range for nz (not including nz = 0).
    for n in range(1, nmax + 1):
        for nx in nr:
            for ny in nr:
                for nz in nzr:
                    ns = nx**2 + ny**2 + nz**2
                    if 3*n**2 >= ns: #Condition to get all possible norms less or equal of 3*n**2.
                        thetan = np.arccos(nz/np.sqrt(ns)) #Calcule the polar angle.
                        if nx == 0:
                            if ny == 0:
                                phin = 0 #Azimuthal angle for this condition.
                                lista.append((np.sqrt(ns), nx, ny, nz, thetan, phin))
                            else:
                                phin = np.pi/2 #Azimuthal angle for this condition.
                                lista.append((np.sqrt(ns), nx, ny, nz, thetan, phin))
                        else:
                            phin = np.arctan(ny/nx) #Azimuthal angle for generic condition.
                            lista.append((np.sqrt(ns), nx, ny, nz, thetan, phin))
    lista = np.array(lista) #Transform this list in a numpy array.
    lista = np.array(sorted(lista, key=lambda norm: norm[0])) #Sort generated list in crescent order of sqrt(n).
    return lista

def unique_rows(a): #This function remove duplicates in the list.
    """
    Function that removes repeated arrays. For more details, check:
    
    https://stackoverflow.com/questions/31097247/remove-duplicate-rows-of-a-numpy-array

    Args:
        a (array): Specify an array that you want to removed duplicated rows.

    Returns:
        Array. 
    """
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# Before you start computing all maps and realizations, you'll need to set up the number of realizations you need
# and the cutoff for nmax that you are interested.

###############################################
# Parameters: Number of realizations and kmax #
###############################################

cutoff = int(input("Enter the value for nmax: "))
NaA = unique_rows(listanorms(cutoff)) #Norms and Angles.
NoR = int(input("Enter the number of realizations: ")) #NoR: Number of Realizations.
RA = np.array(np.random.randn(len(NaA), NoR, 8)) #Random Array for all realizations.

#####################
# FLRW - Background #
#####################

@jit(nogil = True)
def T(L, A, theta, phi, realization):
    """
    Function that evaluate temperature fluctuations in pixel space.

    Args:
        L (float): Lengh of a box.
        A (float): Amplitude of the CMB.
        theta (float): Polar angle in real space.
        phi (float): Azimuthal angle in real space.
        realization (int): Realization that this function will calculate. 

    Returns:
        Float number. 
    """
    inside_sum = [0.]
    for n in range(1, len(NaA) + 1):
        k = 2*np.pi*(NaA[n - 1, 0])/L
        power = np.sqrt(A*k**(- 3)) #sqrt of the FRLW Power Spectrum.
        cosgamma1 = np.cos(theta)*np.cos(NaA[n - 1, 4]) + np.sin(theta)*np.sin(NaA[n - 1, 4])*np.cos(phi - NaA[n - 1, 5])
        cosgamma2 = np.cos(theta)*np.cos(NaA[n - 1, 4]) - np.sin(theta)*np.sin(NaA[n - 1, 4])*np.cos(phi + NaA[n - 1, 5]) 
        cosgamma3 = np.cos(theta)*np.cos(NaA[n - 1, 4]) - np.sin(theta)*np.sin(NaA[n - 1, 4])*np.cos(phi - NaA[n - 1, 5]) 
        cosgamma4 = np.cos(theta)*np.cos(NaA[n - 1, 4]) + np.sin(theta)*np.sin(NaA[n - 1, 4])*np.cos(phi + NaA[n - 1, 5])
        #Phi_Real.
        term1R = RA[n - 1, realization, 0]*np.cos(k*cosgamma1)
        term2R = RA[n - 1, realization, 1]*np.cos(k*cosgamma2)
        term3R = RA[n - 1, realization, 2]*np.cos(k*cosgamma3)
        term4R = RA[n - 1, realization, 3]*np.cos(k*cosgamma4)
        #Phi_Imaginary.
        term1I = - RA[n - 1, realization, 4]*np.sin(k*cosgamma1)
        term2I = - RA[n - 1, realization, 5]*np.sin(k*cosgamma2)
        term3I = - RA[n - 1, realization, 6]*np.sin(k*cosgamma3)
        term4I = - RA[n - 1, realization, 7]*np.sin(k*cosgamma4)
        result = term1R + term1I + term2R + term2I + term3R + term3I + term4R + term4I
        inside_sum.append(2*power*result)    
    return sum(inside_sum)

def angles(nside):
    """
    Function that returns Polar (theta) and Azimuthal (phi) angles for Healpy (hp).

    Args:
        nside (int): NSIDE for Healpy.
        index (int): Recpective index for each pixel.

    Returns:
        Float number. 
    """    
    angles = np.ones((hp.pixelfunc.nside2npix(nside), 2))
    for i in range(hp.pixelfunc.nside2npix(nside)):
        theta, phi = hp.pixelfunc.pix2ang(nside, i)
        angles[i] = theta, phi
    return angles

def Pixels(nside, L, A):
    """
    Function that returns a Map in the Pixel space.

    Args:
        nside (int): NSIDE for Healpy.
        L (float): Lengh of a box.
        A (float): Amplitude of the CMB.
        
    Returns:
        Array. 
    """
    return np.array([[T(L, A, theta, phi, i) for theta, phi in angles(nside)] for i in range(0, NoR)])

def Alms(nside, L, A):
    """
    Alms for Healpy generate maps. Healpy will determine the lmax and mmax based in the value of nside.

    Args:
        nside (int): NSIDE for Healpy.
        L (int): Lengh of a box.
        A (int): Amplitude of the CMB.
        
    Returns:
        Complex array. 
    """
    return np.array([hp.map2alm(np.array([T(L, A, theta, phi, i) for theta, phi in angles(nside)])) for i in range(0, NoR)])
