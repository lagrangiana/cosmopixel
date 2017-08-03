
***
# Readme
***

It is recomended to install most of all Python packages provided by [Anaconda](https://www.continuum.io/downloads). This includes Numba which will acceleted your calculations. The section below provides some mathematical aspects for this code. Aditionaly, to run this code script properly, you will need:

- [Healpy](https://anaconda.org/conda-forge/healpy)

***
# Temperature fluctuations on Pixel in FRLW
***
In this work, we are interested in large scale effects, which is dominated by the Sachs Wolfe effect,
$$\frac{\Delta\mathrm{T}}{\mathrm{T}}\left(\mathbf{r}\right)	=	\frac{1}{3}\Phi\left(\mathbf{r}\right),$$
where $\frac{\Delta\mathrm{T}}{\mathrm{T}}$ is the temperature fluctuations evaluated at $\mathbf{r}$ and $\Phi\left(\mathbf{r}\right)$ is the gravitational potential. Using $$\frac{\Delta\mathrm{T}}{\mathrm{T}}\left(\mathbf{r}\right)=\frac{1}{3\left(2\pi\right)^{3}}\int\mathrm{d^{3}}\mathbf{k}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right),$$
let us suppose now that our system is embedded in a box, with volume $V$ and side $L$, and it has periodicity
$$\Phi\left(\mathbf{r}\right)=\Phi\left(\mathbf{r}+\mathbf{L}\right)=\Phi\left(\mathbf{r}+L\hat{\mathbf{x}}\right)=\Phi\left(\mathbf{r}+L\hat{\mathbf{y}}\right)=\Phi\left(\mathbf{r}+L\hat{\mathbf{z}}\right).$$
Plugging this in the Fourier transform and comparing terms:
$$e^{i\mathbf{k}\cdot\mathbf{r}}=e^{i\mathbf{k}\cdot\left(\mathbf{r}+L\hat{\mathbf{x}}\right)}\Rightarrow e^{iL_{x}k_{x}}=1,$$
and doing the same for $k_{y}$ and $k_{z}$, we have that
$$\mathbf{k}=\frac{2\pi}{L}\left(\mathbf{n}_{x}+\mathbf{n}_{y}+\mathbf{n}_{z}\right).$$
We will use a discretization of the Fourier space. Let us rewrite the Fourier transform as
$$\frac{\Delta\mathrm{T}}{\mathrm{T}}\left(\mathbf{r}\right)=\frac{1}{3V}\sum_{\mathbf{k}}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right)\Rightarrow\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right)=\sum_{\mathbf{k}}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right).$$
For a homogeneous distribution,
$$\left\langle \Phi\left(\mathbf{k}\right)\Phi^{*}\left(\mathbf{q}\right)\right\rangle =P\left(\mathbf{k}\right)\delta\left(\mathbf{k}-\mathbf{q}\right),$$
however, for a Gaussian variables $\phi\left(\mathbf{k}\right)$
$$\left\langle \phi\left(\mathbf{k}\right)\phi^{*}\left(\mathbf{q}\right)\right\rangle =\delta\left(\mathbf{k}-\mathbf{q}\right).$$
Let us assume that
$$\phi\left(\mathbf{k}\right)\equiv\frac{\Phi\left(\mathbf{k}\right)}{\sqrt{P\left(\mathbf{k}\right)}}\Rightarrow\Phi\left(\mathbf{k}\right)=\phi\left(\mathbf{k}\right)\sqrt{P\left(\mathbf{k}\right)}.$$
The function $\Delta\mathrm{T}/\mathrm{T}$ is a real function, for our sum be a real function as well we should impose that
$$\Phi\left(\mathbf{k}\right)=\Phi^{*}\left(-\mathbf{k}\right).$$
Let us break our sum in two hemispheres using the parity condition above using $\left(\theta_{k},\phi_{k}\right)\rightarrow\left(\pi-\theta_{k},\pi+\phi_{k}\right)$ as
$$\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right)	=	\sum_{\mathbf{k}}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right)+\sum_{-\mathbf{k}}e^{-i\mathbf{k}\cdot\mathbf{r}}\Phi\left(-\mathbf{k}\right)=\sum_{\mathbf{k}}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right)+\sum_{-\mathbf{k}}e^{-i\mathbf{k}\cdot\mathbf{r}}\Phi^{*}\left(\mathbf{k}\right)
	=	\sum_{\mathbf{k}}e^{i\mathbf{k}\cdot\mathbf{r}}\Phi\left(\mathbf{k}\right)+e^{-i\mathbf{k}\cdot\mathbf{r}}\Phi^{*}\left(\mathbf{k}\right).$$
Using the redefinition of the Power Spectrum, and specifically for this case, $P\left(\mathbf{k}\right)=P\left(k\right)$ leading that $\sqrt{P\left(\mathbf{k}\right)}$ is already a real function, however
$$\phi\left(\mathbf{k}\right)=\phi^{R}\left(\mathbf{k}\right)+i\phi^{I}\left(\mathbf{k}\right).$$
Plugging the results above in the discretized Fourier transform:
$$\begin{align}
\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right) & =\sum_{\mathbf{k}}\sqrt{P\left(\mathbf{k}\right)}\left\{ e^{i\mathbf{k}\cdot\mathbf{r}}\left[\phi^{R}\left(\mathbf{k}\right)+i\phi^{I}\left(\mathbf{k}\right)\right]+e^{-i\mathbf{k}\cdot\mathbf{r}}\left[\phi^{R}\left(\mathbf{k}\right)-i\phi^{I}\left(\mathbf{k}\right)\right]\right\} ,
\end{align}$$
leading to
$$\boxed{\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right)=2\sum_{\mathbf{k}}\sqrt{P\left(\mathbf{k}\right)}\left[\cos\left(\mathbf{k}\cdot\mathbf{r}\right)\phi^{R}\left(\mathbf{k}\right)-\sin\left(\mathbf{k}\cdot\mathbf{r}\right)\phi^{I}\left(\mathbf{k}\right)\right].}$$
Using $\bf{k}$, and defining
$$\mathbf{n}\equiv\mathbf{n}_{x}+\mathbf{n}_{y}+\mathbf{n}_{z},$$
where $\mathbf{n}\in\mathbb{Z}$, and $$\mathbf{k}=\frac{2\pi}{L}\mathbf{n}\Rightarrow k=\frac{2\pi}{L}n.$$
Let $$\mathbf{r}=R\left(\sin\theta\cos\phi\hat{\mathbf{x}}+\sin\theta\sin\phi\hat{\mathbf{y}}+\cos\theta\hat{\mathbf{z}}\right),$$ and $\mathbf{n}=n\left(\sin\theta_{n}\cos\phi_{n}\hat{\mathbf{x}}+\sin\theta_{n}\sin\phi_{n}\hat{\mathbf{y}}+\cos\theta_{n}\hat{\mathbf{z}}\right)$, $$\mathbf{k}\cdot\mathbf{r}=2\pi n\frac{R}{L}\cos\gamma,$$
where
$$\cos\gamma=\cos\theta\cos\theta_{n}+\sin\theta\sin\theta_{n}\cos\left(\phi-\phi_{n}\right),$$
and
$$\theta_{n}	=\arccos\left(\frac{n_{z}}{n}\right),\qquad\phi_{n}=\arctan\left(\frac{n_{y}}{n_{x}}\right).$$

We can rewrite the equation for the variation of the temperature as
$$\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right)=2\sum_{\mathbf{k}}\sqrt{P\left(\mathbf{k}\right)}\left[\cos\left(2\pi n\frac{R}{L}\cos\gamma\right)\phi^{R}\left(\mathbf{k}\right)-\sin\left(2\pi n\frac{R}{L}\cos\gamma\right)\phi^{I}\left(\mathbf{k}\right)\right].$$
To speed up our numerical code, we can evaluate all possible norms in terms of the first octant, and we do not take care of the terms of the south hemisphere below $k_{z}$ since we'd used the parity relation. Also, we can map all the other points using only the first octant:

![Octants](http://mathworld.wolfram.com/images/eps-gif/Octant_800.gif "Wolfram MathWorld")

| Octant | Relationship | $$\qquad\qquad\qquad\mathrm{\cos}\gamma\qquad\qquad\qquad\qquad\qquad\qquad$$ |
|:-:|:-:|:-:|
| I | $\left(\theta,\phi\right)$ | $\cos\theta\cos\theta_{k}+\sin\theta\sin\theta_{k}\cos\left(\phi-\phi_{k}\right)$ |
| II | $\left(\theta,\pi-\phi\right)$ | $\cos\theta\cos\theta_{k}-\sin\theta\sin\theta_{k}\cos\left(\phi+\phi_{k}\right)$ |
| III | $\left(\theta,\pi+\phi\right)$ | $\cos\theta\cos\theta_{k}-\sin\theta\sin\theta_{k}\cos\left(\phi-\phi_{k}\right)$ |
| IV | $\left(\theta,-\phi\right)$ | $\cos\theta\cos\theta_{k}+\sin\theta\sin\theta_{k}\cos\left(\phi+\phi_{k}\right)$ |

This implies that
$$\begin{eqnarray}
\frac{\Delta\tilde{\mathrm{T}}}{\tilde{\mathrm{T}}}\left(\mathbf{r}\right) & = & 2\sum_{\mathbf{k}}\sqrt{P\left(\mathbf{k}\right)}\left\{ \phi_{1}^{R}\left(\mathbf{k}\right)\cos\left(2\pi n\frac{R}{L}\cos\gamma_{I}\right)-\phi_{1}^{I}\left(\mathbf{k}\right)\sin\left(2\pi n\frac{R}{L}\cos\gamma_{I}\right)\right.\nonumber \\
 &  & +\phi_{2}^{R}\left(\mathbf{k}\right)\cos\left(2\pi n\frac{R}{L}\cos\gamma_{II}\right)-\phi_{2}^{I}\left(\mathbf{k}\right)\sin\left(2\pi n\frac{R}{L}\cos\gamma_{II}\right)\nonumber \\
 &  & +\phi_{3}^{R}\left(\mathbf{k}\right)\cos\left(2\pi n\frac{R}{L}\cos\gamma_{III}\right)-\phi_{3}^{I}\left(\mathbf{k}\right)\sin\left(2\pi n\frac{R}{L}\cos\gamma_{III}\right)\nonumber \\
 &  & \left.+\phi_{4}^{R}\left(\mathbf{k}\right)\cos\left(2\pi n\frac{R}{L}\cos\gamma_{IV}\right)-\phi_{4}^{I}\left(\mathbf{k}\right)\sin\left(2\pi n\frac{R}{L}\cos\gamma_{IV}\right)\right\} .
\end{eqnarray}.$$

We will start first with the invariant scale Harrison-Zel'dovich power spectrum:
$$\mathcal{P}(k)=Ak^{-3},$$
and later, for other geometries. 


```python
__author__ = "Renan Alves de Oliveira, Thiago dos Santos Pereira"
__license__ = "PSFL"
__version__ = "1.0"
__maintainer__ = "Renan Alves de Oliveira"
__email__ = "fisica.renan@gmail.com"

##################
# Basic packages #
##################

import healpy as hp
from numba import jit
import numpy as np
import warnings

get_ipython().magic('matplotlib inline')
warnings.filterwarnings('ignore')

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


@jit
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
```

This part is tricky in terms of speed of my code. I created a function that generate all random numbers inside of some function, and returns an array containing several maps for each realization. However, it is better create an 3-dimentional array containing all random numbers for all realizations and call it inside of the program. This only happens because we are handling with a huge array. The result will be the same if I could implement the generating random array for all possible norms inside a function that generates maps.


```python
###############################################
# Parameters: Number of realizations and kmax #
###############################################

cutoff = int(input("Enter the value for nmax: "))
NaA = unique_rows(listanorms(cutoff)) #Norms and Angles.
NoR = int(input("Enter the number of realizations: ")) #NoR: Number of Realizations.
RA = np.array(np.random.randn(len(NaA), NoR, 8)) #Random Array for all realizations.
```

The following part includes the temperature fluctuations that are computated in pixels. After calculate all of that, you can compute the $\mathrm{a}_{\ell m}$s.


```python
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

@jit(nogil = True)
def thetahp(nside, index):
    """
    Function that returns Polar angle (theta) for Healpy (hp).

    Args:
        nside (int): NSIDE for Healpy.
        index (int): Recpective index for each pixel.

    Returns:
        Float number. 
    """
    return float(hp.pixelfunc.pix2ang(nside,index)[0])

@jit(nogil = True)
def phihp(nside, index):
    """
    Function that returns Azimuthal angle (phi) for Healpy (hp).

    Args:
        nside (int): NSIDE for Healpy.
        index (int): Recpective index for each pixel.

    Returns:
        Float number. 
    """
    return float(hp.pixelfunc.pix2ang(nside,index)[1])

@jit(nogil = True)
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
    out_var = [0.]
    for i in range(0, NoR):
        var = [0.]
        for index in range(hp.pixelfunc.nside2npix(nside)):
            var.append(T(L, A, thetahp(nside,index), phihp(nside,index), i))
        var = np.array(var[1:])
        out_var.append(var)
    return np.array(out_var[1:])

@jit(nogil = True)
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
    out_var = [0.]
    for i in range(0, NoR):
        var = [0.]
        for index in range(hp.pixelfunc.nside2npix(nside)):
            var.append(T(L, A, thetahp(nside,index), phihp(nside,index), i))
        var = np.array(var[1:])
        out_var.append(hp.sphtfunc.map2alm(var))
    return np.array(out_var[1:])
```

The functions $\mathrm{Pixels}$ and $\mathrm{Alms}$ generate arrays. For each row, it gives all pixels or alms for each realization. If you want Healpy shows maps for each specific realization, you can refer for example:

- hp.mollview(Pixels(16, 1., 1.)[0]) # first realization...
- hp.mollview(Pixels(16, 1., 1.)[1]) # second realization...
