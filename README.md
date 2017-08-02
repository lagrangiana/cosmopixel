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


![Octants](http://mathworld.wolfram.com/images/eps-gif/Octant_800.gif)

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
and later, for Bianchi's $I$ and $VII_{0}$.
