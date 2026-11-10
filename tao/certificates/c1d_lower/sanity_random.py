import numpy as np
rng=np.random.default_rng(1)
worst=0
for trial in range(200):
    n=rng.integers(2,40); v=rng.normal(size=n)+rng.uniform(0,2); v/=v.sum()*(0.5/n)   # int f = 1, cells of width 1/(2n)
    hcell=0.5/n
    # g = f*f is piecewise linear with knots at -1/2 + k*hcell, values hcell*(v*v)[k-1]
    c=np.convolve(v,v)*hcell; knots=np.concatenate([[0],c,[0]]); xk=-0.5+hcell*np.arange(len(knots))
    M=np.abs(knots).max()
    xs=np.linspace(-0.5,0.5,40001); g=np.interp(xs,xk,knots); h=M-g
    assert abs(np.trapezoid(g,xs)-1)<1e-4*M, (np.trapezoid(g,xs),M)
    for z in [1+0.3j,1.2+0.34j,0.9-0.2j,1.1+0.0j,1.23-0.34j]:
        E=np.trapezoid(h*np.exp(-2j*np.pi*z*xs),xs); y=abs(z.imag)
        B=M-1 if y==0 else M/(np.pi*y)*(np.exp(np.pi*y)-np.exp(np.pi*y/M))
        worst=max(worst,abs(E)/B)
        # also the Rouche consequence: |M sinc - fhat^2| = |E|
print("max |h^(z)| / Bnd over random tests:",worst, "(must be <= 1)")
