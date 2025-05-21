### Import Packages
import numpy as np
import BLSfunctions as BLS
import time
### Set Input Parameters
wl = np.linspace(0.4, 0.8, 1201) #wavelength range and resolution (in microns)
r = np.linspace(1.0,3.0,21) #radii to use in calculations (in microns)
n = 1.4 #real part of refractive index
k = 0.0 #imaginary part of refractive index
m = complex(n, k) #complex refractive index 
theta_i = np.pi/2 #incident light theta
phi_i = np.pi/2 #incident light phi
NA = 0.5 #numerical aperture of collection lens
### Calculate Spectra
wn = (2*np.pi)/wl #wavenumber of incident light
vi = np.transpose(np.array([np.sin(theta_i)*np.cos(phi_i), np.sin(theta_i)*np.sin(phi_i), np.cos(theta_i)])) #incident wavevector
x = np.zeros((np.size(r), len(wl))) #space for size parameters
for i0, lam in enumerate(wl):
    x[:,i0] = (2*np.pi*r)/lam #calculate size parameters for each radius
omax = np.round(x + (4*(x**(1/3))) + 2).astype(int) #maximum order to which sums need to be evaluated for each size parameter
o = np.linspace(1, np.max(omax), np.max(omax)) #list of orders used to calculate prefactors, coefficients and angular functions
omaxr = np.max(omax, axis=1) #maximum order required in summation for each radius
npts = np.round((0.75*NA*omaxr)) + 18 #number of pts to use for numerical integration
odd = np.mod(npts, 2) == 1 #find odd grid points
npts[odd] += 1 #make odd grid points even
pf = ((2*o)+1)/(o*(o+1)) #prefactors for eventual sum
spectra = np.zeros_like(x) #space for scattering spectra
for i0 in range(np.shape(x)[0]):
    lstart = time.time()
    an = np.zeros((len(o), len(wl)), dtype=complex)
    bn = np.zeros_like(an) #space for coefficients for given radius
    for i1 in range(omaxr[i0]):
        an[i1,:] = pf[i1]*BLS.Mie_an(int(o[i1]), m, x[i0,:])
        bn[i1,:] = pf[i1]*BLS.Mie_bn(int(o[i1]), m, x[i0,:]) #calculate given order of coefficients multiplied by prefactor
    theta_s = np.linspace(0, np.arcsin(NA), int(npts[i0]))#range of scattered theta to look at (angle between scattered wavevector and detector axis)
    phi_s = np.linspace(0, 2*np.pi, int(npts[i0])) #range of scattered phi to use (angle between scattered wavevector and x-axis)
    theta_c = np.zeros((len(theta_s), len(phi_s)))
    for i1 in range(len(theta_s)):
        for i2 in range(len(phi_s)):
            vs = np.transpose(np.array([np.sin(theta_s[i1])*np.cos(phi_s[i2]), np.sin(theta_s[i1])*np.sin(phi_s[i2]), np.cos(theta_s[i1])])) #scattering vector components (minus k, which is not needed for collection angle calculation)
            theta_c[i1,i2] = np.arccos(np.dot(vi,vs)) #collection angle for scattered wavevector
    pin = np.zeros((len(theta_s), len(phi_s), omaxr[i0]))
    taun = np.zeros_like(pin) #space for values of pi_n and tau_n (independent of size parameter except for summation maximum, so can be evaluated before loop and called on later)
    for i1 in range(len(theta_s)):
        pn = BLS.Mie_pin(omaxr[i0], theta_c[i1,:])
        tn = BLS.Mie_taun(omaxr[i0], theta_c[i1,:]) #calculate necessary range of pi_n and tau_n across collection angle range
        pin[i1,:,:] = np.transpose(pn)
        taun[i1,:,:] = np.transpose(tn) #put into matrices
    spectrum = np.zeros_like(wl) #space for spectrum at given radius
    for i1 in range(len(wl)):
        S1 = complex(0,0)
        S2 = complex(0,0) #initial values for summation
        for i2 in range(omax[i0,i1]):
            S1 += (an[i2,i1]*pin[:,:,i2]) + (bn[i2,i1]*taun[:,:,i2])
            S2 += (an[i2,i1]*taun[:,:,i2]) + (bn[i2,i1]*pin[:,:,i2]) #calculate summation terms
        S11 = (1/(2*(wn[i1]**2)))*((np.abs(S2)**2) + (np.abs(S1)**2)) #calculate S11 as function of theta and phi
        theta_int_terms = np.zeros_like(S11)
        for i2 in range(len(phi_s)):
            theta_int_terms[:,i2] = np.sin(theta_s)*S11[:,i2] #calculate terms for theta integration (multiply S11 terms by sin(theta))
        phi_int_terms = ((theta_s[1] - theta_s[0])/2)*(theta_int_terms[0,:] + theta_int_terms[-1,:] + 2*np.einsum('ij->j',theta_int_terms[1:-1,:])) #perform theta integration using Simpson's 1/3 rule for regular spaced points
        spectrum[i1] = ((phi_s[1] - phi_s[0])/2)*(phi_int_terms[0] + phi_int_terms[-1] + 2*np.einsum('i->', phi_int_terms[1:-1])) #perform phi integration using Simpson's 1/3 rule and put into vector
    spectra[i0,:] = spectrum #put spectrum into matrix
    print("Spectrum", str(i0+1), "of", str(np.size(r)), "complete in", time.time()-lstart, "s.")