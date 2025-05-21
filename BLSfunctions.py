def spherical_hn1(n,x):
    # compute value of spherical Hankel function of the first kind of order n at x. Related to Hankel function by sqrt(pi/2x) prefactor and shift of order by + 0.5.
    import numpy as np
    import scipy as sp #necessary packages
    return np.sqrt(np.pi/(2*x))*sp.special.hankel1(n+0.5,x) #compute values(s)

def spherical_hn2(n,x):
    # compute value of spherical Hankel function of the second kind of order n at x. Related to Hankel function by sqrt(pi/2x) prefactor and shift of order by + 0.5.
    import numpy as np
    import scipy as sp #necessary packages
    return np.sqrt(np.pi/(2*x))*sp.special.hankel2(n+0.5,x) #compute values(s)

def spherical_hn1d(n,x):
    #compute value of derivative of spherical Hankel function of first kind. Recurrence relation used comes from https://www.wolframalpha.com/input?i=derivative+of+spherical+Hankel+function.
    return (-1/(2*x))*(spherical_hn1(n,x) + x*(spherical_hn1(n+1,x)-spherical_hn1(n-1,x)))

def spherical_hn2d(n,x):
    #compute value of derivative of spherical Hankel function of second kind. Recurrence relation used comes from https://www.wolframalpha.com/input?i=derivative+of+spherical+Hankel+function+of+second+kind
    return (-1/(2*x))*(spherical_hn2(n,x) + x*(spherical_hn2(n+1,x)-spherical_hn2(n-1,x)))

def RBjn(n,x):
    #compute value of Riccati-Bessel function of first kind with order n at x. Related to spherical Bessel function of first kind through multiplication by argument. 
    import scipy as sp
    return x*sp.special.spherical_jn(n,x)

def RByn(n,x):
    #compute value of Riccati-Bessel function of second kind with order n at x. Related to spherical Bessel function of second kind through multiplication by argument. 
    import scipy as sp
    return x*sp.special.spherical_yn(n,x)

def RBhn1(n,x):
    #compute value of Riccati-Bessel function of third kind with order n at x. Related to spherical Bessel function of third kind (spherical Hankel function of first kind) through multiplication by argument. 
    return x*spherical_hn1(n,x)

def RBhn2(n,x):
    #compute value of Riccati-Bessel function of fourth kind with order n at x. Related to spherical Bessel function of fourth kind (spherical Hankel function of seocnd kind) through multiplication by argument. 
    return x*spherical_hn2(n,x)

def RBjnd(n,x):
    #compute value of derivative of Riccati-Bessel function of first kind with order n at x. 
    import scipy as sp
    return x*sp.special.spherical_jn(n,x,True) + sp.special.spherical_jn(n,x)

def RBynd(n,x):
    #compute value of derivative of Riccati-Bessel function of second kind with order n at x. 
    import scipy as sp
    return x*sp.special.spherical_yn(n,x,True) + sp.special.spherical_yn(n,x)

def RBhn1d(n,x):
    #compute value of derivative of Riccati-Bessel function of third kind with order n at x. 
    return (1/2)*(spherical_hn1(n,x) - x*(spherical_hn1(n+1,x) - spherical_hn1(n-1,x)))

def RBhn2d(n,x):
    #compute value of derivative of Riccati-Bessel function of fourth kind with order n at x. 
    return (1/2)*(spherical_hn2(n,x) - x*(spherical_hn2(n+1,x) - spherical_hn2(n-1,x)))

def Mie_an(n, m, x):
    #compute Mie a coefficient of order n for a particle of relative refractive index m and size parameter x.
    return ((m*RBjn(n,m*x)*RBjnd(n,x))-(RBjn(n,x)*RBjnd(n,m*x)))/((m*RBjn(n,m*x)*RBhn1d(n,x))-(RBhn1(n,x)*RBjnd(n,m*x)))

def Mie_bn(n, m, x):
    #compute Mie b coefficient of order n for a particle of relative refractive index m and size parameter x.
    return ((RBjn(n,m*x)*RBjnd(n,x))-(m*RBjn(n,x)*RBjnd(n,m*x)))/((RBjn(n,m*x)*RBhn1d(n,x))-(m*RBhn1(n,x)*RBjnd(n,m*x)))

def Mie_pin(n, theta):
    #angular function pi of order 1-n evaluated at theta. Defined and calculated as per Bohren and Huffman P 94-95. n must be >= 1.
    import numpy as np
    if n == 1:
        return np.array([[np.ones_like(theta)]]) #pi_1 = 1 by definition
    else: 
        pin = np.zeros((n+1, len(theta))) #space for pi_n values
        pin[0, :] = np.zeros_like(theta) #row for pi_0
        pin[1, :] = np.ones_like(theta) #row for pi_1
        for i1 in range(2, n+1):
            pin[i1,:] = (((2*i1)-1)/(i1-1))*np.cos(theta)*pin[i1-1,:] - (i1/(i1-1))*pin[i1-2,:] #calculate higher order from recurrence relation
    return pin[1:,:] #order 0 not relevant for Mie theory so don't output

def Mie_taun(n, theta):
    #angular function tau of order 1-n evaluated at theta. Defined and calculated as per Bohren and Huffman P 94-95. n must be >= 1.
    import numpy as np
    pin = Mie_pin(n, theta) #first calculate pi_n for each n at each theta
    taun = np.zeros_like(pin) #space for values of tau_n
    for i1 in range(n):
        if i1 == 0:
            taun[i1,:] = np.cos(theta)*pin[i1,:]
        else:
            taun[i1,:] = (i1+1)*np.cos(theta)*pin[i1,:] - (i1+2)*pin[i1-1,:] #compute tau_n for each value of n (indices are a little confusing)
    return taun