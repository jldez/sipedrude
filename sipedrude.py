import numpy as np
import matplotlib.pyplot as plt



def R(epsilon):
    return (epsilon-1)/(epsilon+1)

def F(s):
    return (s**2+1)**0.5-s

def G(s):
    return ((s**2+4)**0.5+s)/2 - (s**2+1)**0.5


def gamma_t(epsilon, f, s):
    return 2*(epsilon-1)/(4*np.pi*(2+(1-f)*(epsilon-1)*(F(s)-R(epsilon)*G(s))))

def gamma_z(epsilon, f, s):
    return (epsilon-1)/(4*np.pi*(epsilon-(1-f)*(epsilon-1)*(F(s)+R(epsilon)*G(s))))


def ts(epsilon, theta):
    return 2*abs(np.cos(theta))/(abs(np.cos(theta)) + (epsilon-(np.sin(theta))**2)**0.5)

def tx(epsilon, theta):
    return 2*(epsilon-(np.sin(theta))**2)**0.5/(epsilon*abs(np.cos(theta)) + (epsilon-(np.sin(theta))**2)**0.5)

def tz(epsilon, theta):
    return 2*np.sin(theta)/(epsilon*abs(np.cos(theta)) + (epsilon-(np.sin(theta))**2)**0.5)


def hss(kappa, epsilon):
    return 2*complex(0,1)/((1-kappa**2)**0.5+(epsilon-kappa**2)**0.5)

def hkk(kappa, epsilon):
    return 2*complex(0,1)*((epsilon-kappa**2)*(1-kappa**2))**0.5/(epsilon*(1-kappa**2)**0.5+(epsilon-kappa**2)**0.5)

def hkz(kappa, epsilon):
    return 2*complex(0,1)*kappa*(epsilon-kappa**2)**0.5/(epsilon*(1-kappa**2)**0.5+(epsilon-kappa**2)**0.5)

def hzk(kappa, epsilon):
    return 2*complex(0,1)*kappa*(1-kappa**2)**0.5/(epsilon*(1-kappa**2)**0.5+(epsilon-kappa**2)**0.5)

def hzz(kappa, epsilon):
    return 2*complex(0,1)*kappa**2/(epsilon*(1-kappa**2)**0.5+(epsilon-kappa**2)**0.5)


def kappa_plus(kappa_x, kappa_y, theta):
    return (kappa_x**2+(np.sin(theta)+kappa_y)**2)**0.5

def kappa_minus(kappa_x, kappa_y, theta):
    return (kappa_x**2+(np.sin(theta)-kappa_y)**2)**0.5


# def v(kappa, epsilon, )


kappa_max = 6
N_kappa = 100

kappa_x, kappa_y = np.meshgrid(np.linspace(-kappa_max,kappa_max,N_kappa), np.linspace(-kappa_max,kappa_max,N_kappa))
kappa_plus(kappa_x, kappa_y, 0)