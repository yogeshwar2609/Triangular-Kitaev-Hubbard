from HamiltonianPy import *
from numpy.linalg import inv
import numpy as np

__all__=['nameb','namef','nnb','parametermap','idfmap','t','h','U','V1','V2','V3','T1']

# The configs of the model
nameb='BTKH'
namef='FTKH'
nnb=3

# parametermap
parametermap=None

# idfmap
idfmap=lambda pid: Fock(norbital=1,nspin=2,nnambu=1)

# kitaev hopping
def kitaev(bond):
    theta=azimuthd(bond.rcoord)
    if abs(theta)<RZERO or abs(theta-180)<RZERO: return sigmax('sp')
    if abs(theta-60)<RZERO or abs(theta-240)<RZERO: return sigmay('sp')
    if abs(theta-120)<RZERO or abs(theta-300)<RZERO: return sigmaz('sp')

# magnetic field
def magnetic(bond,theta,phi):
    theta,phi=theta*np.pi,phi*np.pi*2
    return sigmax('sp')*np.sin(theta)*np.cos(phi)+sigmay('sp')*np.sin(theta)*np.sin(phi)+sigmaz('sp')*np.cos(theta)

# terms
t=lambda statistics,**parameters: Hopping('t',parameters['t'],indexpacks=kitaev,statistics=statistics)
h=lambda statistics,**parameters: Onsite('h',parameters['h'],indexpacks=lambda bond: magnetic(bond,parameters['theta'],parameters['phi']),statistics=statistics)
U=lambda statistics,**parameters: Hubbard('U',parameters['U'],statistics=statistics,modulate=True)
V1=lambda statistics,**parameters: Coulomb('V1',parameters['V1'],neighbour=1,statistics=statistics,modulate=True)
V2=lambda statistics,**parameters: Coulomb('V2',parameters['V2'],neighbour=2,statistics=statistics,modulate=True)
V3=lambda statistics,**parameters: Coulomb('V3',parameters['V3'],neighbour=3,statistics=statistics,modulate=True)

# cluster
T1=Triangle('T1')
