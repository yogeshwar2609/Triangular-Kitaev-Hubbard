from HamiltonianPy import *
from numpy.linalg import inv
import numpy as np

__all__=['nameb','namef','nnb','parametermap','idfmap','t','U','V','T1']

# The configs of the model
nameb='BTKH'
namef='FTKH'
nnb=1

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

# terms
t=lambda statistics,**parameters: Hopping('t',parameters['t'],indexpacks=kitaev,statistics=statistics)
U=lambda statistics,**parameters: Hubbard('U',parameters['U'],statistics=statistics,modulate=True)
V=lambda statistics,**parameters: Coulomb('V',parameters['V'],neighbour=1,statistics=statistics,modulate=True)

# cluster
T1=Triangle('T1')
