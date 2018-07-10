from HamiltonianPy import *
from numpy.linalg import norm
import numpy as np

__all__=['nameb','namef','nnb','parametermap','idfmap','t','tw','U','T1']

# The configs of the model
nameb='BTKH'
namef='FTKH'
nnb=1

# parametermap
def parametermap(parameters):
    data={}
    data['t']=parameters['t']
    data['txp']=parameters['t']*np.exp(+1.0j*parameters['theta1']*np.pi)
    data['txm']=parameters['t']*np.exp(-1.0j*parameters['theta1']*np.pi)
    data['typ']=parameters['t']*np.exp(+1.0j*parameters['theta2']*np.pi)
    data['tym']=parameters['t']*np.exp(-1.0j*parameters['theta2']*np.pi)
    data['U']=parameters['U']
    return data

# idfmap
idfmap=lambda pid: Fock(norbital=1,nspin=2,nnambu=1)

# kitaev hopping
def kitaev(bond):
    theta=azimuthd(bond.rcoord)
    if abs(theta)<RZERO or abs(theta-180)<RZERO: return sigmax('sp')
    if abs(theta-60)<RZERO or abs(theta-240)<RZERO: return sigmay('sp')
    if abs(theta-120)<RZERO or abs(theta-300)<RZERO: return sigmaz('sp')

# intra cluster amplitude
intracell=lambda bond: 1.0 if bond.isintracell() else 0.0

# twisted boundary amplitudes
txp=lambda bond: 1.0 if bond.icoord[0]>0 else 0.0
txm=lambda bond: 1.0 if bond.icoord[0]<0 else 0.0
typ=lambda bond: 1.0 if bond.icoord[1]>0 else 0.0
tym=lambda bond: 1.0 if bond.icoord[1]<0 else 0.0

# terms
t=lambda statistics,**parameters: Hopping('t',parameters['t'],amplitude=intracell,indexpacks=kitaev,statistics=statistics)
tw=[    lambda statistics,**parameters: Hopping('txp',parameters['t']*np.exp(+1.0j*parameters['theta1']*np.pi),amplitude=txp,indexpacks=kitaev,statistics=statistics,modulate=True),
        lambda statistics,**parameters: Hopping('txm',parameters['t']*np.exp(-1.0j*parameters['theta1']*np.pi),amplitude=txm,indexpacks=kitaev,statistics=statistics,modulate=True),
        lambda statistics,**parameters: Hopping('typ',parameters['t']*np.exp(+1.0j*parameters['theta2']*np.pi),amplitude=typ,indexpacks=kitaev,statistics=statistics,modulate=True),
        lambda statistics,**parameters: Hopping('tym',parameters['t']*np.exp(-1.0j*parameters['theta2']*np.pi),amplitude=tym,indexpacks=kitaev,statistics=statistics,modulate=True)
        ]
U=lambda statistics,**parameters: Hubbard('U',parameters['U'],statistics=statistics,modulate=True)

# cluster
T1=Triangle('T1')
