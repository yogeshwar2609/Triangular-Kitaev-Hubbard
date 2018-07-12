import mkl
import numpy as np
from HamiltonianPy import *
from source import *
from collections import OrderedDict

def tbatasks(name,parameters,lattice,terms,boundaries,statistics,jobs=()):
    import HamiltonianPy.FreeSystem as TBA
    tba=tbaconstruct(name,parameters,lattice,terms,boundaries,statistics)
    #tba.generator.view(suspend=True)
    if 'EB' in jobs:
        path=hexagon_gkm(reciprocals=lattice.reciprocals,nk=100)
        tba.register(EB(name='EB',path=path,run=TBA.TBAEB))
        tba.summary()
    if 'GSE' in jobs:
        tba.register(TBA.GSE(filling=0.1,run=TBA.TBAGSE))
        tba.summary()

def edtasks(name,parameters,basis,lattice,terms,boundaries,statistics,jobs=()):
    import HamiltonianPy.ED as ED
    ed=edconstruct(name,parameters,[basis],lattice,terms,boundaries,statistics)
    if 'EL1' in jobs: ed.register(ED.EL(name='EL1',path=BaseSpace(('theta1',np.linspace(0,1.0,21))),ns=15,run=ED.EDEL))
    if 'EL2' in jobs: ed.register(ED.EL(name='EL2',path=BaseSpace(('theta2',np.linspace(0,1.0,21))),ns=10,run=ED.EDEL))
    if 'ELU' in jobs: ed.register(ED.EL(name='EL2',path=BaseSpace(('U',np.linspace(0,20.0,21))),ns=1,nder=2,run=ED.EDEL))
    if 'EIGS' in jobs: ed.register(ED.EIGS(name='EIGS',ne=10,run=ED.EDEIGS))
    ed.summary()

if __name__=='__main__':
    mkl.set_num_threads(1)
    Engine.DEBUG=True
    Engine.MKDIR=False

    # parameters
    parameters=OrderedDict()
    parameters['t']=1.0

    # tba
    #tbatasks(namef,parameters,T1('1P-1P',nnb),[t],None,'f',jobs=['EB'])

    # ed
    parameters['U']=20.0
    parameters['V']=5.0
    parameters['theta1']=0.0
    parameters['theta2']=0.0

    m,n,filling=4,5,10
    lattice=T1('%sP-%sP'%(m,n),nnb)
    boundaries=[('theta1','theta2'),(0.0,0.0),lattice.vectors,twisttransformation]

    # fermi
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'f',jobs=['EL1'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'f',jobs=['EL2'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'f',jobs=['ELU'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'f',jobs=['EIGS'])

    # bose
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'b',jobs=['EL1'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'b',jobs=['EL2'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'b',jobs=['ELU'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V],boundaries,'b',jobs=['EIGS'])
