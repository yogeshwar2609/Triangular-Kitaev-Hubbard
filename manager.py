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
        if len(lattice.vectors)==2:
            path=hexagon_gkm(reciprocals=lattice.reciprocals,nk=100)
            #path=KPath(KMap(lattice.reciprocals,'H:G-K,K-M,M-G'),nk=100)
        elif len(lattice.vectors)==1:
            path=line_bz(reciprocals=lattice.reciprocals,nk=400)
        else:
            path=None
        tba.register(EB(name='EB',path=path,run=TBA.TBAEB))
        tba.summary()
    if 'GR' in jobs:
        path=KSpace(reciprocals=lattice.reciprocals,nk=50)
        tba.register(EB(name='GR',path=path,plot=False,savedata=False,run=TBA.TBAEB))
        gap=tba.records['GR'][:,2].min()-tba.records['GR'][:,1].max()
        bandwidth=tba.records['GR'][:,1].max()-tba.records['GR'][:,1].min()
        tba.log<<'gap,ratio: %s, %s\n'%(gap,gap/bandwidth)
        tba.summary()

def edtasks(name,parameters,basis,lattice,terms,boundaries,statistics,jobs=()):
    import HamiltonianPy.ED as ED
    ed=edconstruct(name,parameters,[basis],lattice,terms,boundaries,statistics)
    if 'EL1' in jobs: ed.register(ED.EL(name='EL1',path=BaseSpace(('theta1',np.linspace(0,1.0,21))),ns=15,run=ED.EDEL))
    if 'EL2' in jobs: ed.register(ED.EL(name='EL2',path=BaseSpace(('theta2',np.linspace(0,1.0,21))),ns=15,run=ED.EDEL))
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
    #tbatasks(namef,parameters,T1('1P-1P',nnb),[t,h],None,'f',jobs=['EB'])
    #tbatasks(namef,parameters,T1('1P-20O',nnb),[t,h],None,'f',jobs=['EB'])
    #tbatasks(namef,parameters,T1('1P-1P',nnb),[t,h],None,'f',jobs=['GAP'])

    # ed
    parameters['U']=120.0
    parameters['V1']=40.0
    parameters['V2']=40.0
    parameters['V3']=40.0
    parameters['theta1']=0.0
    parameters['theta2']=0.0

    m,n,filling=3,5,15
    lattice=T1('%sP-%sP'%(m,n),nnb)
    boundaries=[('theta1','theta2'),(0.0,0.0),lattice.vectors,twisttransformation]

    # fermi
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'f',jobs=['EL1'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2,V3],boundaries,'f',jobs=['EL2'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'f',jobs=['ELU'])
    #edtasks(namef,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'f',jobs=['EIGS'])

    # bose
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'b',jobs=['EL1'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'b',jobs=['EL2'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'b',jobs=['ELU'])
    #edtasks(nameb,parameters,FBasis(m*n*2,m*n*2/filling),lattice,[t,U,V1,V2],boundaries,'b',jobs=['EIGS'])
