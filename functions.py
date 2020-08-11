# -*- coding: utf-8 -*-
import numpy as np
    
def declareSTR(eps):
# Strategy: (sig_PGG, sig_0, A_sig0_T, A_sig0_noT, A_sig1_T, A_sig1_noT )
    nelem=6    
    STRm=np.zeros((2**nelem,nelem))
    nST=0
    i=np.zeros((nelem),dtype=int)
    for i[0] in range (1,-1,-1):
        for i[1] in range (1,-1,-1):
            for i[2] in range (1,-1,-1):
                for i[3] in range (1,-1,-1):
                    for i[4] in range (1,-1,-1):
                        for i[5] in range (1,-1,-1):
                            nST+=1
                            for j in range (0,nelem):
                                STRm[nST-1,j]=(1.-eps)*i[j]+eps*(1-i[j])
                                #if (j==0 or j==1): STRm[nST-1,j]=i[j]  # signaling without error
    return STRm

def declareSTR_CD(eps):
# AllC and AllD  
    STRm=np.zeros((2,6))
    eps1=1.-eps
    STRm[0,:]=[0, 0, eps, eps, eps, eps]
    STRm[1,:]=[0, 0, eps1, eps1, eps1, eps1]
    return STRm

def declareSTR_REC(eps):
# AllC, AllD, TFT, ATFT  
    STRm=np.zeros((4,6))
    eps1=1.-eps
    STRm[0,:]=[0, 0, eps, eps, eps, eps]
    STRm[1,:]=[0, 0, eps1, eps, eps1, eps]
    STRm[2,:]=[0, 0, eps, eps1, eps, eps1]
    STRm[3,:]=[0, 0, eps1, eps1, eps1, eps1]
    return STRm

def declareSTR_SIG(eps):
# Only signalling strategies
    nelem=6    
    STRm=np.zeros((16,nelem))
    eps1=1.-eps
    nST=0
    i=np.zeros((nelem),dtype=int)
    for i[0] in range (1,-1,-1):
        for i[1] in range (1,-1,-1):
            for i[2] in range (1,-1,-1):
                    for i[4] in range (1,-1,-1):
                            nST+=1
                            STRm[nST-1,0]=eps1*i[0]+eps*(1.-i[0])
                            STRm[nST-1,1]=eps1*i[1]+eps*(1.-i[1])
                            STRm[nST-1,2]=eps1*i[2]+eps*(1.-i[2])
                            STRm[nST-1,3]=STRm[nST-1,2]
                            STRm[nST-1,4]=eps1*i[4]+eps*(1.-i[4])
                            STRm[nST-1,5]=STRm[nST-1,4]
    return STRm

    
def declareSTATE(N):
# State of acting
# k individuls from the first strategy 
    nelem=6    
    STATEmat=np.zeros((N+1,int((N/2+1)**2)+1,2))-1
    nSTATEv=np.zeros((N+1),dtype=int)
    i=np.zeros((nelem),dtype=int)
    for k in range(0,N+1):
        nSTATE=0
        for i in range (k,-1,-1):
            for j in range (N-k,-1,-1):
                nSTATE+=1
                STATEmat[k,nSTATE-1,:]=[i, j];
        nSTATEv[k]=nSTATE
    return STATEmat, nSTATEv
    

def calcERS(b,c,cs,lamb,N,Z,M,eps):
    H=calcH(N,Z)
    STRm=declareSTR(eps); nSTR=STRm.shape[0];
    STATEmat,nSTATEv=declareSTATE(N)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    listERS=[]
    for i in range(0,nSTR): # resident
        isERS=1
        SIGi=STRm[i,0:2]; ACTi=STRm[i,2:6]; 
        for j in range(0,nSTR): #mutant
          if i!=j:
            SIGj=STRm[j,0:2]; ACTj=STRm[j,2:6];
            k=N-1; BCi,BCj=calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
            k=N; BCiR,ttt=calcBC2st(SIGi,ACTi,SIGj,ACTj,N,N,M,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
            PAYi=H[N-1,Z-2]*np.sum(coef*BCiR)+H[N-2,Z-2]*np.sum(coef*BCi)
            PAYj=np.sum(coef*BCj)
            if PAYj>PAYi: isERS=0; break
        if isERS==1: listERS=np.append(listERS,i)
    listERS=listERS.astype(int)  
    return listERS

def calcBC2st1nat(SIGi,ACTi,SIGj,ACTj,k,N,M,Q,STATEm,w): ### MODIFICAR
# Calculate BC, where BC[0]*k*c-BC[1]*c-BC[2]*cs.
#  k individuls type i; N-k individuals type j  
    import scipy.sparse.linalg as lin
    from scipy.stats import binom # import scipy.stats.binom as binom
    N2=1.*N/2.; Nk=N-k;
    nSTATE=(N-k+1)*(k+1); MAT=np.zeros((nSTATE,nSTATE))
    ns=k*SIGi+Nk*SIGj
    consS=np.piecewise(ns,[ns<Q,ns>=Q],[0.,1.]) 
    nc= STATEm[:,0]+STATEm[:,1] # in the current state (that has passed)
    consA=np.piecewise(nc,[nc<M,nc>=M],[0.,1.]) 
    Pcoopi=consA*consS*ACTi[2]+consA*(1.-consS)*ACTi[0]+(1.-consA)*consS*ACTi[3]+(1.-consA)*(1.-consS)*ACTi[1]
    Pcoopj=consA*consS*ACTj[2]+consA*(1.-consS)*ACTj[0]+(1.-consA)*consS*ACTj[3]+(1.-consA)*(1.-consS)*ACTj[1]
    for j in range(0,nSTATE):
        MAT[:,j]=binom.pmf(STATEm[j,0],k,Pcoopi)*binom.pmf(STATEm[j,1],Nk,Pcoopj)   
    if w>=1.:
#        val,vect=lin.eigs(np.transpose(MAT),k=1,which='LR'); vect=np.real(vect/np.sum(vect))  # another method     
        from discreteMarkovChain import markovChain
        mc=markovChain(MAT)
        mc.computePi('eigen') # We can use 'linear', 'power', 'krylov' or 'eigen'
        vect=(mc.pi).reshape(-1,1)
    else:
        vect=(1-w)*np.linalg.inv((np.identity(nSTATE)-w*MAT))[nSTATE-1,:]

    BCi=np.zeros((3)); BCj=np.zeros((3))
    benef=np.dot(nc*consA/N,vect)
    if (k!=0):
        BCi[0]=benef
        BCi[1]=np.dot(STATEm[:,0]/k,vect)
        BCi[2]=SIGi
    if(k!=N):
        BCj[0]=benef
        BCj[1]=np.dot(STATEm[:,1]/Nk,vect)
        BCj[2]=SIGj
    return BCi, BCj

def calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,Q,STATEm,w):  
# output: different rows for different states of nature; columns probabilities of benefit, cooperating, signaling
    BCi1,BCj1=calcBC2st1nat(SIGi[0],ACTi,SIGj[0],ACTj,k,N,M,Q,STATEm,w)
    BCi2,BCj2=calcBC2st1nat(SIGi[1],ACTi,SIGj[1],ACTj,k,N,M,Q,STATEm,w)
    BCi=np.stack((BCi1, BCi2)) ;BCj=np.stack((BCj1, BCj2))
    return BCi, BCj
    
def calcH (N,Z):   
    import numpy as np
    from scipy.stats import hypergeom  
    H=np.zeros((N+1,Z+1))
    H[0,0]=1  # H(0,:)=0, H(0,0)=1  Attention!
    for K in range(1,Z+1):         
        for k in range(0,N+1):
            H[k,K]=hypergeom.pmf(k,Z-1,K,N-1)
    return H
    
def calcFIX1vec(STi,STj,STRm,N,Z,M,Q,STATEmat,nSTATEv,H,w):
# i invades j (j->i)
# output: (Z-1,2,3)    
    SIGi=STRm[STi,0:2]; ACTi=STRm[STi,2:6]; SIGj=STRm[STj,0:2]; ACTj=STRm[STj,2:6];

    BCki=np.zeros((N+1,2,3)); BCkj=np.zeros((N+1,2,3));
    for k in range(0,N+1):
        BCi,BCj=calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        BCki[k,...]=BCi; BCkj[k,...]=BCj;
    PIKi=np.einsum('kK,kij->Kij',H[0:N,0:Z-1],BCki[1:N+1,...])
    PIKj=np.einsum('kK,kij->Kij',H[0:N,1:Z],BCkj[0:N,...])
    PIKjI=np.einsum('kK,kij->Kij',H[0:N,0:Z-1],BCkj[::-1,...][1:N+1,...])
    PIKiI=np.einsum('kK,kij->Kij',H[0:N,1:Z],BCki[::-1,...][0:N,...])
    DIFK=PIKi-PIKj
    DIFKI=PIKjI-PIKiI
    sumterm=np.cumsum(DIFK,axis=0)
    sumtermI=np.cumsum(DIFKI,axis=0)
    return sumterm, sumtermI    

def calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w):
#  from i to j (i->j)
    nSTR=STRm.shape[0]; 
    fixMvec=np.zeros((nSTR,nSTR,Z-1,2,3))
    for i in range(0,nSTR):
        for j in range(i+1,nSTR):
            pfixvec,pfixIvec=calcFIX1vec(i,j,STRm,N,Z,M,Q,STATEmat,nSTATEv,H,w)
            fixMvec[j,i,...]=pfixvec; fixMvec[i,j,...]=pfixIvec
    return fixMvec

    
    
def calcFIXM(coef,expb,Z,fixMvec):
#  from i to j (i->j)
    fixM=1./(1.+np.sum(expb**np.einsum('ijmab,ab->ijm',fixMvec,coef),axis=2))
    np.fill_diagonal(fixM, 0.) 
    nSTR=len(fixM)
  
    fixM=fixM/nSTR
    suma=np.sum(fixM,axis=1)
    fixM[range(nSTR),range(nSTR)]=1.-suma

    return fixM


def calcSD(fixM):   
    from discreteMarkovChain import markovChain
    mc=markovChain(fixM)
    mc.computePi('linear') # We can use 'linear', 'power', 'krylov' or 'eigen'
    SD=(mc.pi).reshape(-1,1)   
    return SD
    
def calcHOMO(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w):
# PAYhomo: (nSx1); payoffs homogeneous populations
# COOPhome: (nSx3); 0:cooperation level, 1: cooperation level lambda=1, 2: cooperation level lambda=0
# COOPtotal: scalar; total cooperation level (taking into account SD)      
    nSTR=SD.shape[0];
    STRm=declareSTR(eps);
    PAYhomo=np.zeros((nSTR)); COOPhomo=np.zeros((nSTR,3)); COOPtot=np.zeros((3))
    for i in range(0,nSTR):
        SIG=STRm[i,0:2]; ACT=STRm[i,2:6]; k=1;
        BCi,BCj=calcBC2st(SIG,ACT,SIG,ACT,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        PAYhomo[i]=np.sum(BCi*coef)
        COOPhomo[i,1]=BCi[0,1]
        COOPhomo[i,2]=BCi[1,1]
    COOPhomo[:,0]=COOPhomo[:,1]*lamb+COOPhomo[:,2]*(1-lamb)
    COOPtot[0]=np.dot(COOPhomo[:,0],SD)
    COOPtot[1]=np.dot(COOPhomo[:,1],SD)
    COOPtot[2]=np.dot(COOPhomo[:,2],SD)
    return PAYhomo,COOPhomo,COOPtot

def calcHOMO_CD(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w):
# PAYhomo: (nSx1); payoffs homogeneous populations
# COOPhome: (nSx3); 0:cooperation level, 1: cooperation level lambda=1, 2: cooperation level lambda=0
# COOPtotal: scalar; total cooperation level (taking into account SD)      
    nSTR=SD.shape[0];
    STRm=declareSTR_CD(eps);
    PAYhomo=np.zeros((nSTR)); COOPhomo=np.zeros((nSTR,3)); COOPtot=np.zeros((3))
    for i in range(0,nSTR):
        SIG=STRm[i,0:2]; ACT=STRm[i,2:6]; k=1;
        BCi,BCj=calcBC2st(SIG,ACT,SIG,ACT,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        PAYhomo[i]=np.sum(BCi*coef)
        COOPhomo[i,1]=BCi[0,1]
        COOPhomo[i,2]=BCi[1,1]
    COOPhomo[:,0]=COOPhomo[:,1]*lamb+COOPhomo[:,2]*(1-lamb)
    COOPtot[0]=np.dot(COOPhomo[:,0],SD)
    COOPtot[1]=np.dot(COOPhomo[:,1],SD)
    COOPtot[2]=np.dot(COOPhomo[:,2],SD)
    return PAYhomo,COOPhomo,COOPtot

def calcHOMO_REC(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w):
# PAYhomo: (nSx1); payoffs homogeneous populations
# COOPhome: (nSx3); 0:cooperation level, 1: cooperation level lambda=1, 2: cooperation level lambda=0
# COOPtotal: scalar; total cooperation level (taking into account SD)      
    nSTR=SD.shape[0];
    STRm=declareSTR_REC(eps);
    PAYhomo=np.zeros((nSTR)); COOPhomo=np.zeros((nSTR,3)); COOPtot=np.zeros((3))
    for i in range(0,nSTR):
        SIG=STRm[i,0:2]; ACT=STRm[i,2:6]; k=1;
        BCi,BCj=calcBC2st(SIG,ACT,SIG,ACT,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        PAYhomo[i]=np.sum(BCi*coef)
        COOPhomo[i,1]=BCi[0,1]
        COOPhomo[i,2]=BCi[1,1]
    COOPhomo[:,0]=COOPhomo[:,1]*lamb+COOPhomo[:,2]*(1-lamb)
    COOPtot[0]=np.dot(COOPhomo[:,0],SD)
    COOPtot[1]=np.dot(COOPhomo[:,1],SD)
    COOPtot[2]=np.dot(COOPhomo[:,2],SD)
    return PAYhomo,COOPhomo,COOPtot

def calcHOMO_SIG(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w):
# PAYhomo: (nSx1); payoffs homogeneous populations
# COOPhome: (nSx3); 0:cooperation level, 1: cooperation level lambda=1, 2: cooperation level lambda=0
# COOPtotal: scalar; total cooperation level (taking into account SD)      
    nSTR=SD.shape[0];
    STRm=declareSTR_SIG(eps);
    PAYhomo=np.zeros((nSTR)); COOPhomo=np.zeros((nSTR,3)); COOPtot=np.zeros((3))
    for i in range(0,nSTR):
        SIG=STRm[i,0:2]; ACT=STRm[i,2:6]; k=1;
        BCi,BCj=calcBC2st(SIG,ACT,SIG,ACT,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        PAYhomo[i]=np.sum(BCi*coef)
        COOPhomo[i,1]=BCi[0,1]
        COOPhomo[i,2]=BCi[1,1]
    COOPhomo[:,0]=COOPhomo[:,1]*lamb+COOPhomo[:,2]*(1-lamb)
    COOPtot[0]=np.dot(COOPhomo[:,0],SD)
    COOPtot[1]=np.dot(COOPhomo[:,1],SD)
    COOPtot[2]=np.dot(COOPhomo[:,2],SD)
    return PAYhomo,COOPhomo,COOPtot
    
    
def writefixMvec(fixMvec,file):
    np.save(file,fixMvec)
    return

def readfixMvec(file):
    fixMvec=np.load(file+'.npy')    
    return fixMvec
    
def doINI(N,Z,M,Q,eps,w):
    from pathlib import Path
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
    file = Path(labelfile+'.npy')
    if not file.is_file():
        print(file)
        H=calcH(N,Z)
        STRm=declareSTR(eps); # nSTR=STRm.shape[0]; #print(STm[0],STm[63])
        STATEmat,nSTATEv=declareSTATE(N); #print(STATEmat); print(nSTATEv)
        fixMvec=calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w)
        writefixMvec(fixMvec,labelfile)
        return fixMvec
    
def doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w):
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
    fixMvec=readfixMvec(labelfile)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    SD=calcSD(fixM)
    return SD

def doINI_CD(N,Z,M,Q,eps,w):
    from pathlib import Path
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_CD'
    file = Path(labelfile+'.npy')
    if not file.is_file():
        print(file)
        H=calcH(N,Z)
        STRm=declareSTR_CD(eps) 
        STATEmat,nSTATEv=declareSTATE(N)
        fixMvec=calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w)
        writefixMvec(fixMvec,labelfile)
        return fixMvec    
    
def doREST_CD(b,c,cs,lamb,beta,N,Z,M,Q,eps,w):
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_CD'
    fixMvec=readfixMvec(labelfile)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    SD=calcSD(fixM)
    return SD

def doINI_REC(N,Z,M,Q,eps,w):
    from pathlib import Path
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_REC'
    file = Path(labelfile+'.npy')
    if not file.is_file():
        print(file)
        H=calcH(N,Z)
        STRm=declareSTR_REC(eps)
        STATEmat,nSTATEv=declareSTATE(N)
        fixMvec=calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w)
        writefixMvec(fixMvec,labelfile)
        return fixMvec    
    
def doREST_REC(b,c,cs,lamb,beta,N,Z,M,Q,eps,w):
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_REC'
    fixMvec=readfixMvec(labelfile)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    SD=calcSD(fixM)
    return SD

def doINI_SIG(N,Z,M,Q,eps,w):
    from pathlib import Path
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_SIG'
    file = Path(labelfile+'.npy')
    if not file.is_file():
        print(file)
        H=calcH(N,Z)
        STRm=declareSTR_SIG(eps)
        STATEmat,nSTATEv=declareSTATE(N)
        fixMvec=calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w)
        writefixMvec(fixMvec,labelfile)
        return fixMvec    
    
def doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w):
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_SIG'
    fixMvec=readfixMvec(labelfile)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    SD=calcSD(fixM)
    return SD


def doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w):
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    STATEmat,nSTATEv=declareSTATE(N)
    PAYhomo,COOPhomo,COOPtot=calcHOMO(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w)
    return PAYhomo,COOPhomo,COOPtot

def doHOMO_CD(lamb,eps,N,M,Q,b,c,cs,SD,w):
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    STATEmat,nSTATEv=declareSTATE(N)
    PAYhomo,COOPhomo,COOPtot=calcHOMO_CD(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w)
    return PAYhomo,COOPhomo,COOPtot
def doHOMO_REC(lamb,eps,N,M,Q,b,c,cs,SD,w):
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    STATEmat,nSTATEv=declareSTATE(N)
    PAYhomo,COOPhomo,COOPtot=calcHOMO_REC(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w)
    return PAYhomo,COOPhomo,COOPtot
def doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w):
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    STATEmat,nSTATEv=declareSTATE(N)
    PAYhomo,COOPhomo,COOPtot=calcHOMO_SIG(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w)
    return PAYhomo,COOPhomo,COOPtot
   
    
def doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps):
# cs,lamb,b,M,Q    
    bigmatSD=np.zeros((nSTR,len(csV),len(lambV),len(bV),len(MV),len(QV)))
    bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    for ics in range(0,len(csV)):
        for ilamb in range(0,len(lambV)):
            for ib in range(0,len(bV)):
                for iM in range(0,len(MV)):
                    for iQ in range(0,len(QV)):
                        bigmatCOOP[ics,ilamb,ib,iM,iQ],bigmatSD[:,ics,ilamb,ib,iM,iQ]=doONEALL(beta,Z,N,c1,csV[ics],lambV[ilamb],bV[ib],MV[iM],QV[iQ],w,eps)
    return bigmatCOOP,bigmatSD

def doMATSD_SIG(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps):
# cs,lamb,b,M,Q    
    bigmatSD=np.zeros((nSTR,len(csV),len(lambV),len(bV),len(MV),len(QV)))
    bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    for ics in range(0,len(csV)):
        for ilamb in range(0,len(lambV)):
            for ib in range(0,len(bV)):
                for iM in range(0,len(MV)):
                    for iQ in range(0,len(QV)):
                        bigmatCOOP[ics,ilamb,ib,iM,iQ],bigmatSD[:,ics,ilamb,ib,iM,iQ]=doONEALL_SIG(beta,Z,N,c1,csV[ics],lambV[ilamb],bV[ib],MV[iM],QV[iQ],w,eps)
    return bigmatCOOP,bigmatSD

def doMATSD_REC(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps):
# cs,lamb,b,M,Q    
    bigmatSD=np.zeros((nSTR,len(csV),len(lambV),len(bV),len(MV),len(QV)))
    bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    for ics in range(0,len(csV)):
        for ilamb in range(0,len(lambV)):
            for ib in range(0,len(bV)):
                for iM in range(0,len(MV)):
                    for iQ in range(0,len(QV)):
                        bigmatCOOP[ics,ilamb,ib,iM,iQ],bigmatSD[:,ics,ilamb,ib,iM,iQ]=doONEALL_REC(beta,Z,N,c1,csV[ics],lambV[ilamb],bV[ib],MV[iM],QV[iQ],w,eps)
    return bigmatCOOP,bigmatSD

def doMATSD_CD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps):
# cs,lamb,b,M,Q    
    bigmatSD=np.zeros((nSTR,len(csV),len(lambV),len(bV),len(MV),len(QV)))
    bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    for ics in range(0,len(csV)):
        for ilamb in range(0,len(lambV)):
            for ib in range(0,len(bV)):
                for iM in range(0,len(MV)):
                    for iQ in range(0,len(QV)):
                        bigmatCOOP[ics,ilamb,ib,iM,iQ],bigmatSD[:,ics,ilamb,ib,iM,iQ]=doONEALL_CD(beta,Z,N,c1,csV[ics],lambV[ilamb],bV[ib],MV[iM],QV[iQ],w,eps)
    return bigmatCOOP,bigmatSD

    
def calcBIGPAY(bigmatSD,csV,lambV,MV,QV,b,c,N,eps,w):
# output: bigmatPAY has the average payoff of an specific environment (weighted over all the strategies, taking into account the SD)
# bigmatSD has to be coherent with other inputsm including their dimensions (except b)
    bigmatPAY=np.zeros((len(csV),len(lambV),len(MV),len(QV)))
    
    STATEmat,nSTATEv=declareSTATE(N)
    nSTR=bigmatSD.shape[0]
    PAYhomo=np.zeros((nSTR))
    STRm=declareSTR(eps)
    
    for ics in range(0,len(csV)):
        cs=csV[ics]
        for ilamb in range(0,len(lambV)): 
            lamb=lambV[ilamb]
            coef=np.array([[b*lamb, -c*lamb, -cs*lamb],[0*(1.-lamb), -c*(1.-lamb), -cs*(1.-lamb)]]) # assuming b=0 in nPGG
            print(['ics, ilamb: ',ics,ilamb])
            for iM in range(0,len(MV)):
                for iQ in range(0,len(QV)):
                   
                    for i in range(0,nSTR):
                        k=1
                        BCi,BCj=calcBC2st(STRm[i,0:2],STRm[i,2:6],STRm[i,0:2],STRm[i,2:6],k,N,MV[iM],QV[iQ],STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
                        PAYhomo[i]=np.sum(BCi*coef)
                    bigmatPAY[ics,ilamb,iM,iQ]=np.dot(PAYhomo,bigmatSD[:,ics,ilamb,0,iM,iQ])
                    
    return bigmatPAY   

    
def doONEALL(beta,Z,N,c1,cs1,lamb,b1,M,Q,w,eps):
    import numpy as np
    c=np.array([1., 1.])*c1
    cs=np.array([1., 1.])*cs1
    b=np.array([b1, 0.]) #*c    
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    SD=SD[:,0]
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)     
    print([cs1,lamb,b1,M,Q,COOPtot])
    return COOPtot[0], SD  # Carefull: only first component of COOPtot

def doONEALL_SIG(beta,Z,N,c1,cs1,lamb,b1,M,Q,w,eps):
    import numpy as np
    c=np.array([1., 1.])*c1
    cs=np.array([1., 1.])*cs1
    b=np.array([b1, 0.]) #*c    
    doINI_SIG(N,Z,M,Q,eps,w)
    SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    SD=SD[:,0]
    PAYhomo,COOPhomo,COOPtot=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)     
    print([cs1,lamb,b1,M,Q,COOPtot])
    return COOPtot[0], SD  # Carefull: only first component of COOPtot

def doONEALL_REC(beta,Z,N,c1,cs1,lamb,b1,M,Q,w,eps):
    import numpy as np
    c=np.array([1., 1.])*c1
    cs=np.array([1., 1.])*cs1
    b=np.array([b1, 0.]) #*c    
    doINI_REC(N,Z,M,Q,eps,w)
    SD=doREST_REC(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    SD=SD[:,0]
    PAYhomo,COOPhomo,COOPtot=doHOMO_REC(lamb,eps,N,M,Q,b,c,cs,SD,w)     
    print([cs1,lamb,b1,M,Q,COOPtot])
    return COOPtot[0], SD  # Carefull: only first component of COOPtot

def doONEALL_CD(beta,Z,N,c1,cs1,lamb,b1,M,Q,w,eps):
    import numpy as np
    c=np.array([1., 1.])*c1
    cs=np.array([1., 1.])*cs1
    b=np.array([b1, 0.]) #*c    
    doINI_CD(N,Z,M,Q,eps,w)
    SD=doREST_CD(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    SD=SD[:,0]
    PAYhomo,COOPhomo,COOPtot=doHOMO_CD(lamb,eps,N,M,Q,b,c,cs,SD,w)     
    print([cs1,lamb,b1,M,Q,COOPtot])
    return COOPtot[0], SD  # Carefull: only first component of COOPtot


def plot_COOPcslamb(bigmatCOOP,csV,lambV,bV,MV,QV):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    nr=bigmatCOOP.shape[3]; nc=bigmatCOOP.shape[4]-1 # excluding last column
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    vmin=0;vmax=1;
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):  
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatCOOP[:,:,0,iM,iQ],vmin=vmin,vmax=vmax)
            axs[iM,iQ].set_xticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"])
            axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Cooperation level')
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.savefig('cooperation.png', dpi=300)
    f.clf()
            
    return
    
    
def classST():
    ST=declareSTR(0)
    STsign=STsignonly=STmem=STmemonly=STsignmem=STs00=STs11=STs10=STs01=np.array([],int)
    for i in range(0,ST.shape[0]):
        if (ST[i,2]!=ST[i,4])or(ST[i,3]!=ST[i,5]): # sign
            STsign=np.append(STsign,i)
            if (ST[i,2]==ST[i,3])and(ST[i,4]==ST[i,5]): # no mem
                STsignonly=np.append(STsignonly,i)
            if (ST[i,2]!=ST[i,3])or(ST[i,4]!=ST[i,5]): # mem
                STsignmem=np.append(STsignmem,i)
        if (ST[i,2]!=ST[i,3])or(ST[i,4]!=ST[i,5]): # mem
            STmem=np.append(STmem,i)  
            if (ST[i,2]==ST[i,4])and(ST[i,3]==ST[i,5]): # no sign
                STmemonly=np.append(STmemonly,i) 
        if (ST[i,0]==0 and ST[i,1]==0): # 00
            STs00=np.append(STs00,i)
        if (ST[i,0]==1 and ST[i,1]==1): # 11
            STs11=np.append(STs11,i)
        if (ST[i,0]==1 and ST[i,1]==0): # 10
            STs10=np.append(STs10,i)
        if (ST[i,0]==0 and ST[i,1]==1): # 01
            STs01=np.append(STs01,i)
    return STs00,STs11,STs10,STs01,STsign, STsignonly, STmem, STmemonly, STsignmem
    


def plot_SDcslamb(label,STv,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax):
# STv: array with the strategies to agregate
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    bigmatAGR=np.sum(bigmatSD[STv,...],axis=0)
    nr=bigmatAGR.shape[3]; nc=bigmatAGR.shape[4]-1 # excluding last column
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[:,:,0,iM,iQ],vmin=vmin,vmax=vmax)
            axs[iM,iQ].set_xticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"])
            axs[iM,iQ].set_yticks([0,0.5,1.,1.5])
            axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"])
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability')
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.5, 0.95, label, ha='center',size=16)
    f.savefig(label+'.pdf', dpi=300)
    f.clf()
            
    return
    
    
def plot_SDcslambDIF(label,labup,labdown,STvpos,STvneg,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap,ext):
# STv: array with the strategies to agregate
# comap='RdBu_r' (blue to red), 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #plt.figure(1)
    bigmatAGR=np.sum(bigmatSD[STvpos,...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    nr=bigmatAGR.shape[3]; nc=bigmatAGR.shape[4] #-1 # excluding last column
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[:,:,0,iM,iQ],vmin=vmin,vmax=vmax, cmap=comap)
            axs[iM,iQ].set_xticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"])
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comap)
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)
    f.savefig(label+'.'+ext, dpi=300)
    f.clf()          
    return


def plot_SDcslambDIF_1(label,labup,labdown,STvpos,STvneg,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap):
# STv: array with the strategies to agregate
# comap='RdBu_r' (blue to red), 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    vmin=0.45
    alp=1.
    bigmatAGR=np.sum(bigmatSD[STvpos,...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR2=np.sum(bigmatSD[[52,53,54,55,60, 61, 62, 63],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR3=np.sum(bigmatSD[[48,49,50,51],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR4=np.sum(bigmatSD[[56,57,58,59],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR5=np.sum(bigmatSD[[33,35],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    nr=bigmatAGR.shape[3]; nc=bigmatAGR.shape[4] #-1 # excluding last column
    f,axs=plt.subplots(nrows=2, ncols=2, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            mins=vmin
            step=0.02
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[:,:,0,iM,iQ],np.arange(mins,1.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Reds')
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR2[:,:,0,iM,iQ],np.arange(mins,1.1,step), alpha=alp,vmin=vmin,vmax=vmax, cmap='Blues')
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR3[:,:,0,iM,iQ],np.arange(mins,1.1,step), alpha=alp,vmin=vmin,vmax=vmax, cmap='Purples')
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR4[:,:,0,iM,iQ],np.arange(mins,1.1,step), alpha=alp,vmin=vmin,vmax=vmax, cmap='Greens')
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR5[:,:,0,iM,iQ],np.arange(mins,1.1,step), alpha=alp,vmin=vmin,vmax=vmax, cmap='Greys')        
            axs[iM,iQ].set_xticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"])
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap='Greys')
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)
    f.savefig(label+'.tiff', dpi=300)
    f.savefig(label+'.pdf', dpi=300)
    f.clf()          
    return
    


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_SDcslambDIF_agre(label,groups,comapsV,nameg,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext):
# STv: array with the strategies to agregate
# comap='RdBu_r' (blue to red), 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    alp=1.
    lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
    for i in range(0,len(groups)):
        bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)
    nr=bigmatAGR.shape[4]; nc=bigmatAGR.shape[5] 
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all' )
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    comaps=comapsV
    for i in range(len(groups)):
        comaps[i]=plt.get_cmap(comapsV[i])
        comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.48,"$M=%s$" % str(MV[iM]), size=9 ,va='center')
        for iQ in range(nc-1,-1,-1):
            step=0.02
            if MV[iM]>5:   # to avoid problems with [0010**], which is two places for w=1
                rg=range(len(groups)-1,-1,-1)     
            else:
                rg=range(0,len(groups))
            for i in rg:
                h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
            axs[iM,iQ].set_xticks([0,0.5,1]); axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            axs[iM,iQ].grid(which='both', axis='both',ls='dashed')
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=9 )
                
    margbottom=0.15; margtop=0.87
    f.text(0.0, 0.5, '$c_s$', va='center', rotation='vertical',size=12)
    if nameg==0:
        margleft=0.1; margright=0.75;
        f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
        cbar_ax = f.add_axes([margright+0.1, margbottom, 1.-margleft-margright-0.12, margtop-margbottom])
        hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comaps[-1])
    else:
        margleft=0.09; margright=0.66;
        f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
        for i in range(0,len(groups)):            
            mr=0.06; hh=(margtop-margbottom)/len(groups);  hib=hh-0.11; botb=margtop-hh*(i+1)+0.109-0.027*i
            cbar_ax = f.add_axes([margright+0.11, botb, 1.-margleft-margright-0.06, hib])
            hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
            step=0.25; ti=np.arange(vmin,vmax+step,step); ti_s=["%.2f" % x for x in ti]
            hb.set_ticks(ti)
            hb.set_ticklabels(ti_s)
            cbar_ax.tick_params(labelsize=8)
            cbar_ax.set_title(nameg[i],size=8,color=mpl.cm.get_cmap(comaps[i])(1.)) 
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=12)
    f.savefig(label+'.'+ext, dpi=300)
    f.clf()       
    return




def plot_SDspace_agre(label,groups,comapsV,nameg,bigmatSDlist,yV,xV,iM,iQ,M,labup,labright,vmin,vmax,ext):
# groups: list of list with the strategies to agregate (properties: compas,nameg)
# each panel: list of list (horizontal and vertical distribution): bigmatSDlist, iQ,iM
# careful, dimensions must be coherent
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    alp=1.
    nc=len(bigmatSDlist[0]); nr=len(bigmatSDlist)

    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all' )
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    comaps=comapsV
    for i in range(len(groups)):
        comaps[i]=plt.get_cmap(comapsV[i])
        comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
    for ir in range(nr-1,-1,-1):
        axs[ir,nc-1].text(1.1,0.48,labright[ir], size=9 ,va='center',ha='left')
        for ic in range(nc-1,-1,-1):
          
            bigmatSD=bigmatSDlist[ir][ic]
            lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
            for i in range(0,len(groups)):
                bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)        
            
            step=0.02
            if M>5:   # to avoid problems with [0010**], which is two places for w=1
                rg=range(len(groups)-1,-1,-1)     
            else:
                rg=range(0,len(groups))
            for i in rg:
                h=axs[ir,ic].contourf(xV,yV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
            axs[ir,ic].set_xticks([0,0.5,1])
            axs[ir,ic].set_xticklabels(["0","0.5","1"])
            axs[ir,ic].set_yticks([0,0.5,1]);
            axs[ir,ic].set_yticklabels(["0","0.5","1"]);
            axs[ir,ic].tick_params(axis='both', which='major', labelsize=8)
            axs[ir,ic].grid(which='both', axis='both',ls='dashed')
            if ir==0:
                axs[ir,ic].set_title(labup[ic], size=9 )
                
    margbottomI=0.15; margtopI=0.87
    margbottom=0.15; margtop=0.87
    if nameg==0:
        margleft=0.1; margright=0.75;
        f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
        cbar_ax = f.add_axes([margright+0.1, margbottom, 1.-margleft-margright-0.12, margtop-margbottom])
        hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comaps[-1])
    else:
        margleft=0.15; margright=0.5;
        f.subplots_adjust(right=margright,top=margtopI,bottom=margbottomI, left=margleft)
        for i in range(0,len(groups)):          
            mr=0.06; hh=(margtop-margbottom)/len(groups);  hib=hh-0.11; botb=margtop-hh*(i+1)+0.109-0.027*i; 
            cbar_ax = f.add_axes([margright+0.15, botb, 1.-margleft-margright-0.1, hib])
            hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
            step=0.25; ti=np.arange(vmin,vmax+step,step); ti_s=["%.2f" % x for x in ti]; 
            hb.set_ticks(ti)
            hb.set_ticklabels(ti_s)
            cbar_ax.tick_params(labelsize=8)
            cbar_ax.set_title(nameg[i],size=8,color=mpl.cm.get_cmap(comaps[i])(1.))   
    f.text(margleft-0.1, (margtopI-margbottomI)/2.+margbottomI, '$c_s$', va='center', rotation='vertical',size=12)
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=12)
    f.savefig(label+'.'+ext, dpi=300)
    f.clf()           
    return



def plot_PAYcslamb(label,bigmatPAY,csV,lambV,bV,MV,QV,vmin,vmax):
# bigmatPAY[cs,lamb,M,Q] (no b)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    nr=bigmatPAY.shape[2]; nc=bigmatPAY.shape[3] # excluding last column
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            step=0.02
            PAYplt=bigmatPAY[:,:,iM,iQ]/(bV[0]*lambV) 
            h=axs[iM,iQ].contourf(lambV,csV,PAYplt,vmin=vmin,vmax=vmax, cmap='Greens')
            axs[iM,iQ].set_xticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"])
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            axs[iM,iQ].grid(which='both', axis='both',ls='dashed')
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margbottom=0.15; margtop=0.87
    f.text(0.0, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    margleft=0.1; margright=0.75;
    f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
    cbar_ax = f.add_axes([margright+0.1, margbottom, 1.-margleft-margright-0.12, margtop-margbottom])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label=r'$\overline{W}\ (\lambda rc)^{-1}$',cmap='Greens')
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.savefig(label+'.eps', dpi=300)
    f.clf()
            
    return


def plot_BAR(labup,STv,STvC,axs,beta,Z,N,M,Q,lamb,eps,w,c1,cs1,b1,vmax):

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    c=np.array([1., 1.])  *c1 #*1.* 5.    #*0.3  *0.8
    cs=np.array([1., 1.]) *cs1 #*0.1 *5.  #*0.06 *c  *0.8
    b=np.array([1., 0.]) *b1 #*20. *c   #7*c  
    
    STRmPUR=declareSTR(0)
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)

    bars=SD[STv]; n=len(STv)
    ind=np.arange(n)
    width=0.8
    
    h=axs.bar(ind,bars[:,0],width,align='center',color=STvC)
    axs.set_xlim(-1.,ind[-1]+1.)
    axs.set_ylim(0,vmax)
    axs.set_ylabel(' ')
    axs.set_yticks([0,0.2,0.4,0.6])
    axs.set_yticklabels([0,0.2,0.4,0.6],fontsize=6)
    axs.set_xticks(ind)
    axs.set_xticklabels(labup,rotation=90,fontsize=5.5,ha='center',va='top')
    [t.set_color(i) for (i,t) in zip(STvC,axs.xaxis.get_ticklabels())]
    axs.yaxis.set_ticks_position('left')
    axs.xaxis.set_ticks_position('bottom')
    title="$M=%s, Q=%s, \lambda=%s, c_s=%s$" % (str(M),str(Q),str(lamb),str(cs1/c1))
    axs.text(ind[-1]/2, vmax-0.01, title, va='top', ha='center',size=6)
                    
    return

    
def reduc(M,SD,th):
    ix=np.where(SD>=th)[0]
    iNx=np.where(SD<th)[0]
    SDN=SD[iNx][:,0]
    #print(SDN)
    sumN=np.sum(SDN)
    nSTred=len(ix)+1
    Mred=np.zeros((nSTred,nSTred))
    for i in range(0,nSTred-1):
        Mred[i,0:nSTred-1]=M[ix[i],ix]
        Mred[nSTred-1,i]=np.dot(SDN,M[iNx,ix[i]])/sumN  
        Mred[i,nSTred-1]=np.dot(SDN,M[ix[i],iNx])/sumN 
    labre=ix; labre=np.append(labre,-999999)
    np.fill_diagonal(Mred,0.)
    SDred=np.append(SD[ix],sumN)
    return Mred,SDred,labre 

def groupM(M,SD,gST):
# input: transition probability matrix (M), groups of strategies (list of lists)(gST), SD of strategies
# output: transition probability matrix and SD of groups (sorted as in gST)
    M2=np.empty([len(gST),len(M)])
    for g in range(0,len(gST)): # vertical (groups receive links)
        M2[g,:]=np.sum(M[:,gST[g]],1)  # M2 is transposed  
    M2=M2*SD[:,0]     # horizontal (groups send links)
    M2=np.transpose(M2)
    Mred=np.empty([len(gST),len(gST)])
    SDred=np.array([])
    for g in range(0,len(gST)):            
        sumN=np.sum(SD[gST[g]])
        Mred[g,:]=np.sum(M2[gST[g],:],0)/sumN
        SDred=np.append(SDred,sumN)
    return Mred,SDred 
    
    
def plotNET(b,c,cs,lamb,beta,N,Z,M,Q,eps,w,SD):
    import networkx as nx
    import matplotlib.pyplot as plt
    th=1./len(SD)
    
    STRmPUR=declareSTR(0)
    
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
    fixMvec=readfixMvec(labelfile)
    #print(fixMvec)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    
    fixMred,SDred,labre=reduc(fixM,SD,th)    
    
    labred=np.array(["%i" % x for x in labre])   
    labred[-1]="others"

    SDred[-1]=0.
    G=nx.from_numpy_matrix(fixMred,create_using=nx.DiGraph())
    
    sizenode=SDred*8000  
    
    plt.figure(1)
    plt.axis('off')
    
    pos = nx.spring_layout(G, iterations=500)
    nx.draw_networkx_nodes(G, pos , node_size=sizenode,node_color=SDred)
    
    labels={}; nodesize={}; nodecolor={};
    for inode in range(0,G.number_of_nodes()-1):
        labels[inode]=str(STRmPUR[labre[inode]].astype(int))
        nodesize[inode]=SDred[inode]*50
        if SDred[inode]>=0.1:
            nodecolor[inode]='Red'
        elif SDred[inode]>=0.05 and SDred[inode]<0.1:
            nodecolor[inode]='Blue'
        else:
            nodecolor[inode]='Green' 
    nodesize[G.number_of_nodes()-1]=1; nodecolor[G.number_of_nodes()-1]='Gray50'; labels[G.number_of_nodes()-1]='others'
    H=G;
    for u,v,d in H.edges(data=True):
        if d['weight']>0.011:
            d['c']='Gray80'
        elif d['weight']>0.009 and d['weight']<0.011:
            d['c']='RoyalBlue'
        else:
            d['c']='Gray05'
        d['weight']*=10

    nx.set_node_attributes(H, 'x_fact', nodesize); nx.set_node_attributes(H, 'y_fact', nodesize)
    nx.set_node_attributes(H, 'bc', 'Black'); nx.set_node_attributes(H, 'ic', nodecolor)
    H=nx.relabel_nodes(H,labels)
    nx.write_pajek(H, "net2.net")
    nx.draw_networkx_labels(G,pos,font_size=8,labels=labels)
    edgewidth =np.array( [ d['weight'] for (u,v,d) in G.edges(data=True)] )
    nx.draw_networkx_edges(G, pos, width=edgewidth, edge_color=edgewidth,
                           edge_vmin=0.0001,edge_vmax=0.001, arrows=True)
    
    plt.savefig('net.png', dpi=300)
    plt.clf()
      
    return fixM
    
 
def plotNETgroup(name,M,SD,labg,colg,nSTg,Z):
# create pajek file 
# nSTg: number of strategies in each group (row or column of M or SD)
    import numpy as np
    import networkx as nx

    mu=1./64 # assuming 64 strategies
    neudrift=1./Z

    G=nx.from_numpy_matrix(M,create_using=nx.DiGraph())
    H=G

    sizeg=SD*10.
    sizeg[sizeg<1]=1.

    for u,v,d in H.edges(data=True):
        if u==v:
            d['weight']=0.  
        d['weight']/=(neudrift*mu)
        if d['weight']>=(1.+0.01): 
            d['c']='Gray70'
            d['weight']=np.log10(d['weight'])
        elif d['weight']>(1.-0.01) and d['weight']<(1.+0.01): 
            d['c']='RoyalBlue'
        else:
            d['c']='Gray05'
            d['weight']=0. 

    nx.set_node_attributes(H, 'Black', 'bc'), nx.set_node_attributes(H, dict(enumerate(colg)), 'ic')
    nx.set_node_attributes(H, dict(enumerate(sizeg.astype(str))), 'x_fact'); nx.set_node_attributes(H, dict(enumerate(sizeg.astype(str))), 'y_fact')

    H=nx.relabel_nodes(H,dict(enumerate(labg)))
    nx.write_pajek(H, name+".net")

    return 




def findDRIFTgroup(fixM,Z):
    M=np.copy(fixM)
    small=0.0000001
    th=1./Z/len(fixM)
    M[M<th-small]=0.
    M[M>th+small]=0.
    groups=[]
    g=-1
    for i in range(0,len(M)):
        jg=0
        for j in range(i+1,len(M)):
            if M[i,j]!=0:
                jg+=1
                if jg==1:
                    g+=1
                    groups.append([i])
                groups[g].append(j)
    for i in range(0,len(groups)):
        if groups[i] !=-1:
            for j in range(i+1,len(groups)):
                if groups[j] !=-1:
                    if groups[i][-1] == groups[j][-1]:
                        groups[j]=-1
    groups[:] = [value for value in groups if value!=-1]
    return groups
    
    
def doONLYONE():
# run for one combination of parameters    
    beta=1.
    Z=100
    N=9
    M=5 #5
    Q=4.5 #4.5
    lamb=.5 #0.5 #0.8
    eps=0.01
    w=1.
    #H, L
    c1=0.5 #2.5
    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
    cs=np.array([1., 1.]) *0.01 *c1  #*0.06 *c  *0.8
    b=np.array([1., 0.]) *20. *c1   #7*c  
    
    STRmPUR=declareSTR(0)
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)
    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
    SSDsort=SSD[np.argsort(SSD[..., 8])] 
    for i in range(0,len(SSDsort)):
       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) 
    fixM=plotNET(b,c,cs,lamb,beta,N,Z,M,Q,eps,w,SD)
    
    return fixM, SD

