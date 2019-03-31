# -*- coding: utf-8 -*-

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
    
def declareSTATE(N):
# State of acting
# k individuls from the first strategy 
    nelem=6    
    #STATEm=np.zeros(((N-k+1)*(k+1),2))
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
            #k=0; ttt,BCjR=calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
            PAYi=H[N-1,Z-2]*np.sum(coef*BCiR)+H[N-2,Z-2]*np.sum(coef*BCi)
            PAYj=np.sum(coef*BCj)
            #print([i,j,PAYi,PAYj]) #, H[N-1,Z-2], H[N-2,Z-2]])
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
    #STATEmOPk=k-STATEm; STATEmOPNk=Nk-STATEm;
    #STATEmk=np.ones((nSTATE,2)) if k==0 else STATEm/k
    #STATEmNk=np.ones((nSTATE,2)) if k==N else STATEm/Nk
    ##SIGTTi=SIGi[0]-SIGi[1]; SIGTTj=SIGj[0]-SIGj[1]    
    ##TS=k*(ACTi[1]-ACTi[0])+Nk*(ACTj[1]-ACTj[0])
    ##TNS=k*ACTi[0]+Nk*ACTj[0]   
    #print(SIGi)
    #print(SIGj)
    ns=k*SIGi+Nk*SIGj
    consS=np.piecewise(ns,[ns<Q,ns>=Q],[0.,1.]) 
    #print([ns,Q,consS])
    nc= STATEm[:,0]+STATEm[:,1] # in the current state (that has passed)
    consA=np.piecewise(nc,[nc<M,nc>=M],[0.,1.]) 
    #print([k,consS,consA])
    Pcoopi=consA*consS*ACTi[2]+consA*(1.-consS)*ACTi[0]+(1.-consA)*consS*ACTi[3]+(1.-consA)*(1.-consS)*ACTi[1]
    Pcoopj=consA*consS*ACTj[2]+consA*(1.-consS)*ACTj[0]+(1.-consA)*consS*ACTj[3]+(1.-consA)*(1.-consS)*ACTj[1]
    for j in range(0,nSTATE):
        #print([STATEmk[j,0],STATEmOPk[j,0],consA[1]])
        #MAT[:,j]= BINOm[STATEm[j,0]]*((consA*SIGTTi+SIGi[1])**STATEm[j,0])*((1.-(consA*SIGTTi+SIGi[1]))**STATEmOPk[j,0]) \
        #        *BINOm[STATEm[j,1]]*((consA*SIGTTj+SIGj[1])**STATEm[j,1])*((1.-(consA*SIGTTj+SIGj[1]))**STATEmOPNk[j,1]) 
        ##for i in range(0,nSTATE):
        ##    MAT[i,j]=binom.pmf(STATEm[j,0],k,consA[i]*SIGTTi+SIGi[1])*binom.pmf(STATEm[j,1],Nk,consA[i]*SIGTTj+SIGj[1])
        MAT[:,j]=binom.pmf(STATEm[j,0],k,Pcoopi)*binom.pmf(STATEm[j,1],Nk,Pcoopj) 
        #print([SIGTTi,SIGi[1],SIGTTj,SIGj[1]]) 
        #print([STATEm[j,0],consA[0]*SIGTTi+SIGi[1],k,np.random.binomial(STATEm[j,0],consA[1]*SIGTTi+SIGi[1],k)])  
            #print([i, consA[i], consA[i]*SIGTTi+SIGi[1],SIGTTi,SIGi[1] ])#binom.pmf(STATEm[j,0],k,consA[i]*SIGTTi+SIGi[1]), binom.pmf(STATEm[j,1],Nk,consA[i]*SIGTTj+SIGj[1])])
    #print(MAT)
    #print(STATEm[:,0]); print(STATEm[:,1])    
    if w>=1.:
#        val,vect=lin.eigs(np.transpose(MAT),k=1,which='LR'); vect=np.real(vect/np.sum(vect))
        
        from discreteMarkovChain import markovChain
        mc=markovChain(MAT)
        mc.computePi('eigen') # We can use 'linear', 'power', 'krylov' or 'eigen'
        vect=(mc.pi).reshape(-1,1)
    else:
        vect=(1-w)*np.linalg.inv((np.identity(nSTATE)-w*MAT))[nSTATE-1,:]

    #print(nc)
    #print(consS)
    #print(consA)
    #print(vect)
    #print(nc*consA/N)
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
# modified for taking into account the perception of the benefit
    BCi1,BCj1=calcBC2st1nat(SIGi[0],ACTi,SIGj[0],ACTj,k,N,M,Q,STATEm,w)
    BCi2,BCj2=calcBC2st1nat(SIGi[1],ACTi,SIGj[1],ACTj,k,N,N+999,Q,STATEm,w) # modified: M impossible to reach
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

    #PAYi=np.zeros((N+1)); PAYj=np.zeros((N+1));
    BCki=np.zeros((N+1,2,3)); BCkj=np.zeros((N+1,2,3));
    for k in range(0,N+1):
        BCi,BCj=calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        #print(k); print(BCi); print(BCj);
        BCki[k,...]=BCi; BCkj[k,...]=BCj;
        #print([k, BCki[k,...],BCkj[k,...]])
        #PAYi[k]=lamb*(BCi[0,0]*b[0]-BCi[0,1]*c[0]-BCi[0,2]*cs[0])+(1.-lamb)*(BCi[1,0]*b[1]-BCi[1,1]*c[1]-BCi[1,2]*cs[1])
        #PAYj[k]=lamb*(BCj[0,0]*b[0]-BCj[0,1]*c[0]-BCj[0,2]*cs[0])+(1.-lamb)*(BCj[1,0]*b[1]-BCj[1,1]*c[1]-BCj[1,2]*cs[1])
        #print([k, PAYi[k], PAYj[k] ])
    #PIKi=np.zeros((Z,2,3)); PIKj=np.zeros((Z,2,3));
    #PIKiI=np.zeros((Z,2,3)); PIKjI=np.zeros((Z,2,3));
    PIKi=np.einsum('kK,kij->Kij',H[0:N,0:Z-1],BCki[1:N+1,...])
    PIKj=np.einsum('kK,kij->Kij',H[0:N,1:Z],BCkj[0:N,...])
    PIKjI=np.einsum('kK,kij->Kij',H[0:N,0:Z-1],BCkj[::-1,...][1:N+1,...])
    PIKiI=np.einsum('kK,kij->Kij',H[0:N,1:Z],BCki[::-1,...][0:N,...])
    #for K in range(1,Z-1+1):
    #    PIKi[K,:,:]=np.sum( [H[0:N,K-1]*BCki[1:N+1,:,:]],axis=0 )
    #    PIKj[K,:,:]=np.sum( [H[0:N,K]*BCkj[0:N,:,:]],axis=0 )
        #PIi[K]=np.sum( [H[0:N,K-1]*PAYi[1:N+1]] )
        #PIj[K]=np.sum( [H[0:N,K]*PAYj[0:N]] )
        #print(H[0:N,K-1]); print(np.flipud(PAYj)[1:3]);print([PAYj[1], PAYj[0]])
        #PIjI[K]=np.sum( [H[0:N,K-1]*PAYj[N-1:0+1:-1]] )
        #PIiI[K]=np.sum( [H[0:N,K]*PAYi[N:1+1:-1]] )
    #    PIKjI[K,:,:]=np.sum( [H[0:N,K-1]*np.flip(BCkj,0)[1:N+1]], axis=0 )
    #    PIKiI[K,:,:]=np.sum( [H[0:N,K]*np.flip(BCki,0)[0:N]],axis=0 )
        #print([K, PIi[K], PIj[K] ])
    #EXPK=expb**(PIi-PIj)
    #EXPKI=expb**(PIjI-PIiI)
    DIFK=PIKi-PIKj
    DIFKI=PIKjI-PIKiI
    #print(DIFK); print(DIFKI)
    
    #suma=0.
    #sumaI=0.
    #term=1.
    #termI=1.
    #sumterm=np.zeros((2,3)); sumtermI=np.zeros((2,3))
    sumterm=np.cumsum(DIFK,axis=0); sumtermI=np.cumsum(DIFKI,axis=0)
    #print('-----')
    #print(sumterm)
    #print('-----')
    #for m in range(0,Z-2+1): #range(1,Z-1+1):
    #    sumterm[m,:,:]=DIFK[m,:,:]
    #    sumtermI[m,:,:]=DIFKI[m,:,:]
        #term*=EXPK[m]
        #termI*=EXPKI[m]
        #suma+=term
        #sumaI+=termI
    return sumterm, sumtermI    

def calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H,w):
#  from i to j (i->j)
    nSTR=STRm.shape[0]; 
    fixMvec=np.zeros((nSTR,nSTR,Z-1,2,3))
    for i in range(0,nSTR):
        for j in range(i+1,nSTR):
          #if i==0 and j==63:
            print([i, j])
            pfixvec,pfixIvec=calcFIX1vec(i,j,STRm,N,Z,M,Q,STATEmat,nSTATEv,H,w)
            fixMvec[j,i,...]=pfixvec; fixMvec[i,j,...]=pfixIvec
            #print(fixMvec[i,j,...])
    return fixMvec

    
    
def calcFIXM(coef,expb,Z,fixMvec):
#  from i to j (i->j)
# fixMvec(2,3)   
    #expb=np.array([[np.exp(-beta*r[0]*c[0]*(1.-lamb)), np.exp(-beta*c[0]*(1.-lamb)), np.exp(-beta*cs[0])*(1.-lamb)], \
    #                [np.exp(-beta*r[1]*c[1]*lamb), np.exp(-beta*c[1]*lamb), np.exp(-beta*cs[1])*lamb]])   
    #print(shape)
    fixM=1./(1.+np.sum(expb**np.einsum('ijmab,ab->ijm',fixMvec,coef),axis=2))
    #print(fixMvec[0,63,...])
    #print(coef)
    #print([fixM[0,63], fixM[63,0]])
    np.fill_diagonal(fixM, 0.) 
    nSTR=len(fixM)
    

    fixM=fixM/nSTR
    suma=np.sum(fixM,axis=1)
    fixM[range(nSTR),range(nSTR)]=1.-suma
    
#    suma=np.sum(fixM,axis=1)   
#    fixM[range(nSTR),range(nSTR)]=1.-suma/nSTR
#    
#    print(fixM)
#    np.savetxt('ttt.dat',fixM,delimiter=', ',newline='],\n')
#    

    #print(np.sum(fixM,axis=1))
    
    #nSTR=fixMvec.shape[0];
    #fixM=np.zeros((nSTR,nSTR))
    #for i in range(0,nSTR):
    #    for j in range(0,nSTR):
    #        if(i!=j):
    #            
    #            
    #            np.einsum('ab,kij->Kij',H[0:N,0:Z-1],np.flip(BCkj,0)[1:N+1,:,:])
    #            
    #            suma=0
    #            for m in range(0,Z-2+1): #range(1,Z-1+1):  #bad. it should be 0->Z-2 (+1)
    #                suma+=np.prod(coef**fixMvec[i,j,m,:,:]) 
    #            fixM[i,j]=1./(1.+suma)
    #    fixM[i,i]=1.-np.sum(fixM[i,:])/nSTR
#    print(np.sum(fixM,axis=1))
    return fixM


def calcSD(fixM):
#    import scipy.sparse.linalg as lin
#    vals,vecs=lin.eigs(np.transpose(fixM),k=1,which='LR',tol=1e-12)
#    vecsabs=np.real(np.absolute(vecs))
#    SD=vecsabs/np.sum(vecsabs)
    
    from discreteMarkovChain import markovChain
    mc=markovChain(fixM)
    mc.computePi('eigen') # We can use 'linear', 'power', 'krylov' or 'eigen'
    SD=(mc.pi).reshape(-1,1)
   
    return SD
    
def calcHOMO(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w):
# PAYhomo: (nSx1); payoffs homogeneous populations
# COOPhome: (nSx3); 0:cooperation level, 1: cooperation level lambda=1, 2: cooperation level lambda=0
# COOPtotal: scalar; total cooperation level (taking into account SD)      
    nSTR=SD.shape[0];
    STRm=declareSTR(eps);
    PAYhomo=np.zeros((nSTR)); COOPhomo=np.zeros((nSTR,3))
    for i in range(0,nSTR):
        SIG=STRm[i,0:2]; ACT=STRm[i,2:6]; k=1;
        BCi,BCj=calcBC2st(SIG,ACT,SIG,ACT,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]],w)
        PAYhomo[i]=np.sum(BCi*coef)
        COOPhomo[i,1]=BCi[0,1]
        COOPhomo[i,2]=BCi[1,1]
    COOPhomo[:,0]=COOPhomo[:,1]*lamb+COOPhomo[:,2]*(1-lamb)
    COOPtot=np.dot(COOPhomo[:,0],SD)
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
    #print(fixMvec)
    fixM=calcFIXM(coef,expb,Z,fixMvec)
    #print(fixM)
    SD=calcSD(fixM)
    return SD


def doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w):
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    STATEmat,nSTATEv=declareSTATE(N)
    PAYhomo,COOPhomo,COOPtot=calcHOMO(coef,lamb,eps,N,M,Q,STATEmat,nSTATEv,SD,w)
    return PAYhomo,COOPhomo,COOPtot
    

def doMATCOOP(csV,lambV,bV,MV,QV):
# cs,lamb,b,M,Q    
    bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    for ics in range(0,len(csV)):
        for ilamb in range(0,len(lambV)):
            for ib in range(0,len(bV)):
                for iM in range(0,len(MV)):
                    for iQ in range(0,len(QV)):
                        bigmatCOOP[ics,ilamb,ib,iM,iQ],SD=doONEALL(csV[ics],lambV[ilamb],bV[ib],MV[iM],QV[iQ],w)
    return bigmatCOOP
    
    
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

        #PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD)
    
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

    #H, L
    c=np.array([1., 1.])*c1
    cs=np.array([1., 1.])*cs1
    b=np.array([b1, 0.]) #*c    
    #eps=0.01
    #STRmPUR=declareSTR(0.)
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    SD=SD[:,0]
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)     
    print([cs1,lamb,b1,M,Q,COOPtot])
    return COOPtot, SD

def plot_COOPcs(bigmatCOOP,csV):
    import matplotlib.pyplot as plt
    plt.figure(1)
    for ilamb in range(0,len(lambV)):
        plt.plot(csV,bigmatCOOP[:,ilamb,0,0,0])
    plt.ylabel('Cooperation level')
    plt.xlabel('c_s')
    return

def plot_COOPcslamb(bigmatCOOP,csV,lambV,bV,MV,QV):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #plt.figure(1)
    nr=bigmatCOOP.shape[3]; nc=bigmatCOOP.shape[4]-1 # excluding last column
#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    vmin=0;vmax=1;
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):  
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatCOOP[:,:,0,iM,iQ],vmin=vmin,vmax=vmax)
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Cooperation level')
    #hb.set_ticks(np.linspace(0,1,11))
#    plt.show()
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
    #plt.figure(1)
    bigmatAGR=np.sum(bigmatSD[STv,...],axis=0)
    nr=bigmatAGR.shape[3]; nc=bigmatAGR.shape[4]-1 # excluding last column
#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[:,:,0,iM,iQ],vmin=vmin,vmax=vmax)
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability')
    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.5, 0.95, label, ha='center',size=16)
    f.savefig(label+'.pdf', dpi=300)
    #f.savefig(label+'.png', dpi=300)
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
#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[:,:,0,iM,iQ],vmin=vmin,vmax=vmax, cmap=comap)
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comap)
    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)
    #f.text(0.5, 0.95, label, ha='center',size=16)
    #f.savefig(label+'.png', dpi=300)
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
    #plt.figure(1)
    bigmatAGR=np.sum(bigmatSD[STvpos,...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR2=np.sum(bigmatSD[[52,53,54,55,60, 61, 62, 63],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR3=np.sum(bigmatSD[[48,49,50,51],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR4=np.sum(bigmatSD[[56,57,58,59],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    bigmatAGR5=np.sum(bigmatSD[[33,35],...],axis=0)-np.sum(bigmatSD[STvneg,...],axis=0)
    nr=bigmatAGR.shape[3]; nc=bigmatAGR.shape[4] #-1 # excluding last column
#    f=plt.figure(1,figsize=(20,20))
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
          
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
    margleft=0.13; margright=0.75
    f.subplots_adjust(right=margright,top=0.87,bottom=0.15, left=margleft)
    cbar_ax = f.add_axes([0.85, 0.13, 0.05, 0.77])
    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap='Greys')
    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)
    f.text(0.04, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)
    #f.text(0.5, 0.95, label, ha='center',size=16)
    #f.savefig(label+'.png', dpi=300)
    #f.savefig(label+'.svg', dpi=300)
    f.savefig(label+'.tiff', dpi=300)
    f.savefig(label+'.pdf', dpi=300)
    f.clf()
            
    return
    


def plot_SDcslambDIF_agre(label,groups,comaps,nameg,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext,Md):
# STv: array with the strategies to agregate
# comap='RdBu_r' (blue to red), 
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    alp=1.
    lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
    for i in range(0,len(groups)):
        bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)
    nr=bigmatAGR.shape[4]; nc=bigmatAGR.shape[5] 
#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all' )
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,1,"$M=%s$" % str(MV[iM]), size=10 ,va='center')
        for iQ in range(nc-1,-1,-1):
            step=0.02
            if MV[iM]>Md:   # to avoid problems with [0010**]
                rg=range(len(groups)-1,-1,-1)     
            else:
                rg=range(0,len(groups))
            for i in rg:
                h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
            axs[iM,iQ].grid(which='both', axis='both',ls='dashed')
            if iM==0:
                axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=10 )
                
    margbottom=0.15; margtop=0.87
    f.text(0.0, 0.5, '$c_s$', va='center', rotation='vertical',size=16)
    if nameg==0:
        margleft=0.1; margright=0.75;
        f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
        cbar_ax = f.add_axes([margright+0.1, margbottom, 1.-margleft-margright-0.12, margtop-margbottom])
        hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comaps[-1])
    else:
        margleft=0.08; margright=0.66;
        f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
        for i in range(0,len(groups)):
            
            mr=0.06; hh=(margtop-margbottom)/len(groups);  hib=hh-0.11; botb=margtop-hh*(i+1)+0.11-0.015*i; 
            
            #botb=(margtop-margbottom)/2.+(i-np.floor(len(groups)/2.))*0.2 ; hib=0.03
            cbar_ax = f.add_axes([margright+0.11, botb, 1.-margleft-margright-0.06, hib])
            hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
            step=0.2; ti=np.arange(vmin,vmax+step,step); ti_s=["%.1f" % x for x in ti];  # ti_s[0]='<'+ti_s[0]
            hb.set_ticks(ti)
            hb.set_ticklabels(ti_s)
            cbar_ax.tick_params(labelsize=8)
            cbar_ax.set_title(nameg[i],size=8,color=mpl.cm.get_cmap(comaps[i])(1.))
    
    f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=16)


    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    #f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    #f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)

    #for i in range(0,len(ext)):
    f.savefig(label+'.'+ext, dpi=300)
    f.clf()
            
    return



def plot_PAYcslamb(label,bigmatPAY,csV,lambV,bV,MV,QV,vmin,vmax):
# bigmatPAY[cs,lamb,M,Q] (no b)
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #plt.figure(1)
    nr=bigmatPAY.shape[2]; nc=bigmatPAY.shape[3] # excluding last column
#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    #vmax=10
    #print(csV[15],lambV[15],bV[0]*lambV[15],bigmatPAY[15,15,3,0])

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    for iM in range(nr-1,-1,-1):
        axs[iM,nc-1].text(1.1,0.4,"$M=%s$" % str(MV[iM]), size=10 )
        for iQ in range(nc-1,-1,-1):
            step=0.02
            PAYplt=bigmatPAY[:,:,iM,iQ]/(bV[0]*lambV) #np.transpose(np.array(lambV)[np.newaxis])
            h=axs[iM,iQ].contourf(lambV,csV,PAYplt,vmin=vmin,vmax=vmax, cmap='Greens')
            axs[iM,iQ].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[iM,iQ].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
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
    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    #f.text(0.5, 0.95, label, ha='center',size=16)
    #f.savefig(label+'.png', dpi=300)
    f.savefig(label+'.eps', dpi=300)
    f.clf()
            
    return


def plot_BAR(labup,STv,STvC,axs,beta,Z,N,M,Q,lamb,eps,w,c1,cs1,b1,vmax):

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    c=np.array([1., 1.])  *c1 #*1.* 5.    #*0.3  *0.8
    cs=np.array([1., 1.]) *cs1 #*0.1 *5.  #*0.06 *c  *0.8
    b=np.array([1., 0.]) *b1 #*20. *c   #7*c  
    
    STRmPUR=declareSTR(0); # nSTR=STRmPUR.shape[0];
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)

    #nr=1; nc=1
    #f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all')
    #f.subplots_adjust(hspace=0.2, wspace=0.2)
    ###f = plt.figure()
    ###axs = f.add_subplot(ncol,nrow,npl)

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
    #title="$M=%s, Q=%s, \lambda=%s, c_s=%s, b=%s$" % (str(M),str(Q),str(lamb),str(cs1/c1),str(b1/c1))
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
        Mred[nSTred-1,i]=np.dot(SDN,M[iNx,ix[i]])/sumN  # it may be wrong the /sumN
        Mred[i,nSTred-1]=np.dot(SDN,M[ix[i],iNx])/sumN  # it may be wrong the /sumN
    #labred=np.array(["%i" % x for x in ix])   
    #labred=np.append(labred,'others')
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
    #print(fixMred)
    #print(labre)
    #print(SDred)
    
    labred=np.array(["%i" % x for x in labre])   
    labred[-1]="others"

    SDred[-1]=0.
    G=nx.from_numpy_matrix(fixMred,create_using=nx.DiGraph())
    #print([fixM[28,0],fixM[0,28]])
    #print([G[28][0],G[0][28]])
    
    sizenode=SDred*8000  #np.log(SDred*len(SD)*0.8)*500
    #nx.draw_networkx(G,pos=nx.spring_layout(G,scale=2),
    #                 node_size=sizenode,node_color=SD,
    #                 width=0.05,linewidth=100)
    
    plt.figure(1)
    #plt.subplot(211);
    plt.axis('off')
    
    #pos = nx.circular_layout(G)
    pos = nx.spring_layout(G, iterations=500)
    nx.draw_networkx_nodes(G, pos , node_size=sizenode,node_color=SDred)
    
    labels={}; nodesize={}; nodecolor={};
    for inode in range(0,G.number_of_nodes()-1):
        #labels[inode]=labred[inode]
        labels[inode]=str(STRmPUR[labre[inode]].astype(int))
        nodesize[inode]=SDred[inode]*50
        if SDred[inode]>=0.1:
            nodecolor[inode]='Red'
        elif SDred[inode]>=0.05 and SDred[inode]<0.1:
            nodecolor[inode]='Blue'
        else:
            nodecolor[inode]='Green' 
    nodesize[G.number_of_nodes()-1]=1; nodecolor[G.number_of_nodes()-1]='Gray50'; labels[G.number_of_nodes()-1]='others'
    #edgecolor={};
    #print(fixMvec)
    H=G;
    for u,v,d in H.edges(data=True):
        #edgecolor[iedge]='Red' if d[iedge]>=0.01 else 'Gray';
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
     #print(fixMred)
    #print(labre)
    #print(SDred)
    edgewidth =np.array( [ d['weight'] for (u,v,d) in G.edges(data=True)] )
    nx.draw_networkx_edges(G, pos, width=edgewidth, edge_color=edgewidth,
                           edge_vmin=0.0001,edge_vmax=0.001, arrows=True)
    
    plt.savefig('net.png', dpi=300)
    plt.clf()
    
    
    return fixM
    
 
def plotNETgroup(name,M,SD,labg,colg,nSTg,Z):
# create pajek file 
# nSTg: number of strategies in each group (row or column of M or SD)
    import networkx as nx
    import matplotlib.pyplot as plt

    #print(nSTg)
    ##drift=(1./Z/64)*np.array(nSTg)  # drift changes depending on the number of strategies in each receptor group; assumed 64 strategies in total
    ##M=Mi/drift

    #print(M[4,0])

    mu=1./64 # assuming 64 strategies
    neudrift=1./Z

    G=nx.from_numpy_matrix(M,create_using=nx.DiGraph())
    H=G
    
     #np.log(SDred*len(SD)*0.8)*500
    #nx.draw_networkx(G,pos=nx.spring_layout(G,scale=2),
    #                 node_size=sizenode,node_color=SD,
    #                 width=0.05,linewidth=100)
    

    
    #labels={}; nodesize={}; nodecolor={};
    #for inode in range(0,G.number_of_nodes()-1):
        #labels[inode]=labred[inode]
        #labels[inode]=str(STRmPUR[labre[inode]].astype(int))
        #nodesize[inode]=SDred[inode]*50
        #if SDred[inode]>=0.1:
        #    nodecolor[inode]='Red'
        #elif SDred[inode]>=0.05 and SDred[inode]<0.1:
        #    nodecolor[inode]='Blue'
        #else:
        #    nodecolor[inode]='Green' 
    
    ##nodesize[G.number_of_nodes()-1]=1;# nodecolor[G.number_of_nodes()-1]='Gray50'; labels[G.number_of_nodes()-1]='others'
    #edgecolor={};
    #print(fixMvec)
    sizeg=SD*10.
    sizeg[sizeg<1]=1.


    
    H=G;
    ##drift=1./100/64  # assuming Z=100, and transition probabilities were divided by N_str=64
    #small=0.001 #drift*0.01
    for u,v,d in H.edges(data=True):
        #edgecolor[iedge]='Red' if d[iedge]>=0.01 else 'Gray';
        if u==v:
            d['weight']=0.  
        d['weight']/=(neudrift*mu)
        if d['weight']>=(1.+0.01): #1.+small: #.drift+small:
            d['c']='Gray70'
            #d['weight']=.3
            #d['weight']*=50.
        elif d['weight']>(1.-0.01) and d['weight']<(1.+0.01): #d['weight']>1.-small and d['weight']<1.+small: 
            d['c']='RoyalBlue'
        else:
            d['c']='Gray05'
            d['weight']=0. 
        




    nx.set_node_attributes(H, 'Black', 'bc'), nx.set_node_attributes(H, dict(enumerate(colg)), 'ic')
    nx.set_node_attributes(H, dict(enumerate(sizeg.astype(str))), 'x_fact'); nx.set_node_attributes(H, dict(enumerate(sizeg.astype(str))), 'y_fact')

#    nx.set_node_attributes(H, {1:'ee'}, 'ic')
#    nx.set_node_attributes(H, dict(enumerate(sizeg)), 'x_fact'); nx.set_node_attributes(H, dict(enumerate(sizeg)), 'y_fact')


#    nx.set_node_attributes(H, sizeg, 'x_fact'); nx.set_node_attributes(H, sizeg, 'y_fact')
#    nx.set_node_attributes(H, 'Black', 'bc'); nx.set_node_attributes(H, colg, 'ic')
    #print(labg)
    H=nx.relabel_nodes(H,dict(enumerate(labg)))
    nx.write_pajek(H, name+".net")


#    plt.figure(1)
#    #plt.subplot(211);
#    plt.axis('off')
#    
#    #pos = nx.circular_layout(G)
#    pos = nx.spring_layout(G, iterations=500)
#    nx.draw_networkx_nodes(G, pos , node_size=sizenode,node_color=SDred)    
#    nx.draw_networkx_labels(G,pos,font_size=8,labels=labg)
#     #print(fixMred)
#    #print(labre)
#    #print(SDred)
#    edgewidth =np.array( [ d['weight'] for (u,v,d) in G.edges(data=True)] )
#    nx.draw_networkx_edges(G, pos, width=edgewidth, edge_color=edgewidth,
#                           edge_vmin=0.0001,edge_vmax=0.001, arrows=True)
#    
#    plt.savefig(name+'.png', dpi=300)
#    plt.clf()
    
    
    return fixM




def findDRIFTgroup(fixM,Z):
    M=np.copy(fixM)
    small=0.0000001
    th=1./Z/len(fixM)
    M[M<th-small]=0.
    M[M>th+small]=0.
      
    #print(M)
    
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
    
    beta=1.
    Z=100
    N=9
    M=5 #5
    Q=4.5 #4.5
    lamb=.5 #0.5 #0.8
    eps=0.01
    w=0.9
    #H, L
    c1=1. #2.5
    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
    cs=np.array([1., 1.]) *0.01 *c1  #*0.06 *c  *0.8
    b=np.array([1., 0.]) *20. *c1   #7*c  
    
    STRmPUR=declareSTR(0); # nSTR=STRmPUR.shape[0];
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)
    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
    SSDsort=SSD[np.argsort(SSD[..., 8])] 
    for i in range(0,len(SSDsort)):
       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
    #print(COOPtot)
    fixM=plotNET(b,c,cs,lamb,beta,N,Z,M,Q,eps,w,SD)
    
    return fixM, SD


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plot_SDspace_agre(label,groupsV,comapsV,nameg,bigmatSDlist,yV,xV,iM,iQ,M,labup1,labup2,labright,vmin,vmax,ext):
# groups: list of list with the strategies to agregate (properties: compas,nameg)
# each panel: list of list (horizontal and vertical distribution): bigmatSDlist
# different rows: different iM
# careful, dimensions must be coherent

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    alp=1.
    nc=len(bigmatSDlist[0]); nr=len(bigmatSDlist)

#    f=plt.figure(1,figsize=(20,20))
    f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all' )
    f.subplots_adjust(hspace=0.2, wspace=0.2)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    comaps=comapsV
    for i in range(len(groupsV[0][0])):
        comaps[i]=plt.get_cmap(comapsV[i])
        comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
    for ir in range(nr-1,-1,-1):
        axs[ir,nc-1].text(1.15,0.48,labright[ir], size=12 ,va='center',ha='left',rotation=0)
        for ic in range(nc-1,-1,-1):
          
            groups=groupsV[ir][ic]
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
                h=axs[ir,ic].contourf(xV,yV,bigmatAGR[i,:,:,0,iM[ir],iQ[ir][ic]],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
            axs[ir,ic].set_xticks([0,0.5,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
            axs[ir,ic].set_xticklabels(["0","0.5","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
            axs[ir,ic].set_yticks([0,0.5,1]);
            axs[ir,ic].set_yticklabels(["0","0.5","1"]);
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
            #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
            axs[ir,ic].tick_params(axis='both', which='major', labelsize=10)
            axs[ir,ic].grid(which='both', axis='both',ls='dashed')
            if ir==0:
                axs[ir,ic].set_title(labup2[ic], pad=15,size=12)
                
    margbottomI=0.27; margtopI=0.82;
    margleft=0.1; margright=0.9;

    f.subplots_adjust(right=margright,top=margtopI,bottom=margbottomI, left=margleft)
    
    f.text(0, (margtopI-margbottomI)/2.+margbottomI, '$c_s$', va='center', rotation='vertical',size=16)
    f.text((margright-margleft)/2+margleft, margbottomI-0.12, '$\lambda$', ha='center',size=16)
    f.text((margright-margleft)/4+margleft, 0.95, labup1[0], ha='center',size=12)
    f.text((margright-margleft)*3/4+margleft, 0.95, labup1[1], ha='center',size=12)

    #hb.set_ticks(np.linspace(vmin,vmax,11))
#    plt.show()
    #f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
    #f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)

    #for i in range(0,len(ext)):
    f.savefig(label+'.'+'svg', dpi=300)
    f.clf()
    
    
  
    
    margbottomI=0.25; margtopI=0.75
    margbottom=0.1; margtop=0.95
    if nameg==0:
        margleft=0.1; margright=0.75;

    else:
        #margleft=0.08; margright=0.66;
        margleft=0.1; margright=0.6;
        for i in range(0,len(groups)):
            
            mr=0.06; hh=(margtop-margbottom)/len(groups);  hib=hh-0.09; botb=margtop-hh*(i+1)+0.09-0.02*i; 
            
            #botb=(margtop-margbottom)/2.+(i-np.floor(len(groups)/2.))*0.2 ; hib=0.03
            cbar_ax =f.add_axes([margright+0.11, botb, 1.-margleft-margright-0.047, hib])
            hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
            step=0.25; ti=np.arange(vmin,vmax+step,step); ti_s=["%.2f" % x for x in ti];  # ti_s[0]='<'+ti_s[0]
            hb.set_ticks(ti)
            hb.set_ticklabels(ti_s)
            cbar_ax.tick_params(labelsize=8)
            cbar_ax.set_title(nameg[i],size=9,color=mpl.cm.get_cmap(comaps[i])(1.))
            
    f.savefig(label+'_legend.'+'svg', dpi=300)
    f.clf()            
            
            
    return


if __name__ == "__main__":
    import numpy as np
    import time; import timeit

    #doONLYONE()

    gSC=[20,28] # SC
    gSCm=[22,30] # SC mut
    gSD=[33,35] # SD
    gSDm=[41,43] # SD mut
    gSF=[48,49,50,51] # SF
    gSFm=[56,57,58,59] # SF mut
    gMBc=[52, 53, 54, 55] # MB C
    gMBd=[60, 61, 62, 63] # MB D
    gMB=gMBc+gMBd
    gALL=list(range(0,64))
    gG=gSC+gSCm+gSD+gSDm+gSF+gSFm+gMBc+gMBd
    gNO = [x for x in gALL if x not in gG]
    gST=[gSC,gSCm,gSD,gSDm,gSF,gSFm,gMBc,gMBd,gNO]
    gS11=list(range(0,16))
    gS10=list(range(16,32))
    gS01=list(range(32,48))
    gS00=list(range(48,64))
    #print(gST)
    
    gSC1=[29,28]
    gSD1=[39,35]
    
    gSCt=[20,29,28]
    gSDt=[33,39,35]
    gGt=gSCt+gSDt+gSF+gSFm+gMBc+gMBd
    gNOt=[x for x in gALL if x not in gGt]

    gNO2 = [x for x in gALL if x not in gSC+gSD+gSF+gSFm+gMBc+gMBd ]    
    
    
    
###################  PLOTS summary  ############################################################################    

 
    lambV=np.linspace(0,1,31)   
    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
    bVo=np.array([10.])
    
    beta=1.
    Z=100
    N=9
    eps=0.01   

    c1=1.
    csV= csVo
    bV=c1* bVo

    
    bigmatSD_abun_1=np.load('file_SD_N9_beta1_b10_w1_bper_abun_1.npy') 
    bigmatSD_none_1=np.load('file_SD_N9_beta1_b10_w1_bper_none_1.npy') 
    bigmatSD_self_1=np.load('file_SD_N9_beta1_b10_w1_SELF_NEW_1_b.npy')
    g_abun_1=[gMB,gMBc,gSF,gSFm,[28,24]+[35,34],[22,30, 43,41]]
    g_none_1=[gMB,gMBc,gSF,gSFm,gSC+gSD,[22,30, 43,41]]
    g_self_1=[gMB,[],gSF+gMBc,gSFm,gSCt+gSDt,[22,30, 43,41]]
    g_self_1b=[gMB,[],gSF,gSFm,gSCt+gSDt,[22,30, 43,41]]
    
    labup1=['Benefit-perceived', 'Self-aware']
    labup2=['No resource','Abundance','R=M-1','R=M']
    labright=['M=5','M=7']
    bigmatSDlist=[[bigmatSD_none_1,bigmatSD_none_1],[bigmatSD_abun_1,bigmatSD_abun_1],[bigmatSD_self_1,bigmatSD_self_1],[bigmatSD_self_1,bigmatSD_self_1]]
    bigmatSDlist=list(map(list, zip(*bigmatSDlist))) # transposing list
    groupslist=[[g_none_1,g_none_1],[g_abun_1,g_abun_1],[g_self_1,g_self_1b],[g_self_1b,g_self_1b]]
    groupslist=list(map(list, zip(*groupslist))) # transposing list
    #iQ=np.zeros((3,2),int); iM=np.zeros((3,2),int)
    iQ=[[2,2,1,2],[2,2,3,4]]; iM=[2,3]
    
    
    ext='eps'
    vmin=0.5; vmax=1.
    comaps=['Blues','Oranges','Purples','Greens','Reds','Greys']
    ngroups=[r'FR-D',r'FR-F',r'FR-C',r'FR-O', r'SC + SD',r'SC-O + SD-O']

    plot_SDspace_agre('SD_agre_bper-self_2',groupslist,comaps,ngroups,bigmatSDlist,csV,lambV,iM,iQ,5,labup1,labup2,labright,vmin,vmax,ext)
    
  
#######    
    
    
    
    
    
######  Extracting grpah of invasions ##########################
#
####----- separating strategies ------
#    gST=[[28]]+[[20]]+[[35]]+[[33]]+[[i] for i in gSF]+[[i] for i in gSFm]+[[i] for i in gMBc] +[[i] for i in gMBd] +[gNO2]    
#    #gST=[[28]]+[[20]]+[[29]]+[[35]]+[[33]]+[[39]]+[[i] for i in gSF]+[[i] for i in gSFm]+[[i] for i in gMBc] +[[i] for i in gMBd] #+[gNOt]
#    colg= ['Red']*2+['Grey30']*2+['Mulberry']*len(gSF)+['Green']*len(gSFm)+['Cyan']*len(gMBc)+['Blue']*len(gMBd) +['Gray05']
#
#    labg=['']*len(gST)
#    STRmPUR=declareSTR(0)
#    for i in range(0,len(gST)):
#        labg[i]=str(STRmPUR[gST[i][0]].astype(int))
#    labg[-1]="others"
#####-------------------------------
#
####----- groups of strategies ------
#    gST=[gSC,gSD,gSF,gSFm,gMBc+gMBd,gNO2]
#    colg= ['Red']+['Gray30']+['Purple']+['Green']+['NavyBlue'] +['Gray05']
#    labg= ['SC']+['SD']+['FR-C']+['FR-01']+['FR-D + FR-10'] +['others']
#####-------------------------------    
#    
#####----- test ------
##    gST=[[28],[35]]
##    colg= ['Red']+['Gray30'] #+['Purple']+['Green']+['NavyBlue'] +['Gray05']
##    labg= ['SC']+['SD'] #+['FR-C']+['FR-01']+['FR-D + FR-10'] +['others']
######-------------------------------   
#
#    beta=1.
#    Z=100
#    N=9
#    M=5
#    Q=4.5
#    lamb=0.6 #0.5 #0.8
#    cs0=1.
#    eps=0.01
#    w=0.9
#    #H, L
#    c1=0.5 #2.5
#    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#    cs=np.array([1., 1.]) *cs0 *c1  #*0.06 *c  *0.8
#    b=np.array([1., 0.]) *20. *c1   #7*c  
#    
#    expb=np.exp(-beta)
#    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
#    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
#    fixMvec=readfixMvec(labelfile)
#    fixM=calcFIXM(coef,expb,Z,fixMvec)  
#    doINI(N,Z,M,Q,eps,w)
#    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#    
#    Mred,SDred=groupM(fixM,SD,gST)
#    SDred[SDred<0]=0.; Mred[Mred<0]=0. # correct very small negative values   
#
#    STRmPUR=declareSTR(0)
#    #print(STRmPUR[28],STRmPUR[35])
#    #print(fixM[28,35]*64*100,fixM[35,28]*64*100)
#    print(fixM[:,63]*64*100)
#    print(fixM[:,28]*64*100)
#    print((fixM[:,28]-fixM[:,63])*64*100)
#    print(SD[28],SD[63])
#    print(Mred[4,0]*64*100)
#
#
#    nSTg=[len(g) for g in gST]
#    name='NET_M_'+str(M)+'_Q_'+str(Q)+'_lamb_'+str(lamb)+'_cs_'+str(cs0)
#    plotNETgroup(name,Mred,SDred,labg,colg,nSTg,Z)
#    
#    
#    
################################################################   
    
    
    


########### Limits drifting groups #########
#    from scipy.stats import binom
#    import matplotlib.pyplot as plt
#    from decimal import Decimal
#    
#    Z=100
#    N=9
#    eps=0.01
#    M=np.array(range(1,N))
#    
#    Pr_lessM=binom.cdf(M-1,N,eps)*(1.-eps**M) 
#    Rlim=np.log(1.-1./Z)/np.log(Pr_lessM) +1 # R < R_lim ---> xx10xx equivalent to xx00xx (assuming a tolerance of 1/Z). All start D, and they never enter in Nc>=M by mistake
#    wlim=1.-1./Rlim
#
#    print(Pr_lessM)
#    print(Rlim)
#    print(wlim)
#    
##    f = plt.figure()
##    for i in range(0,len(cs1V)):
##        nrow=len(cs1V)
##        ncol=1
##        npl=i+1     
##        axs = f.add_subplot(nrow,ncol,npl)
##        if npl!=nrow:
##            labx=[]
##        else:
##            labx=labup
##        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i],QV[i],lambV[i],epsV[i],w,c1,cs1V[i],bV[i],vmax)
##    f.text(0.05, 0.5, 'Stationary Distribution', rotation=90, va='center', ha='center',size=8)
##    f.savefig(label+'.eps', dpi=300)
##    f.clf()
##
#
#
#    f = plt.figure()
#    for tol in [1e-2, 1e-3, 1e-6, 1e-9]:
#        Rlim=np.log(1.-tol)/np.log(Pr_lessM) +1 # R < R_lim ---> xx10xx equivalent to xx00xx (assuming a tolerance of 1/Z). All start D, and they never enter in Nc>=M by mistake
#        wlim=1.-1./Rlim    
#        #axs=plt.subplot(221); plt.semilogy(M,wlim); axs.set_xlim(0, N+1); axs.set_ylim(0.9, 1.01); axs.set_yticks([0.9,0.99,1]);  axs.set_xticks(range(0,N+1))
#        #axs=plt.subplot(222); plt.plot(M,wlim); axs.set_xlim(0, N+1); axs.set_ylim(0.9, 1.01); axs.set_yticks([0.9,0.99,1]);  axs.set_xticks(range(0,N+1))
#        #axs=plt.subplot(223); plt.semilogy(M,Rlim); axs.set_xlim(0, N+1); axs.set_ylim(1, 100000); axs.set_yticks([10,100]);  axs.set_xticks(range(0,N+1))
#        #axs=plt.subplot(224); plt.plot(M,Rlim); axs.set_xlim(0, N+1); axs.set_ylim(1, 100000); axs.set_yticks([10,100]);  axs.set_xticks(range(0,N+1))
#        axs=f.add_subplot(111); plt.semilogy(M,Rlim,label='%.0e' % Decimal(tol)); axs.set_xlim(0.5, N); axs.set_ylim(1, 100000); axs.set_ylabel('$R_{lim}=(1-\omega_{lim})^{-1}$'); axs.set_xlabel('M');
#        axs.set_xticks(range(1,N+1)); axs.tick_params(axis='major', which='major', labelsize=8); axs.grid(which='major', axis='both',ls='dashed')
#    axs.legend(loc='upper left', shadow=False, fontsize=10, title='Tolerance')
#    f.savefig('equiv_00-10.tiff', dpi=300)
#    f.clf()       
#    
#    1.-Pr_lessM > 1.-1./Z ; Pr_lessM < 1./Z             #  ---> xx10xx equivalent to xx11xx. Start D and enter Nc>=M by mistakes
#    print(1.-Pr_lessM)  # never commit enough mitakes, unless w=1
#
############################################




#    beta=1.
#    Z=100
#    N=9 21     1   0     1   0   1   0     9.24e-01   0.50     0.00

#    w=0.9
#    eps=0.01
#    
#    c1=5.
#    csV=c1*np.linspace(0,0.3,51)
#    lambV=np.linspace(0,1,51)
#    #bV=np.array([20.,10.,5.])
#    bV=c1*np.array([20.])
#    MV= np.array([5,6,7,8,9]) #np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
##    beta=1.
##    Z=100
##    N=9
##    w=1.
##    eps=0.01
##    c1=5.
##    csV=c1*np.linspace(0,0.3,51)
##    lambV=np.linspace(0,1,51)
##    #bV=np.array([20.,10.,5.])
##    bV=c1*np.array([20.])
##    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
##    QV= np.array([1.0,3.0,4.5,5.0,7.0,9.0]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
##
#    bigmatSD=np.load('file_SD_N9_beta5_b20_w1_X_testQM.npy')
#    vmin=-0.
#    vmax=1.
#    groups=[20,28,33,35,48, 49, 50, 51,52, 53, 54, 55, 60, 61, 62, 63]
#    nogroups= [x for x in range(0,64) if x not in groups]
#    print(groups)
#    print(nogroups)
#    plot_SDcslambDIF('all--','others','groups',nogroups,groups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)

    
    
#    np.set_printoptions(threshold=np.inf)
#    fixM,SD=doONLYONE()
#    Mred,SDred=groupM(fixM,SD[:,0],gST)
#    print(Mred); print(np.sum(SDred))
    
########################## find drift groups ###############################3 
#    beta=1.
#    Z=100
#    N=9
#    M=1
#    Q=4.5
#    lamb=0.25 #0.5 #0.8
#    eps=0.01
#    w=0.9
#    #H, L
#    c1=0.5 #2.5
#    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#    cs=np.array([1., 1.]) *7. *c1  #*0.06 *c  *0.8
#    b=np.array([1., 0.]) *20. *c1   #7*c  
#    
#    STRmPUR=declareSTR(0)
#    expb=np.exp(-beta)
#    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
#    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
#    fixMvec=readfixMvec(labelfile)
#    fixM=calcFIXM(coef,expb,Z,fixMvec)
#    
#    doINI(N,Z,M,Q,eps,w)
#    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)
#    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
#    SSDsort=SSD[np.argsort(SSD[..., 8])] 
#    for i in range(0,len(SSDsort)):
#       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
#    
#    groups=findDRIFTgroup(fixM,100)
#    print(groups)
################################################  
   
   
###########   CREATE BARS  ############################################################################    
#    beta=1.
#    Z=100
#    N=9
#    w=0.9
#    eps=0.01
#
#    c1=0.5
#
###------ vertical ---------------------    
##    csv=np.array([0, 0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6])
##    cs1V=   c1*np.transpose(np.array([csv,csv]))
##    lambV=  np.transpose(np.array([[0.3]*len(cs1V),[0.7]*len(cs1V)]))
##    bV=c1*np.transpose(np.array([[20]*len(lambV)]*2))
##    epsV= np.transpose(np.array([[0.01]*len(lambV)]*2))  
##    QV=    np.transpose( np.array([[4.5]*len(lambV)]*2))
##    MV=    np.transpose( np.array([[3]*len(lambV)]*2))
###--------------------------------------  
##    
#####------ horizontal ---------------------    
#    lv=[0,0.2,0.4,0.5,0.6,0.8,1.]
#    lambV=   np.transpose(np.array([lv,lv]))
#    cs1V=  c1*np.transpose(np.array([[0.01]*len(lambV),[0.5]*len(lambV)])) 
#    bV=c1*np.transpose(np.array([[15]*len(lambV)]*2))
#    epsV= np.transpose(np.array([[0.01]*len(lambV)]*2))  
#    QV=    np.transpose( np.array([[4.5]*len(lambV)]*2))
#    MV=    np.transpose( np.array([[3]*len(lambV)]*2))
####-------------------------------------- 
#
####------ horver ---------------------    
##    lv=[0,0.2,0.4,0.5,0.6,0.7,0.8,1.] # 2 horizontal
##    csv=np.array([0, 0.2,0.4,0.6,0.8,1.,1.2,1.5]) # 1 vertical
##    lambV=   np.transpose(np.array([lv,lv,[0.5]*len(lv)]))
##    cs1V=  c1*np.transpose(np.array([[0.5]*len(lv),[1.5]*len(lv),csv])) 
##    
##    
##    bV=c1*np.transpose(np.array([[20]*len(lambV)]*3))
##    epsV= np.transpose(np.array([[0.01]*len(lambV)]*3))  
##    QV=    np.transpose( np.array([[4.5]*len(lambV)]*3))
##    MV=    np.transpose( np.array([[3]*len(lambV)]*3))
###-------------------------------------- 
#    
##  
###    bV=np.array([20,20,20,20,20,20]) *c1
###    epsV= np.array([0.01,  0.01,  0.01 ,  0.01,       0.01,0.01])
###    lambV=  np.array([0.5,  0.5,  0.5 ,  0.5,       0.5,0.5])
###    cs1Vl=   np.array([0.05, 0.25, 0.05,   0.25,    0.1,0.1])
###    cs1V= cs1Vl*c1       
###    QV=     np.array([4.5,  4.5,  7,      7,        4.5,4.5])
###    MV=     np.array([5,     5,    9,     9,        5,  7])     
##    
##    #lambV=  np.array([0.5,  0.5,  0.5 ,  0.5,       0.5,0.5])
##    #cs1Vl=   np.array([0.05, 0.13, 0.135,   0.14,    0.145,0.15])
##    #cs1V= cs1Vl*c1       
##    #QV=     np.array([4.5,  4.5,  4.5,     4.5,        4.5,4.5])
##    #MV=     np.array([5,     5,    5,     5,        5,  5])   
##   
##    #epsV= np.array([0.01,  0.1,  0.3 ,  0.01,       0.1,0.3])
##    #lambV=  np.array([0.5,  0.5,  0.5 ,  0.5,       0.5,0.5])
##    #cs1Vl=   np.array([0.05, 0.05, 0.05,   0.25,    0.25,0.25])
##    #cs1V= cs1Vl*c1       
##    #QV=     np.array([4.5,  4.5,  4.5,      4.5,        4.5,4.5])
##    #MV=     np.array([5,     5,    5,     5,        5,  5])   
##    
###    lambV=  np.array([0.55,  0.55,  0.65 ,  0.6,       0.65,0.6])
###    cs1Vl=   np.array([0.1, 0.1, 0.12,   0.12,    0.15,0.25])
###    cs1V= cs1Vl*c1       
###    QV=     np.array([4.5,  4.5,  5.,      5.,        5.,5.])
###    MV=     np.array([7,     7,    7,     7,        7,  7])   
#    
#    STv=np.array([20, 22, 28, 30, 33, 35, 43, 51, 53, 55, 57, 58, 59, 61, 63])
#    STv=np.array(range(0,64))
#    STv=np.array([20,28, # SC
#                  22,30, # SC mut
#                  33,35, # SD
#                  41,43, # SD mut
#                  48,49,50,51, # SF
#                  56,57,58,59, # SF mut
#                  52, 53, 54, 55, # MB C
#                  60, 61, 62, 63, # MB D
#                  0, 4, 8, 12,
#                  1, 3, 5, 7, 9, 11, 13, 15,
#                  2, 6, 10, 14,
#                  17, 19,
#                  21, 23, 29, 31,
#                  25, 27,
#                  36, 44,
#                  37, 39, 45, 47,
#                  38, 46,
#                  16, 18, 24, 26, 32, 34, 40, 45 ])  # not in groups
#    STvC=(['xkcd:dark red','xkcd:dark red', # SC
#          'xkcd:salmon','xkcd:salmon', # SC mut
#          'xkcd:dark grey','xkcd:dark grey', # SD
#          'xkcd:grey','xkcd:grey', # SD mut 21     1   0     1   0   1   0     9.24e-01   0.50     0.00
#          #'xkcd:sienna','xkcd:sienna', # SD
#          #'xkcd:tan','xkcd:tan', # SD mut 21     1   0     1   0   1   0     9.24e-01   0.50     0.00
#          #'xkcd:dark green','xkcd:dark green','xkcd:dark green','xkcd:dark green', # SF
#          'xkcd:violet','xkcd:violet','xkcd:violet','xkcd:violet', # SF
#          'xkcd:green','xkcd:green','xkcd:green','xkcd:green', # SF mut
#          'xkcd:medium blue','xkcd:medium blue','xkcd:medium blue','xkcd:medium blue', # MB C
#          'xkcd:dark blue', 'xkcd:dark blue', 'xkcd:dark blue', 'xkcd:dark blue'] # MB D
#          + ['xkcd:tan']*(64-24) ) # all the others
#    
#    STRp=declareSTR(0)
#    labup=['']*len(STv)
#    for i in range(0,len(STv)):
#        sen=str(STRp[STv[i],:])
#        labup[i]=sen.replace(". ","").replace("[","").replace(".]","").replace(" ","")
#    
#    label='BAR_w09_horit_Q4.5_M3_r15' #'BAR_ttt'
#    vmax=0.4
#    ncol=2
#    import matplotlib.pyplot as plt
#    f = plt.figure()
#    j=0
#    for i in range(0,len(cs1V)):
#        nrow=len(cs1V)
#        npl=ncol*i+j+1  
#        axs = f.add_subplot(nrow,ncol,npl)
#        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
#        if i+1!=len(cs1V):
#            labx=[]
#        else:
#            labx=labup
#        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
#    j=1    
#    for i in range(0,len(cs1V)):
#        nrow=len(cs1V)
#        npl=ncol*i+j+1     
#        axs = f.add_subplot(nrow,ncol,npl)
#        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
#        if i+1!=len(cs1V):
#            labx=[]
#        else:
#            labx=labup
#        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
#        axs.yaxis.set_major_locator(plt.NullLocator())
##    j=2    
##    for i in range(0,len(cs1V)):
##        nrow=len(cs1V)
##        npl=ncol*i+j+1     
##        axs = f.add_subplot(nrow,ncol,npl)
##        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
##        if i+1!=len(cs1V):
##            labx=[]
##        else:
##            labx=labup
##        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
##        axs.yaxis.tick_right()
#    f.text(0.5, 0.96, 'Stationary Distribution', rotation=0, va='center', ha='center',size=8)
#    #f.text(0.05, 0.5, 'Stationary Distribution', rotation=90, va='center', ha='center',size=8)
#    f.savefig(label+'.eps', dpi=300)
#    f.clf()
#    
#
##------- only chosen strategies -------------------------------------------------------
##    STv=np.array([28,20,29, # SC
##                  #16,17,18,19,21,22,23,24,25,26,27,30,31, # 10 non SC
##                  35, 33, 39, # SD
##                  #32,34,36,37,38,40,41,42,43,44,45,46,47, # 01 non SD
##                  48,49,50,51, # FR-C
##                  56,57,58,59, # FR-01
##                  52, 53, 54, 55, # FR-10
##                  60, 61, 62, 63, # FR-D
##                  ])  
##    STvC=(['xkcd:dark red']*3  # SC
##          #+['xkcd:salmon']*(16-3) # 10 non SC
##          +['xkcd:dark grey']*3  # SC
##          #+['xkcd:grey']*(16-3) # 10 non SC
##          +['xkcd:violet']*4  # FR-C
##          +['xkcd:green']*4  # FR-01
##          +['xkcd:medium blue']*4  # FR-10
##          +['xkcd:dark blue']*4)  # FR-D          
##
##            
##    STRp=declareSTR(0)
##    labup=['']*len(STv)
##    for i in range(0,len(STv)):
##        sen=str(STRp[STv[i],:])
##        labup[i]=sen.replace(". ","").replace("[","").replace(".]","").replace(" ","")
##    
##    label='BAR_w09_horver_Q4.5_M3' #'BAR_ttt'
##    vmax=0.5
##    ncol=3
##    import matplotlib.pyplot as plt
##    f = plt.figure()
##    j=0
##    for i in range(0,len(cs1V)):
##        nrow=len(cs1V)
##        npl=ncol*i+j+1  
##        axs = f.add_subplot(nrow,ncol,npl)
##        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
##        if i+1!=len(cs1V):
##            labx=[]
##        else:
##            labx=labup
##        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
##    j=1    
##    for i in range(0,len(cs1V)):
##        nrow=len(cs1V)
##        npl=ncol*i+j+1     
##        axs = f.add_subplot(nrow,ncol,npl)
##        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
##        if i+1!=len(cs1V):
##            labx=[]
##        else:
##            labx=labup
##        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
##        axs.yaxis.set_major_locator(plt.NullLocator())
##    j=2    
##    for i in range(0,len(cs1V)):
##        nrow=len(cs1V)
##        npl=ncol*i+j+1     
##        axs = f.add_subplot(nrow,ncol,npl)
##        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
##        if i+1!=len(cs1V):
##            labx=[]
##        else:
##            labx=labup
##        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
##        axs.yaxis.tick_right()
##    f.text(0.5, 0.96, 'Stationary Distribution', rotation=0, va='center', ha='center',size=8)
##    #f.text(0.05, 0.5, 'Stationary Distribution', rotation=90, va='center', ha='center',size=8)
##    f.savefig(label+'.svg', dpi=300)
##    f.clf()
##--------------------------------------------------
#    
################################################################################################################    
   
    
  
    

#    doINI(N,Z,M,Q,eps)
    
    
    
#    M=4; N=5; eps=0.01; STATEmat,nSTATEv=declareSTATE(N);
#    Z=100; H=calcH(N,Z)
    
#    STRm=declareSTR(eps)
#    nSTR=STRm.shape[0]; 
#    for i in range(0,nSTR):
#        for j in range(i+1,nSTR):
#            SGi=STRm[i,0:2]; ACTi=STRm[i,2:6]
#            SGj=STRm[j,0:2]; ACTj=STRm[j,2:6]

            #for k in range(1,N):
            #    print([i, j, k])
            #    Q=2.5; BCiA,BCjA=calcBC2st(SGi,ACTi,SGj,ACTj,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
            #    Q=3.; BCiB,BCjB=calcBC2st(SGi,ACTi,SGj,ACTj,k,N,M,Q,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
            #    if not (np.all(np.abs(BCiA-BCiB)<0.000001) and np.all(np.abs(BCjA-BCjB)<0.000001)):
            #        print(BCiA)
            #        print(BCiB)
                    
#            print([i, j])
#            Q=2.5; stermA,stermIA=calcFIX1vec(i,j,STRm,N,Z,M,Q,STATEmat,nSTATEv,H)
#            Q=3.; stermB,stermIB=calcFIX1vec(i,j,STRm,N,Z,M,Q,STATEmat,nSTATEv,H)
#            if not (np.all(np.abs(stermA-stermB)<0.000001) and np.all(np.abs(stermIA-stermIB)<0.000001)):
#                print(stermA,stermIA)
#                print(stermB,stermIB)
                
#    fixM=doONLYONE()
   
    #beta=1.
    #Z=100
    #N=5
    #M=4
    #Q=2.5
    #lamb=0.8
    #eps=0.01
    #H=calcH(N,Z)
    #STRm=declareSTR(eps); # nSTR=STRm.shape[0]; #print(STm[0],STm[63])
    #STATEmat,nSTATEv=declareSTATE(N); #print(STATEmat); print(nSTATEv)
    #fixMvec=calcFIXMvec(N,Z,M,Q,STRm,STATEmat,nSTATEv,H)
    

    
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#    
##    csV=np.linspace(0,1,51)
##    lambV=np.linspace(0,1,11)
##    #bV=np.array([20.,10.,5.])
##    bV=np.array([10.])
##    MV=np.array([1,2,3,4,5])
##    QV=np.array([1,2,2.5,3,4,5])
##    #MV=np.array([1,3,5,7,9])
##    #QV=np.array([1,3,4.5,5,7,9])

############## SAVE BIG MATRIX for Q-M plots ################### 
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01
#
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#
#    lambV=np.linspace(0,1,31)
#    
#    csVo= np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18])    #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5])    #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
# 
#    ##### test ####
##    lambV=np.linspace(0,1,31)  
##    csVo= np.linspace(0,2,31) #np.linspace(0,0.3,51)
##    bVo=np.array([20.])
##    MV= np.array([5]) #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
##    QV= np.array([4.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
##    w=0.9
##    c1=0.2
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    #bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    #np.save('file_SD_N9_beta02_b20_w09_NEWtest_4.5_5.npy',bigmatSD)
#    ################
#    
#    
##    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
##    
#    bVo=np.array([20.])
##    w=0.9
##    c1=0.5
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N9_beta05_b20_w09_bper.npy',bigmatSD)
###    np.save('file_COOP_N9_beta5_b10_w1_X.npy',bigmatCOOP)
##    w=0.9
##    c1=0.5
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N9_beta1_b20_w09_NEWtest.npy',bigmatSD)
#
##    bVo=np.array([30.])
##    w=0.9
##    c1=0.5
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N18_beta05_b30_w09_NEW.npy',bigmatSD)
##    
#    w=1.
#    c1=0.5
#    csV=c1* csVo #np.linspace(0,0.3,51)
#    bV=c1* bVo
#    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
#    np.save('file_SD_N9_beta05_b20_w1_bper.npy',bigmatSD)
#
#
###################################################################    
  
       

        
        
#    STs00,STs11,STs10,STs01,STsign, STsignonly, STmem, STmemonly, STsignmem =classST()
#    print(STsign)  
#    print(STsignonly) 
#    print(STmem) 
#    print(STmemonly)  
#    print(STsignmem) 
    
#    csV=np.linspace(0,1,51)
#    lambV=np.linspace(0,1,11)
#    bV=np.array([20.])
#    MV=np.array([1,2,3,4,5])
#    QV=np.array([1,2,2.5,3,4,5])

###################  PLOTS Q-M  ############################################################################    
##    bigmatCOOP=np.load('file_COOP_N9_beta5_b20_w1_X.npy')
##    lab='SD_N9_beta05_b20_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW.npy')
##    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD)
# 
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([20.])
#    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #np.array([1,3,5,7,9]) #np.array([5,6,7,8,9])   #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=0.5
#    csV= csVo
#    bV=c1* bVo
###
##    lab='SD_N9_beta05_b20_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW.npy')    
##    ext='tiff'
##    vmin=0.4; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF,gSFm,gSC,gSD]
##    #ngroups=['Blues','Purples','Greens','Reds','Greys']
##    ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011]',r'[01 1100] + [01 1110]']
##    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##    vmin=0; vmax=1.
##    plot_SDcslambDIF_agre('SD_gS00_'+lab,[gS00],['Blues'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##    plot_SDcslambDIF_agre('SD_gS10_'+lab,[gS10],['Reds'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##    plot_SDcslambDIF_agre('SD_gS01_'+lab,[gS01],['Greys'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##
#    
#    lab='SD_N9_beta05_b20_w09'
#    bigmatSD=np.load('file_'+lab+'_bper.npy')    
#    ext='eps'
#    vmin=0.4; vmax=1.
#    comaps=['Blues','Oranges','Purples','Reds','Greys']
#    groups=[gMB,gMBc,gSF,gSCt,gSDt]
#    #ngroups=['Blues','Purples','Greens','Reds','Greys']
#    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
#    ngroups=[r'FR-D + FR-10$^{M\geqslant3}$', r'FR-10$^{M<3}$',r'FR-C',r'SC + SC$_{C}$',r'SD + SD$_{C}$']
#
#    plot_SDcslambDIF_agre('SD_agre_'+lab+'_bper',groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext,1)
#
#    lab='SD_N9_beta05_b20_w1'
#    bigmatSD=np.load('file_'+lab+'_bper.npy')    
#    ext='eps'
#    vmin=0.4; vmax=1.
#    #comaps=['Blues','Purples','Greens','Reds','Greys']
#    #groups=[gMB,gSF+gMBc,gSFm,gSCt,gSDt]
#    comaps=['Blues','Oranges','Purples','Reds','Greys']
#    groups=[gMB,gMBc,gSF,gSCt,gSDt]
#    ##ngroups=['Blues','Purples','Greens','Reds','Greys']
#    ##ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
#    ngroups=[r'FR-D + FR-10$^{M>5}$', r'FR-10$^{M\leqslant5}$',r'FR-C',r'SC + SC$_{C}$',r'SD + SD$_{C}$']
#    #ngroups=[r'FR-D + FR-10$^{M>5}$', r'FR-C + FR-10$^{M<5}$',r'FR-01',r'SC + SC$_{C}$$^{M>5}$ + SC$_{D}$$^{M<5}$',r'SD + SD$_{C}$$^{M>5}$ + SD$_{D}$$^{M<5}$']
#
#    plot_SDcslambDIF_agre('SD_agre_'+lab+'_bper',groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext,5)
#
##
##
##
##
##    
##    vmin=0.1
##    vmax=1.
##    comap='Reds'
##    plot_SDcslambDIF_1('gSC','gSC',' ',gSC,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gSD','gSD',' ',gSD,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gSF','gSF',' ',gSF,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gSCm','gSCm',' ',gSCm,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gSDm','gSDm',' ',gSDm,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gSFm','gSFm',' ',gSFm,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gMB','gMB',' ',gMBc+gMBd,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('ALL-g','All-groups','groups',gALL,gMBc+gMBd+gSC+gSD+gSF+gSCm+gSDm+gSFm,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    #plot_SDcslambDIF_1('gMBc','gMBc',' ',gMBc,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    #plot_SDcslambDIF_1('gMBd','gMBd',' ',gMBd,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    
##    vmin=0
##    vmax=1.
##    comap='Reds'
##    plot_SDcslambDIF_1('gS00','gS00',' ',gS00,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gS11','gS11',' ',gS11,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gS10','gS10',' ',gS10,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    plot_SDcslambDIF_1('gS01','gS01',' ',gS01,[],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,comap)
##    
#################################################################################################################


############ PAYOFSS ####################    
#    STRp=declareSTR(0)
#    for i in range(0,STRp.shape[0]):
#        print([i, STRp[i,:]])
#   
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01
#    w=0.9
#    c1=0.5
#
#    lambV=np.linspace(0,1,31)
#    csVo= np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([20.])
#    MV= np.array([1,3,5,7,9]) #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5])
#    csV= c1*csVo
#    bV=c1* bVo
#    bigmatSD=np.load('file_SD_N9_beta05_b20_w09_NEW.npy')  # CAREFUL: IT HAS TO BE COHERENT WITH PARAMETERS DEFINED ABOVE
#    bigmatPAY=calcBIGPAY(bigmatSD,csV,lambV,MV,QV,bV[0],c1,N,eps,w)
#    np.save('file_PAY_N9_beta05_b20_w09_NEW.npy',bigmatPAY)
#    
#    bigmatPAY=np.load('file_PAY_N9_beta05_b20_w09_NEW.npy')
#    #print(np.amin(bigmatPAY)/bV[0])
#    vmin=0.
#    vmax=1.    
#    plot_PAYcslamb('PAYavg_w09',bigmatPAY,csVo,lambV,bV,MV,QV,vmin,vmax)
############################################# 



'''
    
    
    plot_COOPcslamb(bigmatCOOP,csV,lambV,bV,MV,QV)
    vmin=0; vmax=1
    #plot_SDcslamb('todo',np.array(np.linspace(0,63,64),int),bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal 00',STs00,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal 11',STs11,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal 10',STs10,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal 01',STs01,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal',STsign,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Memory',STmem,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal only',STsignonly,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Memory only',STmemonly,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    plot_SDcslamb('Signal & Memory',STsignmem,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    vmin=0; vmax=0.3
    for i in range(0,STRp.shape[0]):
        plot_SDcslamb(str(i)+'  '+str(STRp[i,:]),[i],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    #i=63; plot_SDcslamb(str(i)+'  '+str(STRp[i,:]),[i],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    #i=28; plot_SDcslamb(str(i)+'  '+str(STRp[i,:]),[i],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    #i=35; plot_SDcslamb(str(i)+'  '+str(STRp[i,:]),[i],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    #i=51; plot_SDcslamb(str(i)+'  '+str(STRp[i,:]),[i],bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax)
    
    
  

    
    csV=np.linspace(0,1,51)
    lambV=np.linspace(0,1,11)
    bV=np.array([20.])
    MV=np.array([1,2,3,4,5])
    QV=np.array([1,2,2.5,3,4,5])
    
    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
    bigmatCOOP,bigmatSD=doMATSD(nSTR,csV,lambV,bV,MV,QV)
    #bigmatCOOP=doMATCOOP(csV,lambV,bV,MV,QV)
    np.save('file_SD_1.npy',bigmatSD)
    
    ##bigmatCOOP=np.zeros((len(csV),len(lambV),len(bV),len(MV),len(QV)))
    ##bigmatCOOP=np.load('file_1.npy')
    #print(bigmatCOOP)
    #plot_COOPcs(bigmatCOOP,csV)
    
    ##plot_COOPcslamb(bigmatCOOP,csV,lambV,bV,MV,QV)
    ##plot_COOPcslamb(bigmatCOOP,csV,lambV,1)
    ##plot_COOPcslamb(bigmatCOOP,csV,lambV,2)
    ##plot_COOPcslamb(bigmatCOOP,csV,lambV,3)
    #plot_COOPcslamb(bigmatCOOP,csV,lambV,4)
    
    


    beta=1.
    Z=100
    N=5
    M=5
    Q=3
    lamb=0.5
    #H, L
    c=np.array([1., 1.])*1.
    cs=np.array([1., 1.])*0.
    b=np.array([20., 0.]) #*c    
 
    eps=0.01

    STRmPUR=declareSTR(0.)
    
    
    #STi=0; STj=63; print(STRm[STi]); print(STRm[STj])
    #t0=time.clock()
    #pfix,pfixI=calcFIX1(expb,b,c,cs,lamb,STi,STj,STRm,N,Z,M,STATEmat,nSTATEv,H); print([pfix,pfixI])
    #t1=time.clock(); print(t1-t0)
    #total_time = timeit.timeit('calcFIX1(expb,b,c,cs,lamb,STi,STj,STRm,N,Z,M,STATEmat,nSTATEv,H)'); print(total_time)
          
    ##t0=time.clock()                    
    ##fixM=calcFIXM(b,c,cs,lamb,beta,N,Z,M,STRm,STATEmat,H)
    ##t1=time.clock(); print(t1-t0)
    ##SD=calcSD(fixM)
    ##t2=time.clock(); print(t2-t1)
    
    #labelfile='N_'+str(N)+'_M_'+str(M)
    #fixMvec=calcFIXMvec(N,Z,M,STRm,STATEmat,H)
    #writefixMvec(fixMvec,labelfile)
    doINI(N,Z,M,Q,eps)
    #coef=np.exp(-beta*np.array([[b[0]*lamb, c[0]*lamb, c[1]*lamb],[b[0]*(1.-lamb), c[0]*(1.-lamb), c[1]*(1.-lamb)]]))
    ###coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    #print(fixMvec)
    ###fixM=calcFIXM(coef,expb,Z,fixMvec)
    ###SD=calcSD(fixM)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps)
    ##SD=np.zeros(64) 
   
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD)

    
    #for i in range(0,len(SD)):
    #    print([STRmPUR[i,:], SD[i]*len(SD)])
    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD*len(SD)),axis=1)
    #print(SSD)
    SSDsort=SSD[np.argsort(SSD[..., 8])] 
    #print(SSDsort)
    for i in range(0,len(SSDsort)):
       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
    print(COOPtot)
    
    #listERS=calcERS(b,c,cs,lamb,N,Z,M,eps)
    #print('----')
    #for i in range(0,len(listERS)):
    #   print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f}'.format(listERS[i],STRmPUR[listERS[i],0],STRmPUR[listERS[i],1],STRmPUR[listERS[i],2],STRmPUR[listERS[i],3],STRmPUR[listERS[i],4],STRmPUR[listERS[i],5]))
    
    #SIGi=[1.-eps,1.-eps,1.-eps,1.-eps]; ACTi=[0.+eps,1.-eps]; SIGj=[0.+eps,0.+eps,0.+eps,0.+eps]; ACTj=[0.+eps,1.-eps];
    #print('------------------')
    #k=1; BCi,BCj=calcBC2st(SIGi,ACTi,SIGj,ACTj,k,N,M,STATEmat[k,0:nSTATEv[k],0:nSTATEv[k]])
    #print(BCi)
    #print(BCj)

    
    #print(SD)
    #print(coop)

'''

