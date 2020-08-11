# -*- coding: utf-8 -*-

from functions import *

######################## MAIN SCRIPTS #################################################################
# to be modified and run changing parameters

if __name__ == "__main__":
    import numpy as np
    import time; import timeit

    #doONLYONE()

#---- group of strategies --------------------------
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
   
    gSCO=[22,30]
    gSDO=[43,41]
       
    gSC1=[29,28]
    gSD1=[39,35]
    
    gSCt=[20,29,28]
    gSDt=[33,39,35]
    gGt=gSCt+gSDt+gSF+gSFm+gMBc+gMBd
    gNOt=[x for x in gALL if x not in gGt]

    gNO2 = [x for x in gALL if x not in gSC+gSD+gSF+gSFm+gMBc+gMBd ]    

    gALL=list(range(0,64))
#----------------------------------------------------------


############### Cooperation level ##########################
#
#    beta=1.
#    Z=100
#    N=9
#    M=5
#    Q=4.5 #4.5
#    lamb=0.5 #0.5 #0.8
#    eps=0.01
#    w=1.
#    #H, L
#    c1=1. #2.5
#    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#    b=np.array([1., 0.]) *10. *c1   #7*c  
#    
#    #csVo= np.linspace(0,2,10) 
#    lambV= np.linspace(0,1,50)
#    
#    expb=np.exp(-beta)
#
#    coop=np.zeros((len(lambV),10,3))
#    coopPGG=np.zeros((len(lambV),10))
#    for i in range(0,len(lambV)):
#        lamb=lambV[i]
#        
#        cs=np.array([1., 1.]) *0. *c1 
#        coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])       
#
#        doINI(N,Z,M,Q,eps,w)
#        SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,0,:]=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
#        doINI_SIG(N,Z,M,Q,eps,w)
#        SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,1,:]=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
#        cs=np.array([1., 1.]) *0.1 *c1 
#        coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])       
#
#        doINI(N,Z,M,Q,eps,w)
#        SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,2,:]=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
#        doINI_SIG(N,Z,M,Q,eps,w)
#        SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,3,:]=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#
#        cs=np.array([1., 1.]) *0.3 *c1 
#        coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])       
#
#        doINI(N,Z,M,Q,eps,w)
#        SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,4,:]=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
#        doINI_SIG(N,Z,M,Q,eps,w)
#        SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,5,:]=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]      
#        
#        cs=np.array([1., 1.]) *0.5 *c1 
#        coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])       
#
#        doINI(N,Z,M,Q,eps,w)
#        SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,6,:]=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
#        doINI_SIG(N,Z,M,Q,eps,w)
#        SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,7,:]=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]  
#
#        
#        doINI_REC(N,Z,M,Q,eps,w)
#        SD=doREST_REC(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,8,:]=doHOMO_REC(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#
#        doINI_CD(N,Z,M,Q,eps,w)
#        SD=doREST_CD(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#        coop[i,9,:]=doHOMO_CD(lamb,eps,N,M,Q,b,c,cs,SD,w)[2]
#        
##        PAYhomo,COOPhomo,COOPtot=doHOMO_REC(lamb,eps,N,M,Q,b,c,cs,SD,w)
##        print(COOPhomo[:,0])
##        SSD=np.concatenate((declareSTR_REC(0),np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
##        SSDsort=SSD[np.argsort(SSD[..., 8])] 
##        for ii in range(0,len(SSDsort)):
##            print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
##        print(i,lambV[i],COOPtot)
#        
#        print(i,lambV[i],coop[i,:,0])
#
#    np.save('coop_w1_beta1_r10',coop)
#    #np.save('coop_w1_beta1_r10_M7',coop)
#    #np.save('coop_w1_beta1_r10_M9',coop)
##
##    coop=np.load('coop_w1_beta1_r10'+'.npy')      
##    coop7=np.load('coop_w1_beta1_r10_M7'+'.npy') 
##    coop9=np.load('coop_w1_beta1_r10_M9'+'.npy') 
##
##    import matplotlib.pyplot as plt        
##    #lab=["S+R, $c_S=0$", "S, $c_S=0$", "S+R, $c_S=0.5$", "S, $c_S=0.5$", "S+R, $c_S=1$", "S, $c_S=1$","S+R, $c_S=1.5$", "S, $c_S=1.5$", "R", "C+D"]
##    lab=["S+R", "S", "S+R", "S", "R", "C+D","S+R", "S","S+R", "S"] #["S+R", "S", "S+R", "S","S+R", "S","S+R", "S", "R", "C+D"]
##    lin=['b-'          ,    'g-',           'b--',          'g--',         'b:',          'g:',          'b-.',          'g-.',     'r-',   'k-']
####    f = plt.figure()
####    for j in range(0,len(coop[0,:,0])):
####        axs=f.add_subplot(111); plt.plot(lambV,coop[:,j,0],lin[j],label=lab[j]); axs.set_xlim(0., 1.); axs.set_ylim(0., 1.); axs.set_ylabel('Level of cooperation'); axs.set_xlabel('$\lambda$');
####        #axs.set_xticks(range(1,N+1)); axs.tick_params(axis='major', which='major', labelsize=8); axs.grid(which='major', axis='both',ls='dashed')
####    axs.legend(loc='best', shadow=False, fontsize=8)
####    f.savefig('mechanisms_w09.eps', dpi=300)
####    f.clf()         
####    f = plt.figure()
####    for j in range(0,len(coop[0,:])):
####        axs=f.add_subplot(111); plt.plot(lambV,coop[:,j,1],lin[j],label=lab[j]); axs.set_xlim(0., 1.); axs.set_ylim(0., 1.); axs.set_ylabel('Level of cooperation'); axs.set_xlabel('$\lambda$');
####        #axs.set_xticks(range(1,N+1)); axs.tick_params(axis='major', which='major', labelsize=8); axs.grid(which='major', axis='both',ls='dashed')
####    axs.legend(loc='best', shadow=False, fontsize=8)
####    f.savefig('mechanisms_PGG_w09.eps', dpi=300)
####    f.clf() 
###
##
###    f = plt.figure()
###    
###    ax=plt.subplot(121)
###    for j in range(0,len(coop[0,:,0])):
###        plt.plot(lambV,coop[:,j,0],lin[j],label=lab[j]); plt.xlim(0., 1.); plt.ylim(0., 1.); plt.ylabel('Level of cooperation'); plt.xlabel('$\lambda$');
###    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(["0","0.25","0.5","0.75","1"]); ax.tick_params(axis='major', which='major', labelsize=8); ax.grid(which='major', axis='both',ls='dashed')
###    plt.title('$G$ + $\hat{G}$')
###    h, l = ax.get_legend_handles_labels()
###    ph = plt.plot([],marker="", ls="")[0]
###    handles = [ph,h[0],h[1],ph,h[2],h[3],ph,h[4],h[5],ph,h[6],h[7],ph,h[8],h[9] ]
###    labels = ["$c_S=0$",lab[0],lab[1],"$c_S=0.1$",lab[2],lab[3],"$c_S=0.3$",lab[4],lab[5],"$c_S=0.5$",lab[6],lab[7]," ",lab[8],lab[9] ]
###    leg=plt.legend(handles, labels, bbox_to_anchor=(0., 1.15, 2.2, .102), loc=8,
###           ncol=5, mode="expand", borderaxespad=0.,fontsize=8,edgecolor='black')
###    for t in leg._legend_handle_box.get_children():
###        for hpack in t.get_children()[0:1]:
###            hpack.get_children()[0].set_width(0)
###
###            
###    for j in range(0,len(coop[0,:])):
###        axs=f.add_subplot(122); plt.plot(lambV,coop[:,j,1],lin[j],label=lab[j]); axs.set_xlim(0., 1.); axs.set_ylim(0., 1.);  axs.set_xlabel('$\lambda$'); #axs.set_ylabel('Level of cooperation');
###        #axs.set_xticks(range(1,N+1)); axs.tick_params(axis='major', which='major', labelsize=8); axs.grid(which='major', axis='both',ls='dashed')
###    axs.set_xticks([0,0.25,0.5,0.75,1]); axs.set_xticklabels(["0","0.25","0.5","0.75","1"]); axs.tick_params(axis='major', which='major', labelsize=8); axs.grid(which='major', axis='both',ls='dashed')
###    plt.title('$G$')
###    plt.subplots_adjust(top=0.7)
###    #axs.legend(loc='best', shadow=False, fontsize=8)
###    f.savefig('coop_mechanisms_w1_r10.eps', dpi=300)
###    f.clf() 
###    
##    
###    
##    f,axs=plt.subplots(nrows=2, ncols=2, sharex='all', sharey='all' )
##
##    print(coop)
##   
##    selcoop=np.array([0,1,4,5,8,9])
##    
##    ax=ax00=axs[0,0]
###    for j in range(0,len(coop[0,:,0])):
##    for jj in range(0,len(selcoop)):
##        j=selcoop[jj]
##        ax.plot(lambV,coop[:,j,0],lin[j],label=lab[j]); ax.set_xlim(0., 1.); ax.set_ylim(0., 1.); #ax.ylabel('Level of cooperation'); ax.xlabel('$\lambda$');
##    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(["0","0.25","0.5","0.75","1"]); ax.tick_params(axis='major', which='major', labelsize=8); ax.grid(which='major', axis='both',ls='dashed')
##    ax.set_title('$G$ + $\hat{G}$',size=10)
##
##    ax=axs[0,1]          
###    for j in range(0,len(coop[0,:,0])):
##    for jj in range(0,len(selcoop)):
##        j=selcoop[jj]
##        ax.plot(lambV,coop[:,j,1],lin[j],label=lab[j]); ax.set_xlim(0., 1.); ax.set_ylim(0., 1.) #;  ax.set_xlabel('$\lambda$'); #axs.set_ylabel('Level of cooperation');
##    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(["0","0.25","0.5","0.75","1"]); ax.tick_params(axis='major', which='major', labelsize=8); ax.grid(which='major', axis='both',ls='dashed')
##    ax.set_title('$G$',size=10)
##    ax.text(1.1,0.48,"$M=5$", size=12 ,va='center')
##
##    #lambV= np.linspace(0,1,20)
##    ax=axs[1,0]          
###    for j in range(0,len(coop[0,:,0])):
##    for jj in range(0,len(selcoop)):
##        j=selcoop[jj]
##        ax.plot(lambV,coop7[:,j,0],lin[j],label=lab[j]); ax.set_xlim(0., 1.); ax.set_ylim(0., 1.);  ax.set_xlabel('$\lambda$',size=12); #axs.set_ylabel('Level of cooperation');
##    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(["0","0.25","0.5","0.75","1"]); ax.tick_params(axis='major', which='major', labelsize=8); ax.grid(which='major', axis='both',ls='dashed')
##    
##    ax=axs[1,1]          
###    for j in range(0,len(coop[0,:,0])):
##    for jj in range(0,len(selcoop)):
##        j=selcoop[jj]
##        ax.plot(lambV,coop7[:,j,1],lin[j],label=lab[j]); ax.set_xlim(0., 1.); ax.set_ylim(0., 1.);  ax.set_xlabel('$\lambda$',size=12); #axs.set_ylabel('Level of cooperation');
##    ax.set_xticks([0,0.25,0.5,0.75,1]); ax.set_xticklabels(["0","0.25","0.5","0.75","1"]); ax.tick_params(axis='major', which='major', labelsize=8); ax.grid(which='major', axis='both',ls='dashed')
##    ax.text(1.1,0.48,"$M=7$", size=12 ,va='center')
##
##    margleft=0.2; margright=0.8; margtop=0.78; margbottom=0.12; wspace=hspace=0.15
##    f.subplots_adjust(hspace=hspace, wspace=wspace, right=margright,top=margtop,bottom=margbottom, left=margleft)
###        for i in range(0,len(groups)):
###            
###            mr=0.06; hh=(margtop-margbottom)/len(groups);  hib=hh-0.11; botb=margtop-hh*(i+1)+0.11-0.015*i; 
###            
###            #botb=(margtop-margbottom)/2.+(i-np.floor(len(groups)/2.))*0.2 ; hib=0.03
###            cbar_ax = f.add_axes([margright+0.13, botb, 0.2, hib])
###            hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
###            step=0.2; ti=np.arange(vmin,vmax+step,step); ti_s=["%.1f" % x for x in ti];  # ti_s[0]='<'+ti_s[0]
###            hb.set_ticks(ti)
###            hb.set_ticklabels(ti_s)
###            cbar_ax.tick_params(labelsize=7)
###            cbar_ax.set_title(nameg[i],size=8,color=mpl.cm.get_cmap(comaps[i])(1.))
##    
##    f.text(margleft-0.13, (margtop-margbottom)/2.+margbottom, 'Cooperation level', va='center', rotation='vertical',size=12)
##    #f.text((margright-margleft)/2+margleft, margbottom-0.1, '$\lambda$', ha='center',size=12)
##
##
##    h, l = ax00.get_legend_handles_labels()
##    ph = ax00.plot([],marker="", ls="")[0]
###    handles = [ph,h[0],h[1],ph,h[2],h[3],ph,h[4],h[5],ph,h[6],h[7],ph,h[8],h[9] ]
###    labels = ["$c_S=0$",lab[0],lab[1],"$c_S=0.1$",lab[2],lab[3],"$c_S=0.3$",lab[4],lab[5],"$c_S=0.5$",lab[6],lab[7]," ",lab[8],lab[9] ]
##    handles = [ph,h[0],h[1],ph,ph,ph,ph,h[2],h[3],ph,ph,ph,ph,h[4],h[5] ]
##    labels = ["$c_S=0$",lab[0],lab[1]," "," "," ","$c_S=0.3$",lab[2],lab[3]," "," "," "," ",lab[4],lab[5] ]
##    leg=ax00.legend(handles, labels, bbox_to_anchor=(0., margtop+0.47, 2+wspace, .102), loc=8,
##           ncol=5, mode="expand", borderaxespad=0.,fontsize=8,edgecolor='black')
##    for t in leg._legend_handle_box.get_children():
##        for hpack in t.get_children()[0:1]:
##            hpack.get_children()[0].set_width(0)
##            
##    f.savefig('coop_mechanisms_w1_r10_reduc.eps', dpi=300)
##    f.clf() 
###        
#########################################################################



    
######  Extracting graph of invasions ##########################
#
#####----- separating strategies ------
##    gST=[[28]]+[[20]]+[[35]]+[[33]]+[[i] for i in gSF]+[[i] for i in gSFm]+[[i] for i in gMBc] +[[i] for i in gMBd] +[gNO2]    
##    #gST=[[28]]+[[20]]+[[29]]+[[35]]+[[33]]+[[39]]+[[i] for i in gSF]+[[i] for i in gSFm]+[[i] for i in gMBc] +[[i] for i in gMBd] #+[gNOt]
##    colg= ['Red']*2+['Grey30']*2+['Mulberry']*len(gSF)+['Green']*len(gSFm)+['Cyan']*len(gMBc)+['Blue']*len(gMBd) +['Gray05']
##
##    labg=['']*len(gST)
##    STRmPUR=declareSTR(0)
##    for i in range(0,len(gST)):
##        labg[i]=str(STRmPUR[gST[i][0]].astype(int))
##    labg[-1]="others"
######-------------------------------
#
####----- groups of strategies ------
#    gN = [x for x in gALL if x not in gSC+gSD+gSF+gSFm+gMBc+gMBd+gSCO+gSDO ] 
#    gST=[gSC,gSD,gSCO,gSDO,gSF,gSFm,gMB,gN]
#    colg= ['Red']+['Red']+['Gray30']+['Gray30']+['Purple']+['Green']+['NavyBlue'] +['Gray05']
#    labg= ['SC']+['SD']+['SC-O']+['SD-O']+['FR-C']+['FR-O']+['FR-D'] +['others']
##    gN = [x for x in gALL if x not in gSC+gSD+gSF+gSFm+gMBc+gMBd] 
##    gST=[gSC,gSD,gSF,gSFm,gMB,gN]
##    colg= ['Red']+['Red']+['Purple']+['Green']+['NavyBlue'] +['Gray05']
##    labg= ['SC']+['SD']+['FR-C']+['FR-O']+['FR-D'] +['others']
#####-------------------------------    
#    
#####----- test ------
##    gST=[[28],[35]]
##    colg= ['Red']+['Gray30'] #+['Purple']+['Green']+['NavyBlue'] +['Gray05']
##    labg= ['SC']+['SD'] #+['SF-C']+['SF-01']+['SF-D + SF-10'] +['others']
######-------------------------------   
#
#    beta=1.
#    Z=100
#    N=9
#    M=5
#    Q=4.5
#    lamb=0.9 #0.5 #0.8
#    cs0=0.5
#    eps=0.01
#    w=1.
#    b0=10.
#    #H, L
#    c1=1. #2.5
#    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#    cs=np.array([1., 1.]) *cs0 *c1  #*0.06 *c  *0.8
#    b=np.array([1., 0.]) *b0 *c1   #7*c  
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
##    STRmPUR=declareSTR(0)
##    #print(STRmPUR[28],STRmPUR[35])
##    #print(fixM[28,35]*64*100,fixM[35,28]*64*100)
##    print(fixM[:,63]*64*100)
##    print(fixM[:,28]*64*100)
##    print((fixM[:,28]-fixM[:,63])*64*100)
##    print(SD[28],SD[63])
##    print(Mred[4,0]*64*100)
#
#
#    nSTg=[len(g) for g in gST]
#    name='NET_M_'+str(M)+'_Q_'+str(Q)+'_lamb_'+str(lamb)+'_cs_'+str(cs0)+'_b_'+str(b0)
#    plotNETgroup(name,Mred,SDred,labg,colg,nSTg,Z)
#
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
#    f.savefig('equiv_00-10.eps', dpi=300)
#    f.clf()       
#    
#    1.-Pr_lessM > 1.-1./Z ; Pr_lessM < 1./Z             #  ---> xx10xx equivalent to xx11xx. Start D and enter Nc>=M by mistakes
#    print(1.-Pr_lessM)  # never commit enough mitakes, unless w=1
#
#############################################

    
####################### find drift groups ###############################3 
    beta=1.
    Z=100
    N=9
    M=5
    Q=4.5
    lamb=0.5 #0.5 #0.8
    eps=0.01
    w=1.
    #H, L
    c1=1. #2.5
    c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
    cs=np.array([1., 1.]) *0.3 *c1  #*0.06 *c  *0.8
    b=np.array([1., 0.]) *10. *c1   #7*c  

    STRmPUR=declareSTR(0)
    
    doINI(N,Z,M,Q,eps,w)
    SD=doREST(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD,w)
    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
    SSDsort=SSD[np.argsort(SSD[..., 8])] 
    for i in range(0,len(SSDsort)):
      print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) 
     
    STRmPUR=declareSTR(0)
    expb=np.exp(-beta)
    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)
    fixMvec=readfixMvec(labelfile)
    fixM=calcFIXM(coef,expb,Z,fixMvec)

    groups=findDRIFTgroup(fixM,100)
    print(groups)

#    STRmPUR=declareSTR_SIG(0)
#    expb=np.exp(-beta)
#    coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
#    labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_SIG'
#    fixMvec=readfixMvec(labelfile)
#    fixM=calcFIXM(coef,expb,Z,fixMvec)
#    
#    doINI_SIG(N,Z,M,Q,eps,w)
#    SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#    PAYhomo,COOPhomo,COOPtot=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)
#    SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
#    SSDsort=SSD[np.argsort(SSD[..., 8])] 
#    for i in range(0,len(SSDsort)):
#       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])



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
######------ horizontal ---------------------    
##    lv=[0,0.2,0.4,0.5,0.6,0.8,1.]
##    lambV=   np.transpose(np.array([lv,lv]))
##    cs1V=  c1*np.transpose(np.array([[0.01]*len(lambV),[0.5]*len(lambV)])) 
##    bV=c1*np.transpose(np.array([[15]*len(lambV)]*2))
##    epsV= np.transpose(np.array([[0.01]*len(lambV)]*2))  
##    QV=    np.transpose( np.array([[4.5]*len(lambV)]*2))
##    MV=    np.transpose( np.array([[3]*len(lambV)]*2))
#####-------------------------------------- 
#
###------ horver ---------------------    
#    lv=[0,0.2,0.4,0.5,0.6,0.7,0.8,1.] # 2 horizontal
#    csv=np.array([0, 0.2,0.4,0.6,0.8,1.,1.2,1.5]) # 1 vertical
#    lambV=   np.transpose(np.array([lv,lv,[0.5]*len(lv)]))
#    cs1V=  c1*np.transpose(np.array([[0.5]*len(lv),[1.5]*len(lv),csv])) 
#    
#    
#    bV=c1*np.transpose(np.array([[20]*len(lambV)]*3))
#    epsV= np.transpose(np.array([[0.01]*len(lambV)]*3))  
#    QV=    np.transpose( np.array([[4.5]*len(lambV)]*3))
#    MV=    np.transpose( np.array([[7]*len(lambV)]*3))
##-------------------------------------- 
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
##    STv=np.array([20, 22, 28, 30, 33, 35, 43, 51, 53, 55, 57, 58, 59, 61, 63])
##    STv=np.array(range(0,64))
##    STv=np.array([20,28, # SC
##                  22,30, # SC mut
##                  33,35, # SD
##                  41,43, # SD mut
##                  48,49,50,51, # SF
##                  56,57,58,59, # SF mut
##                  52, 53, 54, 55, # MB C
##                  60, 61, 62, 63, # MB D
##                  0, 4, 8, 12,
##                  1, 3, 5, 7, 9, 11, 13, 15,
##                  2, 6, 10, 14,
##                  17, 19,
##                  21, 23, 29, 31,
##                  25, 27,
##                  36, 44,
##                  37, 39, 45, 47,
##                  38, 46,
##                  16, 18, 24, 26, 32, 34, 40, 45 ])  # not in groups
##    STvC=(['xkcd:dark red','xkcd:dark red', # SC
##          'xkcd:salmon','xkcd:salmon', # SC mut
##          'xkcd:dark grey','xkcd:dark grey', # SD
##          'xkcd:grey','xkcd:grey', # SD mut 21     1   0     1   0   1   0     9.24e-01   0.50     0.00
##          #'xkcd:sienna','xkcd:sienna', # SD
##          #'xkcd:tan','xkcd:tan', # SD mut 21     1   0     1   0   1   0     9.24e-01   0.50     0.00
##          #'xkcd:dark green','xkcd:dark green','xkcd:dark green','xkcd:dark green', # SF
##          'xkcd:violet','xkcd:violet','xkcd:violet','xkcd:violet', # SF
##          'xkcd:green','xkcd:green','xkcd:green','xkcd:green', # SF mut
##          'xkcd:medium blue','xkcd:medium blue','xkcd:medium blue','xkcd:medium blue', # MB C
##          'xkcd:dark blue', 'xkcd:dark blue', 'xkcd:dark blue', 'xkcd:dark blue'] # MB D
##          + ['xkcd:tan']*(64-24) ) # all the others
##    
##    STRp=declareSTR(0)
##    labup=['']*len(STv)
##    for i in range(0,len(STv)):
##        sen=str(STRp[STv[i],:])
##        labup[i]=sen.replace(". ","").replace("[","").replace(".]","").replace(" ","")
##    
##    label='BAR_w09_horit_Q4.5_M3_r15' #'BAR_ttt'
##    vmax=0.4
##    ncol=2
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
###    j=2    
###    for i in range(0,len(cs1V)):
###        nrow=len(cs1V)
###        npl=ncol*i+j+1     
###        axs = f.add_subplot(nrow,ncol,npl)
###        f.subplots_adjust(wspace=0.0,left=0.05,right=0.95,top=0.93,bottom=0.1)
###        if i+1!=len(cs1V):
###            labx=[]
###        else:
###            labx=labup
###        plot_BAR(labx,STv,STvC,axs,beta,Z,N,MV[i,j],QV[i,j],lambV[i,j],epsV[i,j],w,c1,cs1V[i,j],bV[i,j],vmax)
###        axs.yaxis.tick_right()
##    f.text(0.5, 0.96, 'Stationary Distribution', rotation=0, va='center', ha='center',size=8)
##    #f.text(0.05, 0.5, 'Stationary Distribution', rotation=90, va='center', ha='center',size=8)
##    f.savefig(label+'.pdf', dpi=300)
##    f.clf()
#    
#
##------- only chosen strategies -------------------------------------------------------
#    STv=np.array([28,20,29, # SC
#                  #16,17,18,19,21,22,23,24,25,26,27,30,31, # 10 non SC
#                  35, 33, 39, # SD
#                  #32,34,36,37,38,40,41,42,43,44,45,46,47, # 01 non SD
#                  48,49,50,51, # SF-C
#                  56,57,58,59, # SF-01
#                  52, 53, 54, 55, # SF-10
#                  60, 61, 62, 63, # SF-D
#                  ])  
#    STvC=(['xkcd:dark red']*3  # SC
#          #+['xkcd:salmon']*(16-3) # 10 non SC
#          +['xkcd:dark grey']*3  # SC
#          #+['xkcd:grey']*(16-3) # 10 non SC
#          +['xkcd:violet']*4  # SF-C
#          +['xkcd:green']*4  # SF-01
#          +['xkcd:medium blue']*4  # SF-10
#          +['xkcd:dark blue']*4)  # SF-D          
#
#            
#    STRp=declareSTR(0)
#    labup=['']*len(STv)
#    for i in range(0,len(STv)):
#        sen=str(STRp[STv[i],:])
#        labup[i]=sen.replace(". ","").replace("[","").replace(".]","").replace(" ","")
#    
#    label='BAR_w09_horver_Q4.5_M7' #'BAR_ttt'
#    vmax=0.5
#    ncol=3
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
#    j=2    
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
#        axs.yaxis.tick_right()
#    f.text(0.5, 0.96, 'Stationary Distribution', rotation=0, va='center', ha='center',size=8)
#    #f.text(0.05, 0.5, 'Stationary Distribution', rotation=90, va='center', ha='center',size=8)
#    f.savefig(label+'.eps', dpi=300)
#    f.clf()
##--------------------------------------------------
#    
#################################################################################################################    
   
    

############# SAVE BIG MATRIX for Q-M plots ################### 
#    beta=1.
#    Z=100
#    N=18
#    eps=0.01
#
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#
#    lambV=np.linspace(0,1,31)
#    
#    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
#    MV= np.array([2,6,10,14,18]) #np.array([1,3,5,7,9])  #np.array([2,6,10,14,18]) #np.array([1,3,5,7,9])   #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5])   #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
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
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#    
#    bVo=np.array([20.])
##    w=0.9
##    c1=1. #0.5
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N9_beta5_b10_w09_NEW_1.npy',bigmatSD)
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
#    c1=1.
#    csV=c1* csVo #np.linspace(0,0.3,51)
#    bV=c1* bVo
#    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
#    np.save('file_SD_N18_beta1_b20_w1_NEW_1.npy',bigmatSD)
#
#
###################################################################    
  


###################  PLOTS Q-M  ############################################################################    
##    bigmatCOOP=np.load('file_COOP_N9_beta5_b20_w1_X.npy')
##    lab='SD_N9_beta05_b20_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW.npy')
##    PAYhomo,COOPhomo,COOPtot=doHOMO(lamb,eps,N,M,Q,b,c,cs,SD)
# 
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #n  #np.array([5,6,7,8,9])   #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) # #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=0.5
#    csV= csVo
#    bV=c1* bVo
#
#    lab='SD_N9_beta05_b10_w09'
#    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
#    ext='eps'
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
##    plot_SDcslambDIF_agre('SD_gS11_'+lab,[gS11],['Greys'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
#
##    #STRmPURE=declareSTR(0)
##    STs00,STs11,STs10,STs01,STsign, STsignonly, STmem, STmemonly, STsignmem=classST()
##    #np.set_printoptions(threshold=np.inf)
##    #print(STRmPURE[STmemonly])
##    vmin=0; vmax=0.5
##    plot_SDcslambDIF_agre('REC_'+lab,[STmemonly],['Greys'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##    plot_SDcslambDIF_agre('SIG_'+lab,[STsignonly],['Greys'],0,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
#    
##    lab='SD_N9_beta1_b10_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
##     bigmatSD_abun_09=np.load('file_SD_N9_beta05_b20_w09_bper_abun.npy') 
##    bigmatSD_abun_1=np.load('file_SD_N9_beta05_b20_w1_bper_abun.npy')    ext='eps'
##    vmin=0.4; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF,gSFm,gSC,gSD]
##    #ngroups=['Blues','Purples','Greens','Reds','Greys']
##    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
##    ngroups=[r'FR-D + FR-10$^{M\geqslant5}$', r'FR-C',r'FR-01',r'SC + SC$_{C}$$^{M\geqslant5}$',r'SD + SD$_{C}$$^{M\geqslant5}$']
##
##    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##
##    lab='SD_N9_beta1_b10_w1'
##    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
##    ext='eps'
##    vmin=0.4; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF+gMBc,gSFm,gSCt,gSDt]
##    #ngroups=['Blues','Purples','Greens','Reds','Greys']
##    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
##    ngroups=[r'FR-D + FR-10$^{M>5}$', r'FR-C + FR-10$^{M<5}$',r'FR-01',r'SC + SC$_{C}$$^{M>5}$ + SC$_{D}$$^{M<5}$',r'SD + SD$_{C}$$^{M>5}$ + SD$_{D}$$^{M<5}$']
##
##    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
#
##    lab='SD_N9_beta1_b10_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
##    ext='eps'
##    vmin=0.5; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys','Oranges']
##    groups=[gMB,gSF,gSFm,gSC,[22,30],gSD]
##    #ngroups=['Blues','Purples','Greens','Reds','Greys']
##    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
##    ngroups=[r'FR-D + FR-10$^{M\geqslant5}$', r'FR-C',r'FR-01',r'SC + SC$_{C}$$^{M\geqslant5}$',r'[10 $\ast$001]',r'SD + SD$_{C}$$^{M\geqslant5}$']
##
##    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
#
#    lab='SD_N9_beta1_b10_w1'
#    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
#    ext='eps'
#    vmin=0.5; vmax=1.
#    comaps=['Blues','Purples','Greens','Reds','Greys']
#    groups=[gMB,gSF,gSFm,gSC+gSD,[22,30, 43,41]]
#    #ngroups=['Blues','Purples','Greens','Reds','Greys']
#    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
#    #ngroups=[r'FR-D + FR-10$^{M>5}$', r'FR-C + FR-10$^{M<5}$',r'FR-01',r'SC + SC$_{C}$$^{M>5}$ + SC$_{D}$$^{M<5}$',r'[10 $\ast$001]',r'SD + SD$_{C}$$^{M>5}$ + SD$_{D}$$^{M<5}$']
#    ngroups=[r'FR-D',r'FR-C',r'FR-O',r'SC + SD',r'SC-O + SD-O']
#
#    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##
##
##    lab='SD_N9_beta1_b10_w09'
##    bigmatSD=np.load('file_'+lab+'_NEW_1.npy')    
##    ext='eps'
##    vmin=0.5; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF,gSFm,gSC+gSD,[22,30, 43,41]]
##    #ngroups=['Blues','Purples','Greens','Reds','Greys']
##    #ngroups=[r'[00 00$\ast$$\ast$] + [00 10$\ast$$\ast$]', r'[00 11$\ast$$\ast$] + [00 10$\ast$$\ast$]',r'[00 01$\ast$$\ast$]',r'[10 0011] + [10 1011] + [10 0010]',r'[01 1100] + [01 1110] + [01 1000]']
##    #ngroups=[r'FR-D + FR-10$^{M>5}$', r'FR-C + FR-10$^{M<5}$',r'FR-01',r'SC + SC$_{C}$$^{M>5}$ + SC$_{D}$$^{M<5}$',r'[10 $\ast$001]',r'SD + SD$_{C}$$^{M>5}$ + SD$_{D}$$^{M<5}$']
##    ngroups=[r'FR-D',r'FR-C',r'FR-O',r'SC + SD',r'SC-O + SD-O']
##
##    plot_SDcslambDIF_agre('SD_agre_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
##
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
##################################################################################################################



###################  PLOTS signal action panels  ############################################################################    
#    import matplotlib.pyplot as plt
#    import matplotlib as mpl
# 
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#    M=5
#    
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=1.
#    csV= csVo
#    bV=c1* bVo
#
#    bigmatSD=np.load('file_SD_N9_beta1_b10_w1_NEW_1.npy')
#    STs00,STs11,STs10,STs01,STsign, STsignonly, STmem, STmemonly, STsignmem=classST()
#           
#    alp=1.; step=0.02; iM=2; iQ=2
#    f,axs=plt.subplots(nrows=2, ncols=3, sharex='none', sharey='all' )
#    f.subplots_adjust(hspace=0.7, wspace=0.2)  
#    
#    vmin=-1e-10; vmax=1.
#    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#    axs[0,0].contourf(lambV,csVo,np.sum(bigmatSD[gS00,:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Reds')
#    axs[0,0].set_title(r'[00 $\ast$$\ast$$\ast$$\ast$]', size=8 )
#    axs[0,1].contourf(lambV,csVo,np.sum(bigmatSD[gS01,:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Reds')
#    axs[0,1].set_title(r'[01 $\ast$$\ast$$\ast$$\ast$]', size=8 )
#    axs[0,2].contourf(lambV,csVo,np.sum(bigmatSD[gS10,:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Reds')
#    axs[0,2].set_title(r'[10 $\ast$$\ast$$\ast$$\ast$]', size=8 )
#    f.subplots_adjust(left=0.1,right=0.8)
#    cbar_ax = f.add_axes([0.85, 0.6, 0.03, 0.28])
#    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Signal',cmap='Reds')
#    cbar_ax.tick_params(labelsize=8)
#    
#    print(STsignonly)
#    vmin=-1e-10; vmax=0.5
#    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#    axs[1,0].contourf(lambV,csVo,np.sum(bigmatSD[[0,63],:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Blues')
#    axs[1,0].set_title(r'[$\ast$$\ast$ 0000]+[$\ast$$\ast$ 1111]', size=8)
#    axs[1,1].contourf(lambV,csVo,np.sum(bigmatSD[STmemonly,:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Blues')
#    axs[1,1].set_title(r'[$\ast$$\ast$ 1010]+[$\ast$$\ast$ 0101]', size=8)
#    axs[1,2].contourf(lambV,csVo,np.sum(bigmatSD[STsignonly,:,:,0,iM,iQ],axis=0),np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap='Blues')
#    axs[1,2].set_title(r'[$\ast$$\ast$ 1100]+[$\ast$$\ast$ 0011]', size=8)   
#    cbar_ax = f.add_axes([0.85, 0.123, 0.03, 0.28])
#    hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Action',cmap='Blues')
#    step=0.1; ti=np.arange(0.,vmax+step,step); ti_s=["%.1f" % x for x in ti]; ti_s[-1]='>'+ti_s[-1]
#    hb.set_ticks(ti)
#    hb.set_ticklabels(ti_s)
#    cbar_ax.tick_params(labelsize=8)
#   
#    for i in range(0,len(axs)):
#        axs[i,0].set_ylabel(r'$c_s$', size=12)
#        for j in range(0,len(axs[0])):        
#                axs[i,j].set_xlabel(r'$\lambda$',size=12)
#                axs[i,j].set_xticks([0,0.25,0.5,0.75,1]); #axs[iM,iQ].set_yticks([0,0.5,1])
#                axs[i,j].set_xticklabels(["0","0.25","0.5","0.75","1"]); #axs[iM,iQ].set_yticklabels(["0","0.5","1"])
#                axs[i,j].set_yticks([0,0.25,0.5,0.75,1]);
#                axs[i,j].set_yticklabels(["0","0.25","0.5","0.75","1"]);
#                axs[i,j].tick_params(axis='both', which='major', labelsize=8)
#                axs[i,j].grid(which='both', axis='both',ls='dashed')
#   
#    f.savefig('signal-action_mechanism_w1_M5_1.eps', dpi=300)
#    f.clf()    
#   
####################################################################################################################



###################  PLOTS r-N  ############################################################################    
#
# 
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#    M=5
#    
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=1.
#    csV= csVo
#    bV=c1* bVo
#
#    bigmatSD_5_9=np.load('file_SD_N9_beta1_b5_w1_NEW_1.npy') 
#    bigmatSD_10_9=np.load('file_SD_N9_beta1_b10_w1_NEW_1.npy') 
#    bigmatSD_20_9=np.load('file_SD_N9_beta1_b20_w1_NEW_1.npy') 
#    bigmatSD_30_9=np.load('file_SD_N9_beta1_b30_w1_NEW_1.npy') 
#    bigmatSD_5_18=np.load('file_SD_N18_beta1_b5_w1_NEW_1.npy') 
#    bigmatSD_10_18=np.load('file_SD_N18_beta1_b10_w1_NEW_1.npy') 
#    bigmatSD_20_18=np.load('file_SD_N18_beta1_b20_w1_NEW_1.npy')    
#    bigmatSD_30_18=np.load('file_SD_N18_beta1_b30_w1_NEW_1.npy') 
#    labright=['$r=$5','$r=$10','$r=$20','$r=$30']
#    labup=['$N$=9','$N$=18']
#    bigmatSDlist=[[bigmatSD_5_9,bigmatSD_10_9,bigmatSD_20_9,bigmatSD_30_9],[bigmatSD_5_18,bigmatSD_10_18,bigmatSD_20_18,bigmatSD_30_18]]
#    bigmatSDlist=list(map(list, zip(*bigmatSDlist))) # transposing list
#    #iQ=np.zeros((3,2),int); iM=np.zeros((3,2),int)
#    iQ=2; iM=2
#    
#    ext='eps'
#    vmin=0.5; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF,gSFm,gSC,gSD]
##    ngroups=[r'FR-D + FR-10', r'FR-C',r'FR-01',r'SC + SC$_{C}$',r'SD + SD$_{C}$']
#    comaps=['Blues','Purples','Greens','Reds','Greys']
#    groups=[gMB,gSF,gSFm,gSC+gSD,[22,30, 43,41]]
#    ngroups=[r'FR-D',r'FR-C',r'FR-O',r'SC + SD',r'SC-O + SD-O']
#
#    plot_SDspace_agre('SD_agre_r-N_w1',groups,comaps,ngroups,bigmatSDlist,csV,lambV,iM,iQ,M,labup,labright,vmin,vmax,ext)
##    
##    bigmatSD_15_9=np.load('file_SD_N9_beta05_b15_w09_NEW.npy') 
##    bigmatSD_20_9=np.load('file_SD_N9_beta05_b20_w09_NEW.npy') 
##    bigmatSD_30_9=np.load('file_SD_N9_beta05_b30_w09_NEW.npy') 
##    bigmatSD_15_18=np.load('file_SD_N18_beta05_b15_w09_NEW.npy') 
##    bigmatSD_20_18=np.load('file_SD_N18_beta05_b20_w09_NEW.npy') 
##    bigmatSD_30_18=np.load('file_SD_N18_beta05_b30_w09_NEW.npy')     
##    labright=['$r=$15','$r=$20','$r=$30']
##    labup=['$N$=9','$N$=18']
##    bigmatSDlist=[[bigmatSD_15_9,bigmatSD_20_9,bigmatSD_30_9],[bigmatSD_15_18,bigmatSD_20_18,bigmatSD_30_18]]
##    bigmatSDlist=list(map(list, zip(*bigmatSDlist))) # transposing list
##    #iQ=np.zeros((3,2),int); iM=np.zeros((3,2),int)
##    iQ=2; iM=2
##    
##    ext='eps'
##    vmin=0.4; vmax=1.
##    comaps=['Blues','Purples','Greens','Reds','Greys']
##    groups=[gMB,gSF,gSFm,gSC,gSD]
##    ngroups=[r'FR-D + FR-10', r'FR-C',r'FR-01',r'SC + SC$_{C}$',r'SD + SD$_{C}$']
##
##    plot_SDspace_agre('SD_agre_r-N_w09',groups,comaps,ngroups,bigmatSDlist,csV,lambV,iM,iQ,M,labup,labright,vmin,vmax,ext)
#   
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





########## SAVE BIG MATRIX for Q-M plots - ONLY SIG ################### 
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01
#
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#
#    lambV=np.linspace(0,1,31)
#    
#    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
#    MV= np.array([1,3,5,7,9])  #np.array([2,6,10,14,18]) #np.array([1,3,5,7,9])   #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5])   #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
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
#    STRmPURE=declareSTR_SIG(0); nSTR=STRmPURE.shape[0];
#    
#    bVo=np.array([10.])
##    w=0.9
##    c1=1. #0.5
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N9_beta5_b10_w09_NEW_1_SIG.npy',bigmatSD)
####    np.save('file_COOP_N9_beta5_b10_w1_X.npy',bigmatCOOP)
###    w=0.9
###    c1=0.5
###    csV=c1* csVo #np.linspace(0,0.3,51)
###    bV=c1* bVo
###    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
###    np.save('file_SD_N9_beta1_b20_w09_NEWtest.npy',bigmatSD)
##
###    bVo=np.array([30.])
###    w=0.9
###    c1=0.5
###    csV=c1* csVo #np.linspace(0,0.3,51)
###    bV=c1* bVo
###    bigmatCOOP,bigmatSD=doMATSD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
###    np.save('file_SD_N18_beta05_b30_w09_NEW.npy',bigmatSD)
###    
#    w=1.
#    c1=1.
#    csV=c1* csVo #np.linspace(0,0.3,51)
#    bV=c1* bVo
#    bigmatCOOP,bigmatSD=doMATSD_SIG(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
#    np.save('file_SD_N9_beta1_b10_w1_NEW_1_SIG.npy',bigmatSD)
#
###################

###################  PLOTS Q-M - ONLY SIG ############################################################################    
#
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #n  #np.array([5,6,7,8,9])   #np.array([3,9,15,21,27])     
#    QV= np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) # #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=0.5
#    csV= csVo
#    bV=c1* bVo
#
#
#    lab='SD_N9_beta1_b10_w1'
#    bigmatSD=np.load('file_'+lab+'_NEW_1_SIG.npy')    
#    ext='eps'
#    vmin=0.5; vmax=1.
#    comaps=['Blues','Purples','Reds']
#    groups=[[14, 15],[12, 13],[6, 9]]
#    ngroups=[r'FR-D',r'FR-C',r'SC + SD']
#
#    plot_SDcslambDIF_agre('SD_agre_SIG_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
# 
#####################################

# ######################## find drift groups - ONLY SIG ###############################3 
#     beta=1.
#     Z=100
#     N=9
#     M=5
#     Q=4.5
#     lamb=0.8 #0.5 #0.8
#     eps=0.01
#     w=1.
#     #H, L
#     c1=1. #2.5
#     c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#     cs=np.array([1., 1.]) *0.3 *c1  #*0.06 *c  *0.8
#     b=np.array([1., 0.]) *10. *c1   #7*c  

#     STRmPUR=declareSTR_SIG(0)
#     expb=np.exp(-beta)
#     coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
#     labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_SIG'
#     fixMvec=readfixMvec(labelfile)
#     fixM=calcFIXM(coef,expb,Z,fixMvec)
    
#     doINI_SIG(N,Z,M,Q,eps,w)
#     SD=doREST_SIG(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#     PAYhomo,COOPhomo,COOPtot=doHOMO_SIG(lamb,eps,N,M,Q,b,c,cs,SD,w)
#     SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
#     SSDsort=SSD[np.argsort(SSD[..., 8])] 
#     for i in range(0,len(SSDsort)):
#        print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
# ############


###### SAVE BIG MATRIX for Q-M plots - ONLY REC ################### 
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01
#
#    STRmPURE=declareSTR(0); nSTR=STRmPURE.shape[0];
#
#    lambV=np.linspace(0,1,31)
#    
#    csVo= np.linspace(0,1,31) #np.linspace(0,0.3,51)
#    MV= np.array([1,3,5,7,9])  #np.array([2,6,10,14,18]) #np.array([1,3,5,7,9])   #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
#    QV= np.array([1.]) #np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5])   #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
#    
#    
#    STRmPURE=declareSTR_REC(0); nSTR=STRmPURE.shape[0];
#    
#    bVo=np.array([10.])
#
#    w=1.
#    c1=1.
#    csV=c1* csVo #np.linspace(0,0.3,51)
#    bV=c1* bVo
#    bigmatCOOP,bigmatSD=doMATSD_REC(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
#    np.save('file_SD_N9_beta1_b10_w1_NEW_1_REC.npy',bigmatSD)
#
###################

###################  PLOTS Q-M - ONLY REC ############################################################################    
#
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#    MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #n  #np.array([5,6,7,8,9])   #np.array([3,9,15,21,27])     
#    QV= np.array([1.]) #np.array([2., 4., 9.5,13.5,17.5]) # #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
#
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=0.5
#    csV= csVo
#    bV=c1* bVo
#
#
#    lab='SD_N9_beta1_b10_w1'
#    bigmatSD=np.load('file_'+lab+'_NEW_1_REC.npy')    
#    ext='eps'
#    vmin=0.5; vmax=1.
#    comaps=['Blues','Purples','Greens','Reds']
#    groups=[[0],[3],[2], [1]]
#    ngroups=[r'D',r'C',r'aTFT (O)', r'TFT']
#
#    plot_SDcslambDIF_agre('SD_agre_REC_'+lab,groups,comaps,ngroups,bigmatSD,csV,lambV,bV,MV,QV,vmin,vmax,ext)
# 
###########

# ######################## find drift groups - ONLY REC ###############################3 
#     beta=1.
#     Z=100
#     N=9
#     M=5
#     Q=4.5
#     lamb=0.8 #0.5 #0.8
#     eps=0.01
#     w=1.
#     #H, L
#     c1=1. #2.5
#     c=np.array([1., 1.])  *1. *c1    #*0.3  *0.8
#     cs=np.array([1., 1.]) *0.1 *c1  #*0.06 *c  *0.8
#     b=np.array([1., 0.]) *10. *c1   #7*c  

#     STRmPUR=declareSTR_REC(0)
#     expb=np.exp(-beta)
#     coef=np.array([[b[0]*lamb, -c[0]*lamb, -cs[0]*lamb],[b[1]*(1.-lamb), -c[1]*(1.-lamb), -cs[1]*(1.-lamb)]])
#     labelfile='GRIM_N_'+str(N)+'_M_'+str(M)+'_Q_'+str(Q)+'_eps_'+str(eps)+'_w_'+str(w)+'_REC'
#     fixMvec=readfixMvec(labelfile)
#     fixM=calcFIXM(coef,expb,Z,fixMvec)
    
#     doINI_REC(N,Z,M,Q,eps,w)
#     SD=doREST_REC(b,c,cs,lamb,beta,N,Z,M,Q,eps,w)
#     PAYhomo,COOPhomo,COOPtot=doHOMO_REC(lamb,eps,N,M,Q,b,c,cs,SD,w)
#     SSD=np.concatenate((STRmPUR,np.transpose([PAYhomo]),np.transpose([COOPhomo[:,0]]),SD),axis=1)
#     SSDsort=SSD[np.argsort(SSD[..., 8])] 
#     for i in range(0,len(SSDsort)):
#       print('{0:3.0f} {1:5.0f} {2:3.0f} {3:5.0f} {4:3.0f} {5:3.0f} {6:3.0f} {7:12.2e} {8:6.2f} {9:8.2f}'.format(np.argsort(SSD[..., 8])[i],SSDsort[i,0],SSDsort[i,1],SSDsort[i,2],SSDsort[i,3],SSDsort[i,4],SSDsort[i,5],SSDsort[i,6],SSDsort[i,7],SSDsort[i,8])) #print(SSDsort[i,:])
# ############
#
######## SAVE BIG MATRIX for Q-M plots - ONLY C+D ################### 
##    beta=1.
##    Z=100
##    N=9
##    eps=0.01
##
##
##    lambV=np.linspace(0,1,31)
##    
##    csVo= np.linspace(0,1,1)
##    MV= np.array([1,3,5,7,9])  #np.array([2,6,10,14,18]) #np.array([1,3,5,7,9])   #np.array([5,6,7,8,9])  #np.array([2,6,10,14,18]) #np.array([3,9,15,21,27])     
##    QV= np.array([1.]) #np.array([2., 4., 9.5,13.5,17.5]) #np.array([1., 2.5, 4.5,6.5,8.5])   #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    
##
##    
##    
##    STRmPURE=declareSTR_CD(0); nSTR=STRmPURE.shape[0];
##    
##    bVo=np.array([10.])
##
##    w=1.
##    c1=1.
##    csV=c1* csVo #np.linspace(0,0.3,51)
##    bV=c1* bVo
##    bigmatCOOP,bigmatSD=doMATSD_CD(beta,Z,N,nSTR,c1,csV,lambV,bV,MV,QV,w,eps)
##    np.save('file_SD_N9_beta1_b10_w1_NEW_1_C+D.npy',bigmatSD)
##
##########
#
#
# ############### Plot SIG only and REC, C+D
#     lambV=np.linspace(0,1,31)   
#     csVo= np.linspace(0,1,31) #np.linspace(0,2,31) #np.linspace(0,0.3,51)
#     bVo=np.array([10.])
#     MV= np.array([1,3,5,7,9]) #np.array([2,6,10,14,18]) #n  #np.array([5,6,7,8,9])   #np.array([3,9,15,21,27])     
#     QV= np.array([1.,2.5, 4.5,6.5,8.5]) #np.array([2., 4., 9.5,13.5,17.5]) # #np.array([1., 2.5, 4.5,6.5,8.5]) #np.array([1., 2.5, 4.5,6.5,8.5])  #np.array([4.5,5.5,6.5,7.5,8.5]) #np.array([2.0,6.0,9.0,10.0,14.0,18.0]) #np.array([3.0,9.0,13.5,15.0,21.0,27.0])    

#     beta=1.
#     Z=100
#     N=9
#     eps=0.01   

#     c1=0.5
#     csV= csVo
#     bV=c1* bVo


#     label='SD_N9_beta1_b10_w1'
#     bigmatSD=np.load('file_'+label+'_NEW_1_SIG.npy')    
#     ext='svg'
#     vmin=0.5; vmax=1.
#     comapsV=['Blues','Purples','Reds']
#     groups=[[14, 15],[12, 13],[6, 9]]
#     nameg=[r'FR-D',r'FR-C',r'SC + SD']
    
#     bigmatSDREC=np.load('file_'+label+'_NEW_1_REC.npy') 
#     bigmatSDCD=np.load('file_'+label+'_NEW_1_C+D.npy') 
    
#     import matplotlib.pyplot as plt
#     import matplotlib as mpl
#     alp=1.
#     lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
#     for i in range(0,len(groups)):
#         bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)
#     nr=bigmatAGR.shape[4]; nc=bigmatAGR.shape[5] 
#     f,axs=plt.subplots(nrows=nr, ncols=nc, sharex='all', sharey='all' )
#     f.subplots_adjust(hspace=0.2, wspace=0.2)

#     norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#     comaps=comapsV
#     for i in range(len(groups)):
#         comaps[i]=plt.get_cmap(comapsV[i])
#         comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
#     for iM in range(nr-1,-1,-1):
#         axs[iM,nc-1].text(1.1,0.48,"$M=%s$" % str(MV[iM]), size=9 ,va='center')
#         for iQ in range(nc-1,0,-1):
#             step=0.02
#             if MV[iM]>5:   # to avoid problems with [0010**], which is two places for w=1
#                 rg=range(len(groups)-1,-1,-1)     
#             else:
#                 rg=range(0,len(groups))
#             for i in rg:
#                 h=axs[iM,iQ].contourf(lambV,csV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
#             axs[iM,iQ].set_xticks([0,0.5,1]); axs[iM,iQ].set_yticks([0,0.5,1])
#             axs[iM,iQ].set_xticklabels(["0","0.5","1"]); axs[iM,iQ].set_yticklabels(["0","0.5","1"])
#             #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
#             #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
#             #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
#             #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
#             axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
#             axs[iM,iQ].grid(which='both', axis='both',ls='dashed')
#             axs[iM,iQ].set_xlim([0, 1]); axs[iM,iQ].set_ylim([0, 1])
#             if iM==0:
#                 axs[iM,iQ].set_title("$Q=%s$" % str(QV[iQ]), size=9 )
                
#     margbottom=0.15; margtop=0.87
#   #  f.text(0.0, 0.5, '$c_s$', va='center', rotation='vertical',size=12)
#     f.text(0.0, 0.5, '$D$', va='center', rotation='vertical',size=12)
#     if nameg==0:
#         margleft=0.1; margright=0.75;
#         f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
#         cbar_ax = f.add_axes([margright+0.1, margbottom, 1.-margleft-margright-0.12, margtop-margbottom])
#         hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,label='Probability',cmap=comaps[-1])
#     else:
#         margleft=0.09; margright=0.66;
#         f.subplots_adjust(right=margright,top=margtop,bottom=margbottom, left=margleft)
#         for i in range(0,len(groups)):
            
#             mr=0.06; hh=(margtop-0.45)/len(groups);  hib=hh-0.11; botb=margtop-hh*(i+1)+0.109-0.027*i; 
            
#             #botb=(margtop-margbottom)/2.+(i-np.floor(len(groups)/2.))*0.2 ; hib=0.03
#             cbar_ax = f.add_axes([margright+0.11, botb, 1.-margleft-margright-0.06, hib])
#             hb=mpl.colorbar.ColorbarBase(cbar_ax, norm=norm,cmap=comaps[i],orientation='horizontal')
#             step=0.25; ti=np.arange(vmin,vmax+step,step); ti_s=["%.2f" % x for x in ti];  # ti_s[0]='<'+ti_s[0]
#             hb.set_ticks(ti)
#             hb.set_ticklabels(ti_s)
#             cbar_ax.tick_params(labelsize=8)
#             cbar_ax.set_title(nameg[i],size=8,color=mpl.cm.get_cmap(comaps[i])(1.))
    
#     f.text((margright-margleft)/2+margleft, 0.04, '$\lambda$', ha='center',size=12)

#     iQ=0
#     for iM in range(nr-1,-1,-1):
#         axs[iM,iQ].plot(lambV,bigmatSDCD[1,0,:,0,iM,iQ],linewidth=0.8,color='Black',label='C (B)')
#         axs[iM,iQ].plot(lambV,bigmatSDREC[0,0,:,0,iM,iQ],color='Blue',label='D')
#         axs[iM,iQ].plot(lambV,bigmatSDREC[3,0,:,0,iM,iQ],color='Purple',label='C')
#         axs[iM,iQ].plot(lambV,bigmatSDREC[2,0,:,0,iM,iQ],'--',color='Green',label='O')
#         axs[iM,iQ].plot(lambV,bigmatSDREC[1,0,:,0,iM,iQ],'--',color='Orange',label='F')
#         axs[iM,iQ].set_xticks([0,0.5,1]); axs[iM,iQ].set_yticks([0,0.5,1])
#         axs[iM,iQ].set_xticklabels(["0","0.5","1"]); axs[iM,iQ].set_yticklabels(["0","0.5","1"])
#             #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
#             #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
#             #axs[iM,iQ].set_yticks([0,0.5,1.,1.5]); 
#             #axs[iM,iQ].set_yticklabels(["0","0.1","0.2","0.3"]); 
#         axs[iM,iQ].tick_params(axis='both', which='major', labelsize=8)
#         axs[iM,iQ].grid(which='both', axis='both',ls='dashed')
#         axs[iM,iQ].set_xlim([0, 1]); axs[iM,iQ].set_ylim([0, 1])
#     axs[0,0].set_title("B & R   S", size=10, fontweight='bold' )
#     axs[0,0].legend(loc=1,bbox_to_anchor=(8.5,-3), shadow=False, fontsize=8)
        

#     #hb.set_ticks(np.linspace(vmin,vmax,11))
# #    plt.show()
#     #f.text(0.874, 0.95, labup, va='center', ha='center',color='darkred',size=10)
#     #f.text(0.874, 0.08, labdown, va='center', ha='center',color='darkblue',size=10)

#     #for i in range(0,len(ext)):
#     f.savefig(label+'_SIG_REC_forSIG.'+'svg', dpi=600)
#     f.clf()
            
# ############


############ Plot SIG and SIG+REC
#    lambV=np.linspace(0,1,31)   
#    csVo= np.linspace(0,1,31) #np.linspace(0,2,31) #np.linspace(0,0.3,51)
#    bVo=np.array([10.])
#
#    beta=1.
#    Z=100
#    N=9
#    eps=0.01   
#
#    c1=0.5
#    csV= csVo
#    bV=c1* bVo
#
#    ext='eps'
#    vmin=0.5; vmax=1.
#    
#    import matplotlib.pyplot as plt
#    import matplotlib as mpl
#    alp=1.
#    step=0.02
#    iM=2; iQ=2
#
#    f,axs=plt.subplots(nrows=2, ncols=2, sharex='all')#, sharey='all' )
#    f.subplots_adjust(hspace=0.8, wspace=0.45)
#
#     
#
#    label='SD_N9_beta1_b10_w1'
#    bigmatSD=np.load('file_'+label+'_NEW_1_SIG.npy')    
#    comapsV=['Blues','Purples','Reds']
#    groups=[[14, 15],[12, 13],[6, 9]]
#    nameg=[r'FR-D',r'FR-C',r'SC + SD']
#        
#    lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
#    for i in range(0,len(groups)):
#        bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)
#    nr=bigmatAGR.shape[4]; nc=bigmatAGR.shape[5] 
#
#    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#    comaps=comapsV
#    for i in range(len(groups)):
#        comaps[i]=plt.get_cmap(comapsV[i])
#        comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
#   
#    rg=range(0,len(groups))              
#    for i in rg:
#        h=axs[1,0].contourf(lambV,csV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
#    axs[1,0].set_xticks([0,0.5,1]); axs[1,0].set_yticks([0,0.5,1])
#    axs[1,0].set_xticklabels(["0","0.5","1"]); axs[1,0].set_yticklabels(["0","0.5","1"])
#    axs[1,0].tick_params(axis='both', which='major', labelsize=8)
#    axs[1,0].grid(which='both', axis='both',ls='dashed')
#    axs[1,0].set_xlim([0, 1]); axs[1,0].set_ylim([0, 1])
#    axs[1,0].set_xlabel('$\lambda$',size=18)
#    axs[1,0].set_ylabel('$c_s$', rotation='vertical',size=18)
#    axs[1,0].set_title('S',size=18,pad=10)
#    
#    label='SD_N9_beta1_b10_w1'
#    bigmatSD=np.load('file_'+label+'_NEW_1.npy')    
#    comapsV=['Blues','Purples','Greens','Reds','Greys']
#    groups=[gMB,gSF,gSFm,gSC+gSD,[22,30, 43,41]]
#    nameg=[r'FR-D',r'FR-C',r'SC + SD']
#        
#    lAGR=list(bigmatSD.shape); del lAGR[0]; lAGR.insert(0,len(groups)); bigmatAGR=np.empty(lAGR)
#    for i in range(0,len(groups)):
#        bigmatAGR[i,:]=np.sum(bigmatSD[groups[i],...],axis=0)
#    nr=bigmatAGR.shape[4]; nc=bigmatAGR.shape[5] 
#
#    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#    comaps=comapsV
#    for i in range(len(groups)):
#        comaps[i]=plt.get_cmap(comapsV[i])
#        comaps[i]= truncate_colormap(comaps[i], 0.25, 1)
#  
#    rg=range(0,len(groups))              
#    for i in rg:
#        h=axs[1,1].contourf(lambV,csV,bigmatAGR[i,:,:,0,iM,iQ],np.arange(vmin,vmax+0.1,step),alpha=alp,vmin=vmin,vmax=vmax, cmap=comaps[i])
#    axs[1,1].set_xticks([0,0.5,1]); axs[1,1].set_yticks([0,0.5,1])
#    axs[1,1].set_xticklabels(["0","0.5","1"]); axs[1,1].set_yticklabels(["0","0.5","1"])
#    axs[1,1].tick_params(axis='both', which='major', labelsize=8)
#    axs[1,1].grid(which='both', axis='both',ls='dashed')
#    axs[1,1].set_xlim([0, 1]); axs[1,1].set_ylim([0, 1])
#    axs[1,1].set_xlabel('$\lambda$',size=18)    
#    axs[1,1].set_ylabel('$c_s$', rotation='vertical',size=18) 
#    axs[1,1].set_title('S + R',size=18,pad=10)
#
#    axs[1,0].tick_params(labelsize=14)
#    axs[1,1].tick_params(labelsize=14)    
#    f.subplots_adjust(top=1.5, bottom=0.15, left=0.13, right=0.95)
#    
#    f.savefig(label+'_SIG_vs_REC_1panel.'+ext, bbox_inches=mpl.transforms.Bbox([[0,0], [6, 3]]), dpi=600)
#    f.clf()
#            
########



