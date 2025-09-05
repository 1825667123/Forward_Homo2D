import numpy as np
import scipy.special as sp
import global_para_var as gv

def Compute_SctField(Eytot, Eysct):
    
    print("compute scattering fields...")
    Jy=np.empty(np.shape(Eytot),dtype=complex)
    for ifre in range(gv.nfreq):
        for itrans in range(gv.ntrtot):
            Jy[ifre,itrans,::]=1j*gv.xomega[ifre]*gv.cer[ifre]*gv.epsilon0*gv.xkai[ifre,::]*Eytot[ifre,itrans,::]
    
    rho=np.empty([gv.nrectot,gv.mz,gv.mx],dtype=float)
    for irece in range(gv.nrectot):
        rho[irece,:,:]=np.sqrt((gv.xrr[irece]-gv.xxf)**2 + (gv.zrr[irece]-gv.xzf)**2)
    for ifre in range(gv.nfreq):
        for itrans in range(gv.ntrtot):
            for irece in range(gv.nrectot):
                Eysct[ifre, itrans, irece]=-0.25*gv.xomega[ifre]*gv.mu0*gv.mur*np.sum((sp.jv(0, gv.xk[ifre]*rho[irece,:,:])-1j*sp.yv(0, gv.xk[ifre]*rho[irece,:,:]))*Jy[ifre,itrans,:,:])*gv.dxz

    return


