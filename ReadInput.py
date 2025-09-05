import numpy as np
import global_para_var as gv

def Read_Input():

    #(1) 读取背景介质参数
    print("read background parameters...")
    with open("input/background.inp", 'r') as fid:
        next(fid)  
        gv.ttcer, gv.xsig, gv.ttmiur = list(map(float,fid.readline().strip().split()))
    
    #(2) 读取发射机配置
    print("read transmitter config...")
    with open("input/tloc.inp", 'r') as fid:
        next(fid)  
        gv.ntrtot = int(fid.readline().strip())
        gv.xptr=np.empty(gv.ntrtot,dtype=float)
        gv.zptr=np.empty(gv.ntrtot,dtype=float)  
        next(fid)
        for m in range(gv.ntrtot):
            gv.xptr[m], gv.zptr[m] = list(map(float,fid.readline().strip().split()))

    #(3) 读取接收机配置
    print("read receiver config...")
    with open("input/rloc.inp", 'r') as fid:
        next(fid)  
        gv.nrectot = int(fid.readline().strip())
        gv.xrr=np.empty(gv.nrectot,dtype=float)
        gv.zrr=np.empty(gv.nrectot,dtype=float)  
        next(fid)
        for m in range(gv.nrectot):
            gv.xrr[m], gv.zrr[m] = list(map(float,fid.readline().strip().split()))

   #(4) 读取网格配置
    print("read cell config...")
    with open("input/cell.inp", 'r') as fid:
        next(fid)  
        gv.mx, gv.mz = list(map(int,fid.readline().strip().split()))
        next(fid) 
        gv.xx1, gv.zz1 = list(map(float,fid.readline().strip().split()))
        next(fid)
        gv.dx, gv.dz = list(map(float,fid.readline().strip().split()))
    gv.dxz=gv.dx*gv.dz
    gv.xxf=np.tile(np.arange(gv.xx1,gv.xx1+gv.dx*gv.mx,gv.dx),[gv.mz,1])
    gv.xzf=np.tile(np.arange(gv.zz1,gv.zz1+gv.dz*gv.mz,gv.dz).reshape(-1,1),[1,gv.mx])

   #(5) 读取频率配置
    print("read frequency number...")
    with open("input/freq.inp", 'r') as fid:
        next(fid)  
        gv.nfreq = int(fid.readline().strip())
        gv.freq=np.empty(gv.nfreq,dtype=float)
        next(fid) 
        for m in range(gv.nfreq):
            gv.freq[m] = float(fid.readline().strip())
    gv.xomega=2*gv.pai*gv.freq
    gv.cer=gv.ttcer+gv.xsig/(1j*gv.xomega*gv.epsilon0)
    gv.mur=complex(gv.ttmiur)
    gv.xk = gv.xomega*np.sqrt(gv.cer*gv.epsilon0*gv.mur*gv.mu0)

  #(6) 读取散射体信息
    print("read scatterer config...")
    with open("input/scatterer.inp", 'r') as fid:
        next(fid)  
        startx,endx,startz,endz=list(map(int,fid.readline().strip().split()))
        tempscatepsr,tempscatmur,tempscatsigma=list(map(float,fid.readline().strip().split()))
    gv.scatepsr=np.full([gv.mz, gv.mx], gv.ttcer,dtype=float)
    gv.scatsigma=np.full([gv.mz, gv.mx], gv.xsig,dtype=float)
    gv.scatmur=np.full([gv.mz, gv.mx], gv.ttmiur,dtype=float)
    gv.scatepsr[startz-1:endz,startx-1:endx]=np.full([endz-startz+1, endx-startx+1], tempscatepsr,dtype=float)
    gv.scatsigma[startz-1:endz,startx-1:endx]=np.full([endz-startz+1, endx-startx+1], tempscatsigma,dtype=float)        
    gv.scatmur[startz-1:endz,startx-1:endx]=np.full([endz-startz+1, endx-startx+1], tempscatmur,dtype=float)
    gv.xkai=np.empty([gv.nfreq,gv.mz, gv.mx],dtype=complex)
    for ifre in range(gv.nfreq):
        gv.xkai[ifre,:,:]=(gv.scatepsr+gv.scatsigma/(1j*gv.xomega[ifre]*gv.epsilon0))/gv.cer[ifre]-1

   #(6) 读取前向散射配置
    print("read forward calculation configuration...")
    with open("input/forwardconfig.inp", 'r') as fid:
        next(fid)  # 跳过第一行
        gv.maxistep = int(fid.readline().strip())
        next(fid)  # 跳过第一行
        gv.maxresidualerror = float(fid.readline().strip())

    return