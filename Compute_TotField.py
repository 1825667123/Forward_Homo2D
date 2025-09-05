import numpy as np
from scipy.fft import ifftn, fftshift, ifftshift, fftn
import global_para_var
from ReadInput import read_input1, read_input2, compute_cofe, read_gejyymfre_from_file
from Compute_IncField import Compute_IncField
def Compute_TotField(Eyinc, Eytot, gejyymfre):

    for ifre in range(global_para_var.nfreq):
        for itrans in range(global_para_var.ntrtot):
            print(
                f"Perform the BCGSFFT iteration for frequency {ifre + 1} and transmitter {itrans + 1}, please wait...")
            BCGSFFT(
                Eyinc[:, :, itrans, ifre],
                Eytot[:, :, itrans, ifre],
                gejyymfre[:, :, ifre],
                ifre,
                itrans
            )


def BCGSFFT(Eyinc, Eytot, gejyymfre, ifre, itrans):

    rho = complex(1.0, 0.0)
    beta = complex(1.0, 0.0)
    dabu = complex(1.0, 0.0)
    rho1 = complex(0.0, 0.0)
    c1, c2 = complex(0.0, 0.0), complex(0.0, 0.0)

    itabnormal = 500
    memm = 0
    rerror0 = 1.0e+09
    rerror = 0.0
    rnormb = 0.0

    # 初始化数组
    xr0 = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)
    xrn = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)
    xpn = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)
    xvn = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)
    xxn = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)
    xtn = np.zeros((global_para_var.mx, global_para_var.mz), dtype=np.complex_)

    # 设置初始值
    xr0[:] = Eyinc[:]
    xrn[:] = xr0[:]

    rnormb = np.sum(np.conj(xr0) * xr0).real
    print(f"rnornmb= {rnormb:.3f}   {global_para_var.maxresidualerror:.7E}")

    # 迭代求解
    for istep in range(1, global_para_var.maxistep + 1):
        rho1 = np.sum(np.conj(xr0) * xrn)

        beta = (rho1 / rho) * (beta / dabu)

        xpn[:] = xrn + beta * (xpn - dabu * xvn)

        ComputeLJ(
            ifre=ifre,
            gejyymfre=gejyymfre,
            Eyinc=xvn,
            Eytot=xpn
        )

        beta = np.sum(np.conj(xr0) * xvn)
        beta = rho1 / beta

        xrn[:] = xrn - beta * xvn

        rnorms = np.sum(np.conj(xrn) * xrn).real

        ComputeLJ(
            ifre=ifre,
            gejyymfre=gejyymfre,
            Eyinc=xtn,
            Eytot=xrn
        )

        c1 = np.sum(np.conj(xtn) * xrn)
        c2 = np.sum(np.conj(xtn) * xtn)

        if c2 == complex(0.0, 0.0):
            dabu = complex(1.0, 0.0)
        else:
            dabu = c1 / c2

        xxn[:] = xxn + beta * xpn + dabu * xrn
        xrn[:] = xrn - dabu * xtn
        rerror = np.sqrt(rnorms / rnormb)

        if rerror < rerror0:
            memm += 1
            rerror0 = rerror

        print(f"{istep:5d} {istep:5d} {rerror:.7f}")

        if rerror <= global_para_var.maxresidualerror or memm >= itabnormal:
            break

    Eytot[:] = xxn[:]


def ComputeLJ(ifre, gejyymfre, Eyinc, Eytot):

    Jymfre = np.zeros((2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1), dtype=np.complex_)


    JFFT(Eytot, Jymfre, ifre)


    copen = Jymfre * gejyymfre


    copen = ifftshift(copen)
    copen = ifftn(copen, axes=(0, 1))
    copen = fftshift(copen)


    Eysct = copen[0:global_para_var.mx, 0:global_para_var.mz] * global_para_var.acx * global_para_var.acz * global_para_var.dxz


    Eyinc[:] = Eytot - Eysct


def JFFT(Eytot, Jymfre, ifre):


    extended_Eytot = np.zeros((2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1), dtype=np.complex_)
    extended_Eytot[0:global_para_var.mx, 0:global_para_var.mz] = Eytot


    extended_Eytot = ifftshift(extended_Eytot)
    Jymfre[:] = fftn(extended_Eytot, axes=(0, 1))
    Jymfre = fftshift(Jymfre)


    Jymfre[:] = Jymfre * complex(0.0, 1.0) * 2.0 * global_para_var.pai * global_para_var.freq[ifre]



if __name__ == "__main__":

    read_input1()
    read_input2()
    compute_cofe()

    Eyinc = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq),dtype=np.complex_)
    Eytot = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq), dtype=np.complex_)
    gejyymfre = read_gejyymfre_from_file()

    Compute_IncField(Eyinc)
    print("入射场计算完成,Eyinc[ 3, 0, 0, 0]:")
    print(Eyinc[ 3, 0, 0, 0])


    Compute_TotField(Eyinc, Eytot, gejyymfre)

    print("入射场计算完成,Eyinc[ 3, 0, 0, 0]:")
    print(Eyinc[ 0, 0, 0, 0])
    print(Eyinc[ 1, 0, 0, 0])
    print(Eyinc[ 2, 0, 0, 0])
    print(Eyinc[ 3, 0, 0, 0])
    print(Eyinc[ 4, 0, 0, 0])

    print("总场计算完成,Eytot[ 3, 0, 0, 0]:")
    print(Eytot[ 0, 0, 0, 0])
    print(Eytot[ 1, 0, 0, 0])
    print(Eytot[ 2, 0, 0, 0])
    print(Eytot[ 3, 0, 0, 0])
    print(Eytot[ 4, 0, 0, 0])



