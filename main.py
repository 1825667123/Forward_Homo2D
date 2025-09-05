import numpy as np
import global_para_var as gv
from ReadInput import Read_Input
from Compute_IncField import Compute_IncField
from Compute_Gii import Compute_Gii
from FFT_2DGreen import Gej_FFT
from Compute_TotField import Compute_TotField
from Compute_SctField import Compute_SctField
from SaveData import Save_IncField,Save_TotField,Save_SctField

def main():
    # 读取输入参数
    Read_Input()

    #计算入射场
    Eyinc = np.empty((gv.nfreq,gv.ntrtot,gv.mz,gv.mx),dtype=complex)
    Compute_IncField(Eyinc)
    Save_IncField(Eyinc)

    #计算自身到自身格林函数，然后实施FFT
    gejyym=np.empty((gv.nfreq,2*gv.mz-1,2*gv.mx-1),dtype=complex)
    Compute_Gii(gejyym)
    gejyymfre=np.empty(np.shape(gejyym),dtype=complex)
    Gej_FFT(gejyym, gejyymfre)

    #计算总场，通过BCGS-FFT
    Eytot=np.empty(np.shape(Eyinc),dtype=complex)
    Compute_TotField(Eyinc, Eytot, gejyymfre)
    Save_TotField(Eytot)

    #计算散射场
    Eysct=np.empty((gv.nfreq,gv.ntrtot,gv.nrectot),dtype=complex)
    Compute_SctField(Eytot, Eysct)
    Save_SctField(Eysct)

if __name__ == "__main__":
    main()