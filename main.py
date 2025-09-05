import numpy as np
from ReadInput import read_input1, read_input2, compute_cofe
from Compute_IncField import Compute_IncField
from Compute_Gii import compute_gii
from Compute_TotField import Compute_TotField
from Compute_SctField import compute_sct_field
from SaveData import SavePreparedData1st, SavePreparedData2nd, SaveTruePara, SaveSctField
import global_para_var


def main():
    # 读取输入参数
    read_input1()

    # 分配数组空间并读取详细参数
    global_para_var.ttcer = np.zeros(global_para_var.nlayer + 1)
    global_para_var.ttmiur = np.zeros(global_para_var.nlayer + 1)
    global_para_var.xsig = np.zeros(global_para_var.nlayer + 1)
    global_para_var.xloc = np.zeros(global_para_var.nlayer + 1)
    global_para_var.theta_inc = np.zeros(global_para_var.ntrtot)
    global_para_var.xrr = np.zeros(global_para_var.nrectot)
    global_para_var.zrr = np.zeros(global_para_var.nrectot)
    global_para_var.freq = np.zeros(global_para_var.nfreq)
    global_para_var.scatepsr = np.zeros((global_para_var.mx, global_para_var.mz))
    global_para_var.scatmur = np.zeros((global_para_var.mx, global_para_var.mz))
    global_para_var.scatsigma = np.zeros((global_para_var.mx, global_para_var.mz))
    global_para_var.epslonbr = np.zeros(global_para_var.nfreq, dtype=complex)
    global_para_var.ck = np.zeros(global_para_var.nlayer, dtype=complex)
    global_para_var.cer = np.zeros(global_para_var.nlayer, dtype=complex)
    global_para_var.mur = np.zeros(global_para_var.nlayer, dtype=complex)
    global_para_var.xk = np.zeros(global_para_var.nfreq, dtype=complex)
    global_para_var.xkai = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.nfreq), dtype=complex)

    read_input2()

    # 保存真实参数
    SaveTruePara()

    # 计算系数
    compute_cofe()

    # 计算入射场
    Eyinc = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq),
                     dtype=np.complex_)
    Compute_IncField(Eyinc)
    print("入射场计算完成")

    # 计算格林函数Gii
    gejyymfre = compute_gii()
    print("格林函数Gii计算完成")

    # 保存准备数据1
    SavePreparedData1st(Eyinc, gejyymfre)
    print("准备数据1保存完成")

    # 计算总场
    Eytot = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq),
                     dtype=np.complex_)
    Compute_TotField(Eyinc, Eytot, gejyymfre)
    print("总场计算完成")

    # 计算格林函数Gri (这里使用零数组占位，实际需要实现Compute_Gri函数)
    gejriyy = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.nrectot, global_para_var.nfreq),
                       dtype=np.complex_)
    # TODO: 实现Compute_Gri函数

    # 保存准备数据2
    SavePreparedData2nd(gejriyy)
    print("准备数据2保存完成")

    # 计算散射场
    Eysct = compute_sct_field(Eytot, gejriyy)
    print("散射场计算完成")

    # 保存散射场
    SaveSctField(Eysct)
    print("散射场保存完成")


if __name__ == "__main__":
    main()
