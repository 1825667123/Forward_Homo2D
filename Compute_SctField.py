import numpy as np
import global_para_var


def compute_sct_field(Eytot, gejriyy):
    """
    计算散射场

    参数:
    Eytot: 总场 (mx, mz, ntrtot, nfreq)
    gejriyy: 格林函数 (mx, mz, nrectot, nfreq)

    返回:
    Eysct: 散射场 (nrectot, ntrtot, nfreq)
    """

    # 初始化电流密度Jy
    Jy = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq),
                  dtype=np.complex_)

    # 计算电流密度 Jy = i*2*π*freq*ε_b*ε0*χ*Eytot
    for ifre in range(global_para_var.nfreq):
        for itrans in range(global_para_var.ntrtot):
            Jy[:, :, itrans, ifre] = (1j * 2 * global_para_var.pai * global_para_var.freq[ifre] *
                                      global_para_var.epslonbr[ifre] * global_para_var.epsilon0 *
                                      global_para_var.xkai[:, :, ifre] *
                                      Eytot[:, :, itrans, ifre])

    print('Compute scattered field, please wait')

    # 初始化散射场
    Eysct = np.zeros((global_para_var.nrectot, global_para_var.ntrtot, global_para_var.nfreq),
                     dtype=np.complex_)

    # 计算散射场 Eysct = Σ(gejriyy * Jy) * dxz
    for ifre in range(global_para_var.nfreq):
        for itrans in range(global_para_var.ntrtot):
            for irece in range(global_para_var.nrectot):
                Eysct[irece, itrans, ifre] = (Eysct[irece, itrans, ifre] +
                                              np.sum(gejriyy[:, :, irece, ifre] *
                                                     Jy[:, :, itrans, ifre]) *
                                              global_para_var.dxz)

    return Eysct
