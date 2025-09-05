import numpy as np
import global_para_var
from SpatialGreenEHJ_I2R import greens_function_comp_ehj  # 需要实现这个模块


def compute_gri():
    """
    计算接收点格林函数

    返回:
    gejriyy: 格林函数 (mx, mz, nrectot, nfreq)
    """

    # 初始化格林函数数组
    gejriyy = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.nrectot, global_para_var.nfreq),
                       dtype=np.complex_)

    # 创建观测点坐标数组
    xxf = np.zeros(global_para_var.mx * global_para_var.mz)
    xzf = np.zeros(global_para_var.mx * global_para_var.mz)

    m = 0
    for k in range(global_para_var.mz):
        for i in range(global_para_var.mx):
            xxf[m] = global_para_var.xx1 + i * global_para_var.dx
            xzf[m] = global_para_var.zz1 + k * global_para_var.dz
            m += 1

    print('Compute Gri, please wait')

    # 频率循环
    for ifre in range(global_para_var.nfreq):
        print(f'Compute Gri for frequency= {global_para_var.freq[ifre]}')

        global_para_var.xomega = 2.0 * global_para_var.pai * global_para_var.freq[ifre]

        # 更新每个层的介质参数
        for i in range(global_para_var.nlayer):
            global_para_var.mur[i] = complex(global_para_var.ttmiur[i], 0.0)
            global_para_var.cer[i] = complex(global_para_var.ttcer[i],
                                             -(global_para_var.xsig[i]) /
                                             (global_para_var.xomega * global_para_var.epsilon0))
            global_para_var.ck[i] = global_para_var.xomega * np.sqrt(global_para_var.cer[i] * global_para_var.mur[i] *
                                                                     global_para_var.epsilon0 * global_para_var.mu0)
            process_ckz(i)  # 处理ckz

        # 接收点循环
        for irece in range(global_para_var.nrectot):
            # 计算格林函数 (需要实现greens_function_comp_ehj函数)
            green_output_e = greens_function_comp_ehj(xxf, xzf, global_para_var.xrr[irece],
                                                      global_para_var.zrr[irece],
                                                      global_para_var.mx * global_para_var.mz)

            # 将结果存储到gejriyy中
            m = 0
            for k in range(global_para_var.mz):
                for i in range(global_para_var.mx):
                    gejriyy[i, k, irece, ifre] = green_output_e[m]
                    m += 1

    return gejriyy


def process_ckz(layer_index):
    """
    处理波数的z分量
    """
    # 获取当前层的ckz值
    ckz = global_para_var.ck[layer_index]

    # 获取ckz的虚部
    xikz = np.imag(ckz)

    # 如果虚部大于0，则取负值
    if xikz > 0:
        global_para_var.ck[layer_index] = -1 * ckz
