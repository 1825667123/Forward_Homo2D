import numpy as np
import global_para_var
from SpatialGreenEJ_I2I import greens_function_comp_ej_i2i


def compute_gii():

    gejm = np.zeros((2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq),
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

    # 频率循环
    for ifre in range(global_para_var.nfreq):
        print(f'Compute Gii for frequency= {global_para_var.freq[ifre]}')

        xomega = 2.0 * global_para_var.pai * global_para_var.freq[ifre]

        # 更新每个层的介质参数
        for i in range(global_para_var.nlayer):
            global_para_var.mur[i] = complex(global_para_var.ttmiur[i], 0.0)
            global_para_var.cer[i] = complex(global_para_var.ttcer[i],
                                             -(global_para_var.xsig[i]) / (xomega * global_para_var.epsilon0))
            global_para_var.ck[i] = xomega * np.sqrt(global_para_var.cer[i] * global_para_var.mur[i] *
                                                     global_para_var.epsilon0 * global_para_var.mu0)
            process_ckz(i)  # 处理ckz

        # 计算格林函数
        green_output_e1st = greens_function_comp_ej_i2i(xxf, xzf, global_para_var.xx1,
                                                      global_para_var.zz1, global_para_var.mx * global_para_var.mz)
        green_output_e2nd = greens_function_comp_ej_i2i(xxf, xzf, global_para_var.xx1,
                                                      global_para_var.zz1 + (
                                                                  global_para_var.mz - 1) * global_para_var.dz,
                                                      global_para_var.mx * global_para_var.mz)

        # 组装格林函数二维数组
        gejm[:, :, ifre] = assemble_green(green_output_e1st, green_output_e2nd)

    return gejm


def assemble_green(green_e1st, green_e2nd):
    """
    将一维的格林函数数据重新排列成二维数组格式
    """
    # 初始化输出数组
    geyym = np.zeros((2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1), dtype=np.complex_)

    # 重塑输入数据
    green_e1st_2d = green_e1st.reshape((global_para_var.mx, global_para_var.mz))
    green_e2nd_2d = green_e2nd.reshape((global_para_var.mx, global_para_var.mz))

    # 第一层数据处理
    geyym[0:global_para_var.mx, 0:global_para_var.mz] = green_e1st_2d
    geyym[global_para_var.mx:2 * global_para_var.mx - 1, 0:global_para_var.mz] = \
        np.flip(green_e1st_2d[1:, :], axis=0)  # 镜像处理

    # 第二层数据处理
    geyym[0:global_para_var.mx, global_para_var.mz:2 * global_para_var.mz - 1] = \
        green_e2nd_2d[:, 0:global_para_var.mz - 1]
    geyym[global_para_var.mx:2 * global_para_var.mx - 1, global_para_var.mz:2 * global_para_var.mz - 1] = \
        np.flip(green_e2nd_2d[1:, 0:global_para_var.mz - 1], axis=0)  # 镜像处理

    return geyym


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


