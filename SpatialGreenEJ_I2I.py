import numpy as np
import global_para_var
import scipy.special as sp

def greens_function_comp_ej_i2i(xx, xz, xxs, xzs, nxx):
    green_output_logging_2m = np.zeros(nxx, dtype=np.complex128)
    dx = global_para_var.dx

    for i in range(nxx):
        xrho = abs(xx[i] - xxs)
        xrhoz = abs(xz[i] - xzs)

        if np.sqrt(xrho**2 + xrhoz**2) < 0.001 * dx:
            directw = gejs2s()
        else:
            ctx = global_para_var.ck[global_para_var.layerofscat] * np.sqrt(xrho**2 + xrhoz**2)
            cbj, ch = bslsdr(2, ctx, 2)
            directw = -global_para_var.xomega * global_para_var.mur[global_para_var.layerofscat] * \
                      global_para_var.mu0 / 4.0 * (2.0 * cbj[0] - ch[0])

        green_output_logging_2m[i] = directw

    return green_output_logging_2m


def gejs2s():
    a = np.sqrt(global_para_var.dx * global_para_var.dz / np.pi)
    directw = global_para_var.xomega * global_para_var.mur[global_para_var.layerofscat] * \
              global_para_var.mu0 / 4.0 * ((0.0 + 1.0j) / np.pi * (2 * np.log(0.8905 * global_para_var.ck[global_para_var.layerofscat] * a) - 1) - 1)

    return directw


def bslsdr(order, ctx, mode):
    """
    计算0阶贝塞尔函数和汉克尔函数及其导数（针对实际使用需求优化）

    参数:
    order: 阶数
    ctx: 变量参数
    mode: 模式

    返回:
    cbj: 贝塞尔函数值数组 [J_0(x)]
    ch: 汉克尔函数值数组 [H_0(x)]
    cjp: 贝塞尔函数导数值数组 [J'_0(x)]
    chp: 汉克尔函数导数值数组 [H'_0(x)]
    """
    # 初始化
    cbj = np.zeros(1, dtype=np.complex128)
    ch = np.zeros(1, dtype=np.complex128)
    cjp = np.zeros(1, dtype=np.complex128)
    chp = np.zeros(1, dtype=np.complex128)

    # 贝塞尔函数和汉克尔函数
    n = 0
    # 贝塞尔函数
    cbj[0] = sp.jv(n, ctx)
    # 汉克尔函数
    # ch[0] = sp.jv(0, ctx) + sp.yv(0, ctx) * 1j
    ch[0] = sp.hankel1(n, ctx)

    # # 贝塞尔函数导数
    # cjp[0] = -sp.jv(1, ctx)
    # # 汉克尔函数导数
    # chp[0] = -sp.hankel1(1, ctx)

    return cbj, ch

