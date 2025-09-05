import numpy as np
from ReadInput import read_input1, read_input2, compute_cofe
import global_para_var


def Compute_IncField(Eyinc):

    print(f"\nCompute_IncField 函数启动:")
    print(f"  网格大小: mx={global_para_var.mx}, mz={global_para_var.mz}")
    print(f"  发射器数量: ntrtot={global_para_var.ntrtot}")
    print(f"  频率点数: nfreq={global_para_var.nfreq}")


    xxf = np.zeros(global_para_var.mx * global_para_var.mz, dtype=np.float64)
    xzf = np.zeros(global_para_var.mx * global_para_var.mz, dtype=np.float64)


    m = 0
    for k in range(global_para_var.mz):
        for i in range(global_para_var.mx):
            xxf[m] = global_para_var.xx1 + i * global_para_var.dx
            xzf[m] = global_para_var.zz1 + k * global_para_var.dz
            m += 1

    print('Compute inc field, please wait')  # 打印提示信息


    # 遍历所有频率和发射器，计算入射场
    for ifre in range(global_para_var.nfreq):
        xomega = 2 * global_para_var.pai * global_para_var.freq[ifre]  # 计算角频率
        kkk = xomega * np.sqrt(global_para_var.epsilon0 * global_para_var.mu0)  # 计算真空中波数

        for itrans in range(global_para_var.ntrtot):
            m = 0  # 重置计数器
            kkx = kkk * np.cos(global_para_var.theta_inc[itrans])  # 计算波数x分量
            kkz = kkk * np.sin(global_para_var.theta_inc[itrans])  # 计算波数z分量

            # 计算每个网格点的入射场
            for k in range(global_para_var.mz):
                for i in range(global_para_var.mx):
                    phase = -1j * (kkx * xxf[m] + kkz * xzf[m])
                    Eyinc[i, k, itrans, ifre] = np.exp(phase)
                    m += 1


if __name__ == "__main__":
    read_input1()
    read_input2()
    compute_cofe()

    Eyinc = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.ntrtot, global_para_var.nfreq), dtype=np.complex_)
    Compute_IncField(Eyinc)
    print("入射场计算完成")

    # 打印部分值

    print(Eyinc[ 1, 0, 0, 0])