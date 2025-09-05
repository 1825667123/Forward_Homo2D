import numpy as np
import global_para_var


def read_input1():

    print("开始读取输入参数...")

    # (1) 读取分层数量
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\LayerConfig.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.nlayer = int(f.readline().strip())
            print(f"成功读取分层数量: nlayer={global_para_var.nlayer}")
    except Exception as e:
        print(f"读取分层数量失败: {e}")
        return False

    # (2) 读取发射点数量
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\TLoc.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.ntrtot = int(f.readline().strip())
            print(f"成功读取发射点数量: ntrtot={global_para_var.ntrtot}")
    except Exception as e:
        print(f"读取发射点数量失败: {e}")
        return False

    # (3) 读取接收点数量
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\RLoc.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.nrectot = int(f.readline().strip())
            print(f"成功读取接收点数量: nrectot={global_para_var.nrectot}")
    except Exception as e:
        print(f"读取接收点数量失败: {e}")
        return False

    # (4) 读取网格数量
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\Cell.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.mx, global_para_var.mz = map(int, f.readline().split())
            print(f"成功读取网格数量: mx={global_para_var.mx}, mz={global_para_var.mz}")
    except Exception as e:
        print(f"读取网格数量失败: {e}")
        return False

    # (5) 读取频率点数量
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\Freq.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.nfreq = int(f.readline().strip())
            print(f"成功读取频率点数量: nfreq={global_para_var.nfreq}")
    except Exception as e:
        print(f"读取频率点数量失败: {e}")
        return False

    # (6) 读取正演计算配置
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\ForwardConfig.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.maxistep = int(f.readline().strip())
            next(f)  # 跳过空行
            global_para_var.maxresidualerror = float(f.readline().strip())
            next(f)  # 跳过空行
            global_para_var.mxyunit = float(f.readline().strip())
            next(f)  # 跳过空行
            global_para_var.mxymax = float(f.readline().strip())
            print("成功读取正演计算配置")
    except Exception as e:
        print(f"读取正演计算配置失败: {e}")
        return False

    return True  # 所有参数读取成功


def read_input2():
    """读取详细仿真参数（分层信息、收发位置、频率值等）"""
    print("开始读取详细参数")

    # (1) 读取分层配置
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\LayerConfig.inp", 'r') as f:
            for _ in range(3):  # 跳过前三行
                next(f)
            print('The (epsr,sigma,mur,boundary location) for each layer are:')
            global_para_var.ttcer = np.zeros(global_para_var.nlayer + 1)  # +1 是为了与Fortran索引保持一致
            global_para_var.xsig = np.zeros(global_para_var.nlayer + 1)
            global_para_var.ttmiur = np.zeros(global_para_var.nlayer + 1)
            global_para_var.xloc = np.zeros(global_para_var.nlayer + 1)

            for i in range(1, global_para_var.nlayer):  # 1到nlayer-1
                line = f.readline().split(',')
                global_para_var.ttcer[i] = float(line[0])
                global_para_var.xsig[i] = float(line[1])
                global_para_var.ttmiur[i] = float(line[2])
                global_para_var.xloc[i] = float(line[3])
                print(i, global_para_var.ttcer[i], global_para_var.xsig[i], global_para_var.ttmiur[i],
                      global_para_var.xloc[i])

            # 最后一层（没有位置信息）
            line = f.readline().split(',')
            global_para_var.ttcer[global_para_var.nlayer] = float(line[0])
            global_para_var.xsig[global_para_var.nlayer] = float(line[1])
            global_para_var.ttmiur[global_para_var.nlayer] = float(line[2])
            print(global_para_var.nlayer, global_para_var.ttcer[global_para_var.nlayer],
                  global_para_var.xsig[global_para_var.nlayer], global_para_var.ttmiur[global_para_var.nlayer])
    except Exception as e:
        print(f"读取分层配置失败: {e}")
        return False

    # (2) 初始化发射点角度（均匀分布）
    try:
        global_para_var.theta_inc = np.zeros(global_para_var.ntrtot)
        for m in range(global_para_var.ntrtot):
            global_para_var.theta_inc[m] = (2 * global_para_var.pai / global_para_var.ntrtot) * m
        print(f"发射角度计算完成")
    except Exception as e:
        print(f"发射角度计算失败: {e}")
        return False

    # (3) 读取接收点位置
    try:
        print('Read the location of receivers')
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\RLoc.inp", 'r') as f:
            for _ in range(3):  # 跳过前三行
                next(f)
            global_para_var.xrr = np.zeros(global_para_var.nrectot)
            global_para_var.zrr = np.zeros(global_para_var.nrectot)
            for m in range(global_para_var.nrectot):
                global_para_var.xrr[m], global_para_var.zrr[m] = map(float, f.readline().split())
    except Exception as e:
        print(f"读取接收点位置失败: {e}")
        return False

    # (4) 读取网格位置和步长
    try:
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\Cell.inp", 'r') as f:
            for _ in range(3):  # 跳过前三行
                next(f)
            global_para_var.xx1, global_para_var.zz1 = map(float, f.readline().split())
            print(f'The 1st cell center location is {global_para_var.xx1}, {global_para_var.zz1}')
            next(f)  # 跳过空行
            global_para_var.dx, global_para_var.dz = map(float, f.readline().split())
            print(f'The cell increasing steps (dx dz) are {global_para_var.dx}, {global_para_var.dz}')
    except Exception as e:
        print(f"读取网格位置和步长失败: {e}")
        return False

    # (5) 读取频率值
    try:
        print('read frequency')
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\Freq.inp", 'r') as f:
            for _ in range(3):  # 跳过前三行
                next(f)
            global_para_var.freq = np.zeros(global_para_var.nfreq)
            for m in range(global_para_var.nfreq):
                global_para_var.freq[m] = float(f.readline().strip())
    except Exception as e:
        print(f"读取频率值失败: {e}")
        return False

    # (6) 设置散射体信息
    try:
        print('Read the information of scatterer')
        with open(r"D:\Users\guanxin's pc\PycharmProjects\全区域\InPut\Scatterer.inp", 'r') as f:
            next(f)  # 跳过第一行
            global_para_var.layerofscat = int(f.readline().strip())
            print(f'The scatterer is in layer {global_para_var.layerofscat}')
            next(f)  # 跳过空行
            tempscatepsr, tempscatmur, tempscatsigma = map(float, f.readline().split())
            next(f)  # 跳过空行

            # 初始化散射体参数
            global_para_var.scatepsr = np.full((global_para_var.mx, global_para_var.mz),
                                               global_para_var.ttcer[global_para_var.layerofscat])
            global_para_var.scatmur = np.full((global_para_var.mx, global_para_var.mz),
                                              global_para_var.ttmiur[global_para_var.layerofscat])
            global_para_var.scatsigma = np.full((global_para_var.mx, global_para_var.mz),
                                                global_para_var.xsig[global_para_var.layerofscat])

            # 读取散射体区域
            for m in range(global_para_var.mz):
                for i in range(global_para_var.mx):
                    bit = int(f.readline().strip())
                    if bit == 1:
                        global_para_var.scatepsr[i, m] = tempscatepsr
                        global_para_var.scatsigma[i, m] = tempscatsigma
    except Exception as e:
        print(f"读取散射体信息失败: {e}")
        return False

    return True


def compute_cofe():

    print("开始compute_cofe")

    try:
        global_para_var.xloc[0] = global_para_var.xloc[1] - 20

        # 计算背景介电常数（考虑电导率的复数形式）
        global_para_var.epslonbr = np.zeros(global_para_var.nfreq, dtype=complex)
        for iFre in range(global_para_var.nfreq):
            global_para_var.epslonbr[iFre] = global_para_var.ttcer[global_para_var.layerofscat] - \
                                             1j * global_para_var.xsig[global_para_var.layerofscat] / \
                                             (2 * global_para_var.pai * global_para_var.freq[
                                                 iFre] * global_para_var.epsilon0)

        # 计算对比度和波数
        global_para_var.xkai = np.zeros((global_para_var.mx, global_para_var.mz, global_para_var.nfreq), dtype=complex)
        global_para_var.xk = np.zeros(global_para_var.nfreq, dtype=complex)
        for iFre in range(global_para_var.nfreq):
            global_para_var.xkai[:, :, iFre] = (
                                                           global_para_var.scatepsr - 1j * global_para_var.scatsigma /
                                                           (2 * global_para_var.pai * global_para_var.freq[
                                                               iFre] * global_para_var.epsilon0)) / \
                                               global_para_var.epslonbr[iFre] - 1
            global_para_var.xk[iFre] = 2 * global_para_var.pai * global_para_var.freq[iFre] * np.sqrt(
                global_para_var.epslonbr[iFre] * global_para_var.ttmiur[global_para_var.layerofscat] *
                global_para_var.epsilon0 * global_para_var.mu0)

        # 计算网格面积和归一化系数
        global_para_var.dxz = global_para_var.dx * global_para_var.dz
        global_para_var.acx = 1.0 / (2 * global_para_var.mx - 1)
        global_para_var.acz = 1.0 / (2 * global_para_var.mz - 1)

        print("compute_cofe完成")
        return True
    except Exception as e:
        print(f"compute_cofe失败: {e}")
        return False


def read_gejyymfre_from_file(filename=r'D:\桌面\全区域\全区域\Forward\OutPut\gejyymfre.bin'):
    """
    从二进制文件中读取gejyymfre数据
    """
    try:
        # 计算期望的数据大小和形状
        expected_size = (2 * global_para_var.mx - 1) * (2 * global_para_var.mz - 1) * global_para_var.nfreq
        expected_shape = (2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq)

        # 读取二进制数据
        with open(filename, 'rb') as f:
            file_content = f.read()

        # 使用complex64类型（32位复数）解析数据
        gejyymfre_data = np.frombuffer(file_content, dtype=np.complex64)

        # 验证数据大小
        if len(gejyymfre_data) != expected_size:
            raise ValueError(f"数据大小不匹配: 预期{expected_size}，实际{len(gejyymfre_data)}")

        # 重塑为期望的形状
        gejyymfre = gejyymfre_data.reshape(expected_shape, order='F')

        print(f"成功从 {filename} 读取 gejyymfre 数据")
        return gejyymfre

    except FileNotFoundError:
        print(f"警告: 文件 {filename} 不存在，使用随机数据")
        return np.random.rand(2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq) + \
            1j * np.random.rand(2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq)
    except Exception as e:
        print(f"读取文件时出错: {e}")
        print("使用随机数据作为替代")
        return np.random.rand(2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq) + \
            1j * np.random.rand(2 * global_para_var.mx - 1, 2 * global_para_var.mz - 1, global_para_var.nfreq)
