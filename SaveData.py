import numpy as np
from global_para_var import mx, mz, scatepsr,scatsigma,nfreq,ntrtot,nrectot

def SavePreparedData1st(Eyinc, gejyymfre):
    """保存入射电场和格林函数频域数据（二进制格式）"""
    print('Save Eyinc and Gii')  # 打印提示信息，与Fortran一致

    # 以二进制写模式打开文件
    with open('../OutPut/PreparedData1st.bin', 'wb') as f:
        # Fortran数组默认按"列优先"存储，NumPy默认"行优先"，需用order='F'展平数组以匹配二进制格式
        # 写入Eyinc：展平为1维字节流（按列优先）
        f.write(Eyinc.ravel(order='F').tobytes())
        # 写入gejyymfre：同样按列优先展平
        f.write(gejyymfre.ravel(order='F').tobytes())
    # with语句自动关闭文件，无需手动close


def SavePreparedData2nd(gejriyy):
    """保存格林函数gejriyy（二进制格式）"""
    print('Save Gri')  # 打印提示信息

    # 二进制写模式打开文件
    with open('../OutPut/PreparedData2nd.bin', 'wb') as f:
        # 按列优先展平数组（匹配Fortran存储顺序），转换为字节流写入
        f.write(gejriyy.ravel(order='F').tobytes())
    # 自动关闭文件


def SaveTruePara():
    """保存介电常数和电导率到文本文件"""
    filename_epsr = '../OutPut/True_Epsr.dat'
    filename_sigma = '../OutPut/True_Sigma.dat'

    # 同时打开两个文件（'w'为文本写模式），with语句确保自动关闭
    with open(filename_epsr, 'w') as f_epsr, open(filename_sigma, 'w') as f_sigma:
        # Fortran循环顺序：先z（iz=1到mz），再x（ix=1到mx），Python对应0基索引（0到mz-1，0到mx-1）
        for iz in range(mz):  # z方向循环（原iz=1→Python iz=0）
            for ix in range(mx):  # x方向循环（原ix=1→Python ix=0）
                # 提取介电常数（scatepsr是复数，取实部，原Fortran输出默认实部）
                epsr_val = scatepsr[ix, iz].real
                # 提取电导率（同理取实部）
                sigma_val = scatsigma[ix, iz].real
                # 写入文件，格式与Fortran默认输出一致（科学计数法）
                f_epsr.write(f'{epsr_val:.6e}\n')
                f_sigma.write(f'{sigma_val:.6e}\n')


def SaveSctField(Eysct):
    """按频率保存散射场到文本文件"""
    print('Save Eysct')  # 打印提示信息

    # 频率循环：原ifre=1到nfreq→Python ifre=0到nfreq-1（+1对应原序号）
    for ifre in range(nfreq):
        # 生成频率序号字符串：原Fortran用I5格式（5位），Python补零为5位（如1→"00001"）
        tempstr = f'{ifre + 1:05d}'  # ifre+1对应原ifre=1
        # 拼接文件名（与原路径一致）
        filename = f'../OutPut/Escty_Fre{tempstr}.dat'

        # 打开当前频率的文件（文本写模式）
        with open(filename, 'w') as f:
            # 入射角循环：原itrans=1到ntrtot→Python itrans=0到ntrtot-1
            for itrans in range(ntrtot):
                # 接收点循环：原irece=1到nrectot→Python irece=0到nrectot-1
                for irece in range(nrectot):
                    # 提取散射场值（复数）
                    sct_val = Eysct[irece, itrans, ifre]
                    # 按原格式输出：6个科学计数法数值（实部+虚部，宽度11，4位小数）
                    # 原格式(6(1x,e11.4))，这里输出实部和虚部（各占1个位置，其余补空格）
                    f.write(f'{sct_val.real:11.4e} {sct_val.imag:11.4e}      \n')