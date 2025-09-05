# ============================== 常用常数参数 ==========================
pai = 3.1415926535897931  # π
mu0 = 4.0 * pai * 1e-07    # 真空磁导率
c_vac = 2.99792458e8       # 真空中光速
epsilon0 = 1.0 / (c_vac**2 * mu0)  # 真空介电常数

# ============================== 网格参数 ==============================
mx = None  # x方向的cell个数
mz = None  # z方向的cell个数
xx1 = None  # 第一个网格中心的x坐标
zz1 = None  # 第一个网格中心的z坐标
dx = None   # x方向的网格步长
dz = None   # z方向的网格步长
xxf = None  # 网格中心点位置
xzf = None  # 网格中心点位置

# ============================== 频率参数 ==============================
nfreq = None    # 频率点的个数
freq=None   #频率值（real数组）
xomega=None  #角度频率值（real数组）

# ============================== 发射机参数 ==============================
ntrtot = None  # 发射点源的个数
xptr= None # 发射点的x坐标（real数组）
zptr= None # 发射点的z坐标（real数组）

# ============================== 接收机参数 ==============================
nrectot = None  # 发射点源的个数
xrr= None # 接收点的x坐标（real数组）
zrr= None # 接收点的z坐标（real数组）

# ============================背景介质参数 ==============================
ttcer = None        # 背景的相对介电常数（real数组）
ttmiur = None      # 每层的相对磁导率（real数组）
xsig = None        # 每层的电导率（real数组）
cer = None         # 当前频率下的复数相对介电常数（complex数组）
mur = None        # 当前频率下的复数相对磁导率（complex数组）

# ============================散射体参数 ==============================
startx = None       # 散射体x起始位置
endx = None         # 散射体x结束位置
startz = None       # 散射体z起始位置
endz = None         # 散射体z结束位置

scatepsr = None            # 散射体的相对介电常数（real数组，mx×mz）
scatmur = None             # 散射体的相对磁导率（real数组，mx×mz）
scatsigma = None           # 散射体的电导率（real数组，mx×mz）

xkai = None                # 对比度（complex数组，mx×mz×nfreq）
xk = None                  # 背景波数（complex数组，nfreq）

# ============================== BCGS迭代参数 ========================
maxistep = None               # 最大迭代步数
maxresidualerror = None       # 残差最大值

# ============================== FFT加速参数 ==========================
dxz = None      # 网格面积
#acx = None      # 归一化系数 不再需要，python中的ifft自带这个系数了
#acz = None      # 归一化系数
