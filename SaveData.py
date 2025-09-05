import global_para_var as gv

def Save_IncField(Eyinc):

    for ifre in range(gv.nfreq):
        tempstr = f'{ifre + 1:03d}'
        filename = f'output/Eyinc_Fre{tempstr}.dat'
        with open(filename, 'w') as fid:
            for itrans in range(gv.ntrtot):
                for k in range(gv.mz):
                    for i in range(gv.mx):
                        fid.write(f'{Eyinc[ifre,itrans,k,i].real:11.7e} {Eyinc[ifre,itrans,k,i].imag:11.7e}      \n')
    return

def Save_TotField(Eytot):

    for ifre in range(gv.nfreq):
        tempstr = f'{ifre + 1:03d}'
        filename = f'output/Eytot_Fre{tempstr}.dat'
        with open(filename, 'w') as fid:
            for itrans in range(gv.ntrtot):
                for k in range(gv.mz):
                    for i in range(gv.mx):
                        fid.write(f'{Eytot[ifre,itrans,k,i].real:11.7e} {Eytot[ifre,itrans,k,i].imag:11.7e}      \n')
    return


def Save_SctField(Eysct):

    for ifre in range(gv.nfreq):
        tempstr = f'{ifre + 1:03d}'
        filename = f'output/Eysct_Fre{tempstr}.dat'
        with open(filename, 'w') as fid:
            for itrans in range(gv.ntrtot):
                for irece in range(gv.nrectot):
                    fid.write(f'{gv.xptr[itrans]} {gv.zptr[itrans]} {gv.xrr[irece]} {gv.zrr[irece]} {Eysct[ifre,itrans,irece].real:11.4e} {Eysct[ifre,itrans,irece].imag:11.7e}      \n')