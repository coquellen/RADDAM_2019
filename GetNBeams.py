import numpy as np
from scipy.optimize import curve_fit



def beam_Nth(x, a0, n, a=0.5):
    return a / (2 ** n) * np.exp(- a0  * x / (2 ** n))

from matplotlib import pyplot as plt
fig, ax1 = plt.subplots(1,1, figsize=(6,10))

for choice in [0]:
    if choice == 1:
        DR = 2.1
        data = np.load('ILD_decay_04_2019.npy').T
        DATA = np.load('fitALL_LD.npy')
        X0 = np.arange(0., 180., 2) / 1000 * DR

        t0 = np.arange(0., 180., 2)

    if choice == 0:
        DR = 28.7
        data = np.load('IHD_decay_04_2019.npy').T
        DATA = np.load('fitALL_HD.npy')
        X0 = np.arange(0., 80., 2) / 1000 * DR
        t0 = np.arange(0., 80., 2)

    if choice == 2:
        DR = 28.7
        data = np.load('IC_decay_04_2019.npy').T
        DATA = np.load('fitALL_C.npy')
        X0 = np.arange(0., 180., 2) / 1000 * DR
        t0 = np.arange(0., 180., 2)

    normalized = data / data[0,:]

    N = normalized.shape[1]

    markersize=4

    y = np.nanmean(normalized,axis=1)
    err = np.nanstd(normalized, axis=1)
    X1 = X0 + 1.079 * X0[1]
    Ntot = 50
    WeigthedDose = np.zeros((X0.size, Ntot - 1))
    Popts = np.zeros((Ntot-1,))
    for N in range(1,Ntot):
        def beamn(x, a0, a=0.5, n=N):

            if n > 1: return a / (2 ** n) * np.exp(-  x / (a0 * (2 ** n))) + beamn(x, a0, a, n - 1)
            if n == 1: return a / (2 ** n) * np.exp(-  x / (a0 * (2 ** n))) + a / (2 ** n - 1) * np.exp(
                - x / (a0 * (2 ** (n - 1))))


        popt, pcov = curve_fit(beamn,X0 ,y ,p0=[1.])#, bounds=(0.01, np.inf))# 1, 1]
        #print popt
        if choice == 0:
            marker = '^'

        Popts[N-1] = popt[0]
        beams = []
        SUM = np.zeros_like(X0)
        #Nbeams = 30
        for k in range(N):
            b = 0.5 / (2 ** k) * np.exp(- (X1 / (popt[0] * 2. ** k)))
            beams.append(b)
            SUM += b

        #print beams

        for k in range(N):
            WeigthedDose[:,N-1] += beams[k] * DR * (t0 + 2) / 1000.

    #diff = np.diff(WeigthedDose,n=1, axis =1)
    diff = np.diff(Popts)
    print Popts
    print diff
    ax1.semilogy(np.arange(1,Ntot-1,1), np.abs(diff), marker='o')#,label='%i'%i)
    ax1.set_ylabel("Difference in final DWD value")
    ax1.set_xlabel("Number of beams used in the model")
    #ax1.plot(np.arange(0,29,1), WeigthedDose[-1,:])#, label='%i' % i)
    ax1.legend()
    #print diff[-1,:]
    for i in range(0,Ntot-1):
        if np.abs(diff[i]) < 1e-6:
            print i+1
            break
    #print diff[-1,20]
    #plt.hlines(1e-6,1,29)
    plt.show()
    #print WeigthedDose[:,0]
#plt.savefig('Figure3.eps', format='eps', dpi=1200)
#plt.show()
#
#     beams = []
#     SUM = np.zeros_like(X)
#     W = []
#     WeigthedDose = np.zeros_like(X)
#
#     for k in range(i):
#             b = 0.5 / (2 ** k) *  np.exp(- popt[0] * (X /2.** k))
#             beams.append(b)
#             SUM += b
#
#     for k in range(i):
#             W.append(beams[k] / SUM)
#
#     for k in range(i):
#             WeigthedDose += beams[k] * DR * t / 1000.
#
#     def SpecN(x, A, b, a0, k=1, n=i ):
# #        K = (k - (1. / (2 ** n)))
# #        if n > 1: return A * (W[n-1]  * K *   np.exp(- a0  * x / (2 ** (n-1)))) + SpecN(x, A, b, a0,  k, n=n-1)
# #        if n == 1: return A * (W[n-1] * K * np.exp( - a0 * x / (2 ** (n - 1))))  + b
#         K = 1/n
#         if n > 1: return A * (W[n - 1]  * np.exp(- K * a0 * x / (2 ** (n - 1)))) + SpecN(x, A, b, a0, k=1, n=n - 1)
#         if n == 1: return A * (W[n - 1] * np.exp(- K *a0 * x / (2 ** (n - 1)))) + b
#
#
#
#
#     popt, pcov = curve_fit(SpecN, X ,SS0,p0=[ 1, 1, 1])#,sigma=dataHDrefLD2[2, idxHD2[add1], :20])
#     print "%i %4.2f +/- %4.2f" % (i, 1. / popt[2], pcov[2][2] ** 0.5)
#     ax3.plot(X , SpecN(X , *popt),'^')
#     #ax3.plot(WeigthedDose, y, 'o')
#
#     ax2.plot(X , WeigthedDose, 'o',label='WD %i beams'%i)#, color='r',label='W1')
#
#
#
#
#
#
# ax2.set_ylabel("DWD (MGy)")
# ax1.errorbar(X , y, yerr=err, marker='o', color='k', markersize=3, linestyle='',linewidth=0.5)
# ax1.set_xlabel("Nominal Dose ()")
# ax2.set_xlabel("Nominal Dose (MGy)")
# ax3.set_xlabel("DWD (MGy)")
#
# ax1.legend()
#
# ms=5
# ls = ''
# errev = 5
#
# add = [(0,7), (1,6), (2,4), (3,5)]
# idxHD =   np.array([idx for idx, i in enumerate(idHD) if 'CYS' in i and 'SG' in i])
# idxHD2 =   np.array([idx for idx, i in enumerate(idHD2) if 'CYS' in i and 'SG' in i])
# add0, add1 = add[0]
# label = labelsCRYO[add0][0:7] +' - '+ labelsCRYO[add1][0:7]
#
#
#
# ax3.errorbar( X , SS0, marker='x', linestyle='', markersize=ms, linewidth=1, errorevery=errev,
#              yerr=dataHDrefLD2[2, idxHD2[add1], :], label=label)
#
#
# #ax3.set_ylabel("I/I1")
# ax3.legend()
#
# plt.show()
#
