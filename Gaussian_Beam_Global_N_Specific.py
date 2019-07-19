import os, sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


data = np.loadtxt(sys.argv[1])
shape = data.shape

if shape[1] == 2:
    print("Only the global damage will be analyzed")
elif shape[1] == 3:
    print("Both global and specific damage will be analyzed")
else:
    exit()



ms=3
ls = ''
lw= 1.5
errev = 5



fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(3,8))

for i in [20]:

    def beamn(x, a0, a, n=i):

        if n > 1: return a / (2 ** n) * np.exp(- a0  * x / (2 ** n)) + beamn(x, a0, a, n - 1)
        if n == 1: return a / (2 ** n) * np.exp(- a0  * x / (2 ** n)) + a / (2 ** n - 1) * np.exp( - a0 * x / (2 ** (n - 1)))


    popt, pcov = curve_fit(beamn,X[:-12] ,y[:-12] ,p0=[1., 1])#, bounds=(0.01, np.inf))# 1, 1]
    print 1/ popt[0]
    b0 = popt[0]
    ax1.set_ylabel("I/I1")
    ax1.plot(Dose0,  beamn(Dose0 ,*popt), color='k', linewidth = lw)

    beams = []
    SUM = np.zeros_like(Dose0)
    W = []
    WeigthedDose = np.zeros_like(Dose0)

    for k in range(i):
            b = 0.5 / (2 ** k) *  np.exp(- popt[0] * (Dose0/2.** k))
            beams.append(b)
            SUM += b

    for k in range(i):
            W.append(beams[k] / SUM)

    for k in range(i):
            WeigthedDose += beams[k] * Dose0 #/ (2 ** (k+1))

    print W[1]
    def SpecN(x, A, b, a0, k=1, n=i):
        if n > 1: return A * (W[n - 1] * np.exp( a0 * x / (2 ** (n - 1)))) + SpecN(x, A, b, a0, k=1, n=n - 1)
        if n == 1: return A * (W[n - 1] * np.exp( a0 * x / (2 ** (n - 1)))) + b

    popt, pcov = curve_fit(SpecN, Dose0, SS0, p0=[1, 1, 5])  # ,sigma=dataHDrefLD2[2, idxHD2[add1], :20])
    print "%i %4.2f +/- %4.2f" % (i, 1. / popt[2], pcov[2][2] ** 0.5)
    ax2.plot(Dose0,WeigthedDose, 's',color='k', markersize=ms)
    ax3.plot(Dose0, SpecN(Dose0, *popt),color='k', linewidth = lw)


# ax1 I/I1 versus nominal dose
ax1.errorbar(X[:-12] , y[:-12], marker='s', color='k', markersize=ms, linestyle='',linewidth=0.5)
ax1.set_ylim([0,1.05])

ax1.set_xlabel("Nominal Dose (MGy)")
ax2.set_xlabel("Nominal Dose (MGy)")
ax2.set_ylabel("DWD (MGy)")

ax3.set_xlabel("Nominal Dose (MGy)")
ax3.set_ylabel("I/I1")

ax1.legend()



ax3.plot(Dose0,SS0, 's',color='k', markersize=ms)
ax3.set_ylabel("Integrated density (e-)")

print Dose0.shape
DATA = np.zeros((Dose0.shape[0] , 2))
print DATA.shape
DATA[:,0] = Dose0
DATA[:,1] = SpecN(Dose0, *popt)

np.savetxt("DATA_FIT_CRYO_%s.txt"%Z, DATA)
plt.savefig("MultiBeam_Cryo.eps",dpi=300)
plt.show()

