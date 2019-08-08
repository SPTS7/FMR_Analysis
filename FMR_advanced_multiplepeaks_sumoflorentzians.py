#%%
import math, scipy, pylab
import glob, os
import matplotlib.pyplot as plt
import re, operator
import time
from lmfit.models import LorentzianModel, Model, ExponentialModel
from numpy import loadtxt
from tqdm import tqdm
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from scipy import signal


#%%

## Material variables

G = 3e6
Ms = 1222
Hu = 0
Ha = 0
dH0 = 0  # dH[0]
alpha = 0.005
erro = 1

## Choosing Anisotropy (In-Plane (IP) or Out-of-Plane (OP) )

Anisotropy = "IP"

## max errors of fit

erromaxlorentzfwhm = 0.1
erromaxlorentzpeak = 0.01
erromaxlorentzamp = 0.1
maxphase = 10


## Program variables

peakmag = []
fwhmmag = []
dH = []
dHa = []
field = []
fieldsmag = []
name1 = []
modsup = LorentzianModel()
zeta = []
paramos = []


## Opening files and getting field values from filenames

names = glob.glob("*.dat")

arr = [re.split(r"(\d+)", s) for s in names]
arr.sort(key=lambda x: int(x[1]))

for x in arr:
    field.append(int(x[-2]))

## Fitting each graphic to a lorentzian


def Lorexponential(x, simetricL, asymetricL, center, sigma, C):
    return (
        simetricL * (sigma ** 2 / ((x - center) ** 2 + sigma ** 2))
        + asymetricL * ((sigma * (x - center)) / ((x - center) ** 2 + sigma ** 2))
        + C
    )


def make_model(num):
    pref = "f{0}_".format(num)

    model = Model(Lorexponential, prefix=pref)

    model.set_param_hint(pref + "simetricL", value=0.01)
    model.set_param_hint(pref + "asymetricL", value=0.01)
    model.set_param_hint(pref + "center", value=paramos[num]["center"])
    model.set_param_hint(pref + "sigma", value=paramos[num]["sigma"])
    model.set_param_hint(pref + "C", value=0.001)

    return model


nn = 0

for files in tqdm(names, ncols=100):
    # print(files)

    paramos = []
    data = loadtxt(files)
    x = data[:, 0]
    y = data[:, 1]
    y = -y
    # y[scipy.where(y<0)]=0

    """ find peaks and determine number of peaks"""
    # Meus dados
    peakind = signal.find_peaks(y, width=2, height=0.004)
    # Dados Ana
    # peakind = signal.find_peaks(y , width=10 , height= 0.01 )

    # print(peakind[0])
    wps = x[peakind[0]]

    for pp in range(len(peakind[0])):

        larguraajuste = 2

        ypeak = y[peakind[0][pp] - larguraajuste : peakind[0][pp] + larguraajuste]
        xpeak = x[peakind[0][pp] - larguraajuste : peakind[0][pp] + larguraajuste]

        paramos.append(modsup.guess(ypeak, x=xpeak))

    mod = None
    for i in range(len(peakind[0])):
        this_mod = make_model(i)
        if mod is None:
            mod = this_mod
        else:
            mod = mod + this_mod

    outy = mod.fit(y, x=x)
    # print(outy.fit_report())

    for ind in range(len(peakind[0])):
        try:
            peakmag[ind].append(outy.best_values["f" + str(ind) + "_center"])
            fwhmmag[ind].append(outy.best_values["f" + str(ind) + "_sigma"])
            fieldsmag[ind].append(field[nn])
        except:
            peakmag.append([])
            fwhmmag.append([])
            fieldsmag.append([])
            peakmag[ind].append(outy.best_values["f" + str(ind) + "_center"])
            fwhmmag[ind].append(2 * outy.best_values["f" + str(ind) + "_sigma"])
            fieldsmag[ind].append(field[nn])

    # Code to see quality of fits
    plt.figure()
    plt.plot(x, y)
    plt.plot(x, outy.best_fit)
    plt.show()

    nn += 1

    ##____   Code to review peaks  ______
    """
    plt.figure()
    plt.plot(-y)
    plt.plot(peakind[0],-y[peakind[0]], 'ro')
    """

## Kittel fit check G and Ms

for u in range(len(peakmag)):
    peak = peakmag[u]
    fwhm = fwhmmag[u]
    fields = fieldsmag[u]
    dH = []
    dHa = []

    print("Fit peak number:" + str(u + 1))

    def KittelIP(field, G, Hu):
        return G * ((field) * (field - Hu)) ** 0.5

    def KittelOP(field, G, Hu):
        return (G) * (field + Hu)

    if Anisotropy == "IP":
        Kit = Model(KittelIP)
    else:
        Kit = Model(KittelOP)

    ntira = 0  # numero de pontos a tirar

    newfields = fields[ntira:]
    newpeak1 = peak[ntira:]

    paramK = Kit.make_params()

    paramK.add("Hu", value=Hu)
    paramK.add("G", value=G, min=2.8e6, max=3.1e6)

    Kittelfit = Kit.fit(newpeak1, paramK, field=newfields)
    print(Kittelfit.fit_report())

    G1 = Kittelfit.best_values["G"]

    ## Valor de G a usar no resto dos fits G = 3E6  --  G1 = resultado do fit de kittel

    G_a_usar = G1

    ## dF to dH conversion

    for v in range(len(fwhm)):
        dh = fwhm[v] / (
            G_a_usar * np.sqrt(1 + ((G_a_usar * 4 * np.pi * Ms) / (2 * peak[v])) ** 2)
        )
        dH.append(dh)

    ## damping

    def damping(x, dH0, G, alpha):
        return dH0 + ((2 * x) / G) * alpha

    damp = Model(damping)

    ntirar = 3  # numero de pontos a tirar

    newdH = dH[ntirar:]
    newpeak = peak[ntirar:]
    # print(newpeak)
    parameters2 = damp.make_params(dH0=dH0, G=G_a_usar, alpha=alpha)

    # parameters2['dH0'].vary=False
    parameters2["G"].vary = False

    myfit = damp.fit(newdH, parameters2, x=newpeak)
    print(myfit.fit_report())

    ## app Damping

    for j in range(len(fwhm)):
        dha = fwhm[j] / (G * (2 * (abs(fields[j]) + Ha) + Ms * 4 * math.pi))
        dHa.append(dha)

    minapp = min([t for t in dHa if t > 0])
    num = 0
    soma = 0

    for t in range(len(dH)):
        if fields[t] > 800 and dHa[t] < 1.2 * minapp:
            soma += dHa[t]
            num += 1
            # print(soma, num)

    media = soma / num
    print("App damping=" + str(media))
    baseline = np.full((1, len(dH)), media)

    ## Plotting the result

    plt.figure("Fit_Peak" + str(u + 1))

    plt.subplot(2, 2, 1)
    plt.plot(fields, peak, "ro")
    plt.plot(newfields, Kittelfit.best_fit)
    plt.xlabel("Magnetic Field (Oe)")
    plt.ylabel("Peak (Hz)")

    plt.subplot(2, 2, 2)
    plt.plot(fields, fwhm)
    plt.xlabel("Magnetic Field (Oe)")
    plt.ylabel("FWHM (Hz)")

    plt.subplot(2, 2, 3)
    plt.plot(peak, dH, "ro")
    plt.plot(newpeak, myfit.best_fit)
    plt.ylabel("dH")
    plt.xlabel("Frequency (GHz)")
    valor = "%.6f" % float(myfit.best_values["alpha"])
    plt.text(
        peak[0],
        dH[-2],
        "Gilbert Damping =" + str(valor),
        bbox=dict(facecolor="red", alpha=0.5),
    )

    plt.subplot(2, 2, 4)
    plt.plot(fields, dHa, "ro")
    plt.plot(fields, baseline[0])
    plt.ylabel("app Damp")
    plt.xlabel("Magnetic Field (Oe)")
    valorapp = "%.6f" % float(media)
    plt.text(
        fields[-5],
        dHa[0],
        "Gilbert Damping =" + str(valorapp),
        bbox=dict(facecolor="red", alpha=0.5),
    )

    plt.savefig("Fit_Peak" + str(u + 1) + ".png")

    plt.show()

    ## Making a csv to save the results

    le = len(fields)
    file = open("fields_fwhm_peak_" + str(u + 1) + ".csv", "w")

    file.write("Field,Peak,FWHM" + "\n")
    for i in range(le):
        line = str(fields[i]) + "," + str(peak[i]) + "," + str(fwhm[i])
        file.write(line + "\n")

    file.close()

    serie = myfit.best_fit.tolist()
    for ss in range(ntirar):
        newpeak.insert(0, 0)
        serie.insert(0, 0)

    file = open("Dampings" + str(u + 1) + ".csv", "w")
    file.write(
        "x1_Peak, y1_dH, x2_newpeak, y2_dH_Fit, x3_Field, y3_dHa, y3_baseline" + "\n"
    )
    for k in range(le):
        line = (
            str(peak[k])
            + ","
            + str(dH[k])
            + ","
            + str(newpeak[k])
            + ","
            + str(serie[k])
            + ","
            + str(fields[k])
            + ","
            + str(dHa[k])
            + ","
            + str(baseline[0][k])
        )
        file.write(line + "\n")
    file.close()

    file = open("Fit_Values_Peak" + str(u + 1) + ".csv", "w")
    file.write(Kittelfit.fit_report())
    file.write("\n")
    file.write(myfit.fit_report())
    file.close()
