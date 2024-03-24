import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
import imageio
import qupled.classic as qpc
import qupled.quantum as qpq

darkmode = False
if darkmode:
    theme = "dark_background"
    colormap = cm["plasma"]
    animationFile = "qupled_animation_dark.gif"
else :
    theme = "ggplot"
    colormap = cm["viridis"].reversed()
    animationFile = "qupled_animation_light.gif"
nIterations = 33
color = colormap(1.0)    
labelsz = 16
ticksz = 14
width = 2.0

def solve_qstls(i):
    qstls = qpq.Qstls(15.0, 1.0,
                      mixing = 0.3,
                      resolution = 0.1,
                      cutoff = 10,
                      matsubara = 16,
                      threads = 16,
                      iterations = 0)
    if (i > 0):
        qstls.setGuess("rs15.000_theta1.000_QSTLS.h5")
    qstls.inputs.fixed = "adr_fixed_rs15.000_theta1.000_QSTLS.bin"
    qstls.compute()
    return [qstls.scheme.wvg,
            qstls.scheme.adr,
            qstls.scheme.idr,
            qstls.scheme.ssf,
            qstls.scheme.error]

def plot_ssf(plt, wvg, ssf):
    plt.subplot(2, 2, 4)
    plt.plot(wvg, ssf, color=color, linewidth=width)
    plt.xlim(0, 6)
    # plt.ylim(0, 1.4)
    plt.xlabel("Wave-vector", fontsize=labelsz)
    plt.title("Static structure factor", fontsize=labelsz, fontweight="bold")
    plt.xticks(fontsize=ticksz)
    plt.yticks(fontsize=ticksz)

def plot_density_response(plt, wvg, adr, idr):
    idr[idr == 0.0] = 1.0
    dr = np.divide(adr, idr)
    plt.subplot(2, 2, 3)
    parameters = np.array([0, 1, 2, 3, 4])
    numParameters = parameters.size
    for i in np.arange(numParameters):
        if (i == 0) :
            label = r"$\omega = 0$"
        else:
            label = r"$\omega = {}\pi/\beta\hbar$".format(parameters[i]*2)
        color = colormap(1.0 - 1.0*i/numParameters)
        plt.plot(wvg, dr[:,parameters[i]], color=color, linewidth=width, label=label)
    plt.xlim(0, 6)
    # plt.ylim(0, 1.4)
    plt.xlabel("Wave-vector", fontsize=labelsz)
    plt.title("Density response", fontsize=labelsz, fontweight="bold")
    plt.legend(fontsize=ticksz, loc="lower right")
    plt.xticks(fontsize=ticksz)
    plt.yticks(fontsize=ticksz)

def plot_error(plt, iteration, errorList, error):
    errorList.append(error)
    horizontalLineColor = mpl.rcParams["text.color"]
    plt.subplot(2, 1, 1)
    plt.plot(range(iteration+1), errorList, color=color, linewidth=width)
    plt.scatter(iteration, error, color="red", s=150, alpha=1)
    plt.axhline(y=1.0e-5, color=horizontalLineColor, linestyle="--")
    plt.text(3, 1.5e-5, "Convergence", horizontalalignment="center",
             fontsize=ticksz)
    plt.xlim(0, nIterations)
    plt.ylim(1.0e-6, 1.1e1)
    plt.yscale("log")
    plt.xlabel("Iteration", fontsize=labelsz)
    plt.title("Residual error", fontsize=labelsz, fontweight="bold")
    plt.xticks(fontsize=ticksz)
    plt.yticks(fontsize=ticksz)
    
def create_plot(i, errorList):
    [wvg, adr, idr, ssf, error] = solve_qstls(i)
    plt.figure(figsize=(12, 8))
    plt.style.use(theme)
    plot_density_response(plt, wvg, adr, idr)
    plot_ssf(plt, wvg, ssf)
    plot_error(plt, i, errorList, error)
    plt.tight_layout()
    plotName = "plot" + str(i) + ".png"
    plt.savefig(plotName, dpi=300)
    plt.close()
    return plotName

fig, ax = plt.subplots()
images = []
error = []
for i in range(nIterations):
    plotName = create_plot(i, error)
    images.append(imageio.v2.imread(plotName)) 
    os.remove(plotName)
imageio.mimsave(animationFile, images, fps=4)
