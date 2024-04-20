import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
import imageio
import qupled.classic as qpc
import qupled.quantum as qpq

def main():
    darkmode = False
    nIterations = 33
    svg_files = create_all_svg_files(nIterations, darkmode)
    combine_svg_files(svg_files, darkmode)


def create_all_svg_files(nFiles, darkmode):
    fig, ax = plt.subplots()
    images = []
    error = []
    file_names = []
    for i in range(nFiles):
        file_names.append(create_one_svg_file(i, error, darkmode))
    return file_names

def combine_svg_files(svg_files, darkmode):
    svg_template = """
    <svg width="{}" height="{}"
     xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink">
     {}
    </svg>
    """
    image_template = """
    <image xlink:href="{}" width="100%" height="100%" visibility="hidden">
      <animate attributeName="visibility" values="hidden;visible" begin="{}s"/>
      <animate attributeName="visibility" values="visible;hidden" begin="{}s"/>
    </image>
    """
    image_width = 1600 # in pixels
    image_height = 900 # in pixels
    image_duration = 0.18 # in seconds
    svg_image = ""
    animation_file = "qupled_animation_light.svg"
    if (darkmode):
        animation_file = "qupled_animation_dark.svg"
    for i in range(len(svg_files)):
        begin_visible = i * image_duration
        begin_hidden = begin_visible + image_duration
        svg_image += image_template.format(svg_files[i], begin_visible, begin_hidden)
    with open(animation_file, "w") as f:
        f.write(svg_template.format(image_width, image_height, svg_image))


def create_one_svg_file(i, errorList, darkmode):
    # Solve scheme
    [wvg, adr, idr, ssf, error] = solve_qstls(i)
    # Get plot settings
    settings = PlotSettings(darkmode)
    plt.figure(figsize=settings.figure_size)
    plt.style.use(settings.theme)
    # Plot quantities of interest
    plot_density_response(plt, wvg, adr, idr, settings)
    plot_ssf(plt, wvg, ssf, settings)
    plot_error(plt, i, errorList, error, settings)
    # Combine plots
    plt.tight_layout()
    # Save figure
    od = settings.output_directory
    if not os.path.exists(od):
        os.makedirs(od)
    file_name = os.path.join(od, f"plot{i:03}.svg")
    plt.savefig(file_name)
    plt.close()
    return file_name


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
    qstls.inputs.fixed = "adr_fixed_theta1.000_matsubara16.bin"
    qstls.compute()
    return [qstls.scheme.wvg,
            qstls.scheme.adr,
            qstls.scheme.idr,
            qstls.scheme.ssf,
            qstls.scheme.error]
    

def plot_density_response(plt, wvg, adr, idr, settings):
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
        color = settings.colormap(1.0 - 1.0*i/numParameters)
        plt.plot(wvg, dr[:,parameters[i]], color=color,
                 linewidth=settings.width, label=label)
    plt.xlim(0, 6)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Density response", fontsize=settings.labelsz,
              fontweight="bold")
    plt.legend(fontsize=settings.ticksz, loc="lower right")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)



def plot_ssf(plt, wvg, ssf, settings):
    plt.subplot(2, 2, 4)
    plt.plot(wvg, ssf, color=settings.color, linewidth=settings.width)
    plt.xlim(0, 6)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Static structure factor", fontsize=settings.labelsz,
              fontweight="bold")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)



def plot_error(plt, iteration, errorList, error, settings):
    errorList.append(error)
    horizontalLineColor = mpl.rcParams["text.color"]
    plt.subplot(2, 1, 1)
    plt.plot(range(iteration+1), errorList,
             color=settings.color, linewidth=settings.width)
    plt.scatter(iteration, error, color="red", s=150, alpha=1)
    plt.axhline(y=1.0e-5, color=horizontalLineColor, linestyle="--")
    plt.text(3, 1.5e-5, "Convergence", horizontalalignment="center",
             fontsize=settings.ticksz)
    plt.xlim(0, 33)
    plt.ylim(1.0e-6, 1.1e1)
    plt.yscale("log")
    plt.xlabel("Iteration", fontsize=settings.labelsz)
    plt.title("Residual error", fontsize=settings.labelsz,
              fontweight="bold")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


class PlotSettings():

    def __init__(self,
                 darkmode):
        self.labelsz = 16
        self.ticksz = 14
        self.width = 2.0
        self.theme = "ggplot"
        self.colormap = cm["viridis"].reversed()
        self.output_directory = "ligth_mode"
        if darkmode:
            self.theme = "dark_background"
            self.colormap = cm["plasma"]
            self.output_directory = "dark_mode"
        self.color = self.colormap(1.0)
        self.figure_size = (12,8)


if __name__ == "__main__":
    main()
