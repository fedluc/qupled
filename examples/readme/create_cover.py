import os
import tempfile
import shutil
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
from qupled.util import HDF
from qupled.quantum import Qstls


def main():
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        os.chdir(temp_dir)
        try:
            svg_files = create_svg_files()
            for svg_file in svg_files:
                shutil.copy(svg_file, original_dir)
        finally:
            os.chdir(original_dir)


def create_svg_files():
    plot_data, error = get_plot_data()
    if plot_data is not None:
        svg_files = []
        for darkmode in [True, False]:
            svg_files.append(create_one_svg_file(darkmode, plot_data, error))
        return svg_files
    else:
        return []


def get_plot_data():
    plot_data = None
    error = []
    while plot_data is None or plot_data.error > 1e-5:
        plot_data = solve_qstls()
        error.append(plot_data.error)
    return plot_data, error


def create_one_svg_file(darkmode, plot_data, error):
    # Get plot settings
    settings = PlotSettings(darkmode)
    plt.figure(figsize=settings.figure_size)
    plt.style.use(settings.theme)
    # Clip plot data
    plot_data.clip(settings)
    # Plot quantities of interest
    plot_density_response(plt, plot_data, settings)
    plot_ssf(plt, plot_data, settings)
    plot_error(plt, error, settings)
    # Combine plots
    plt.tight_layout()
    # Save figure
    file_name = settings.figure_name
    plt.savefig(file_name)
    plt.close()
    return file_name


def solve_qstls():
    qstls = Qstls()
    rs = 15.0
    theta = 1.0
    inputs = Qstls.Input(rs, theta)
    inputs.mixing = 0.3
    inputs.resolution = 0.1
    inputs.cutoff = 10
    inputs.matsubara = 16
    inputs.threads = 16
    inputs.iterations = 0
    guess_file = f"rs{rs:.3f}_theta{theta:.3f}_QSTLS.h5"
    adr_file = f"adr_fixed_theta{theta:.3f}_matsubara{inputs.matsubara}_QSTLS.bin"
    inputs.guess = (
        Qstls.getInitialGuess(guess_file)
        if os.path.exists(guess_file)
        else inputs.guess
    )
    inputs.fixed = adr_file if os.path.exists(adr_file) else inputs.fixed
    qstls.compute(inputs)
    results = HDF().read(qstls.hdf_file_name, ["wvg", "adr", "ssf", "idr", "error"])
    return QStlsData(
        results["wvg"], results["adr"], results["idr"], results["ssf"], results["error"]
    )


def clip_data(wvg, adr, idr, ssf, error, settings):
    mask = np.less_equal(wvg, settings.xlim)
    wvg = wvg[mask]
    adr = adr[mask]
    idr = idr[mask]
    ssf = ssf[mask]


def plot_density_response(plt, plot_data, settings):
    plot_data.idr[plot_data.idr == 0.0] = 1.0
    dr = np.divide(plot_data.adr, plot_data.idr)
    plt.subplot(2, 2, 3)
    parameters = np.array([0, 1, 2, 3, 4])
    numParameters = parameters.size
    for i in np.arange(numParameters):
        if i == 0:
            label = r"$\omega = 0$"
        else:
            label = r"$\omega = {}\pi/\beta\hbar$".format(parameters[i] * 2)
        color = settings.colormap(1.0 - 1.0 * i / numParameters)
        plt.plot(
            plot_data.wvg,
            dr[:, parameters[i]],
            color=color,
            linewidth=settings.width,
            label=label,
        )
    plt.xlim(0, settings.xlim)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Density response", fontsize=settings.labelsz, fontweight="bold")
    plt.legend(fontsize=settings.ticksz, loc="lower right")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def plot_ssf(plt, plot_data, settings):
    plt.subplot(2, 2, 4)
    plt.plot(
        plot_data.wvg, plot_data.ssf, color=settings.color, linewidth=settings.width
    )
    plt.xlim(0, settings.xlim)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Static structure factor", fontsize=settings.labelsz, fontweight="bold")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def plot_error(plt, error, settings):
    iterations = range(len(error))
    horizontalLineColor = mpl.rcParams["text.color"]
    plt.subplot(2, 1, 1)
    plt.plot(iterations, error, color=settings.color, linewidth=settings.width)
    plt.scatter(iterations[-1], error[-1], color="red", s=150, alpha=1)
    plt.axhline(y=1.0e-5, color=horizontalLineColor, linestyle="--")
    plt.text(
        3, 1.5e-5, "Convergence", horizontalalignment="center", fontsize=settings.ticksz
    )
    plt.xlim(0, 33)
    plt.ylim(1.0e-6, 1.1e1)
    plt.yscale("log")
    plt.xlabel("Iteration", fontsize=settings.labelsz)
    plt.title("Residual error", fontsize=settings.labelsz, fontweight="bold")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def optimise_svg(file_name):
    tmp_file_name = "tmp.svg"
    options = parse_args()
    options.enable_viewboxing = True
    options.strip_ids = True
    options.remove_titles = True
    options.remove_descriptions = True
    options.remove_metadata = True
    options.remove_descriptive_elements = True
    options.indent_type = None
    options.strip_comments = True
    options.strip_xml_space_attribute = True
    options.strip_xml_prolog = True
    options.infilename = file_name
    options.outfilename = tmp_file_name
    (infile, outfile) = getInOut(options)
    start(options, infile, outfile)
    os.rename(tmp_file_name, file_name)


class QStlsData:
    def __init__(self, wvg, adr, idr, ssf, error):
        self.wvg = wvg
        self.adr = adr
        self.idr = idr
        self.ssf = ssf
        self.error = error

    def clip(self, settings):
        mask = np.less_equal(self.wvg, settings.xlim)
        self.wvg = self.wvg[mask]
        self.adr = self.adr[mask]
        self.idr = self.idr[mask]
        self.ssf = self.ssf[mask]


class PlotSettings:
    def __init__(self, darkmode):
        self.labelsz = 16
        self.ticksz = 14
        self.width = 2.0
        self.theme = "dark_background" if darkmode else "ggplot"
        self.colormap = cm["plasma"] if darkmode else cm["viridis"].reversed()
        self.xlim = 6
        self.color = self.colormap(1.0)
        self.figure_size = (12, 8)
        self.figure_name = (
            "qupled_animation_dark.svg" if darkmode else "qupled_animation_light.svg"
        )


if __name__ == "__main__":
    main()
