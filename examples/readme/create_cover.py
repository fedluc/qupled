import os
import shutil
import tempfile

import matplotlib as mpl
import numpy as np
from matplotlib import colormaps as cm
from matplotlib import pyplot as plt

import qupled.qstls as qstls


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
    scheme, error = get_plot_data()
    if scheme is not None:
        svg_files = []
        for darkmode in [False, True]:
            svg_files.append(create_one_svg_file(darkmode, scheme, error))
        return svg_files
    else:
        return []


def get_plot_data():
    scheme: qstls.Solver = None
    error: list[float] = []
    while scheme is None or scheme.results.error > 1e-5:
        scheme = solve_qstls(scheme.run_id if scheme is not None else scheme)
        error.append(scheme.results.error)
    return scheme, error


def create_one_svg_file(darkmode: bool, scheme: qstls.Solver, error: np.array):
    # Get plot settings
    settings = PlotSettings(darkmode)
    # Set style
    plt.rcdefaults()
    plt.style.use(settings.theme)
    # Create figure
    plt.figure(figsize=settings.figure_size)
    # Plot quantities of interest
    plot_density_response(plt, scheme, settings)
    plot_ssf(plt, scheme, settings)
    plot_error(plt, error, settings)
    # Combine plots
    plt.tight_layout()
    # Save figure
    file_name = settings.figure_name
    plt.savefig(file_name)
    plt.show()
    plt.close()
    return file_name


def solve_qstls(guess_run_id: int):
    scheme = qstls.Solver()
    inputs = qstls.Input(coupling=15.0, degeneracy=1.0)
    inputs.mixing = 0.3
    inputs.resolution = 0.1
    inputs.cutoff = 10
    inputs.matsubara = 16
    inputs.threads = 16
    inputs.iterations = 0
    inputs.guess = (
        scheme.get_initial_guess(guess_run_id)
        if guess_run_id is not None
        else inputs.guess
    )
    scheme.compute(inputs)
    return scheme


def plot_density_response(plt: plt, scheme: qstls.Solver, settings: PlotSettings):
    results = scheme.results
    inputs = scheme.inputs
    wvg_squared = results.wvg[:, np.newaxis] ** 2
    denominator = wvg_squared + inputs.coupling * (results.idr - results.adr)
    dr = results.idr * wvg_squared / denominator
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
            results.wvg,
            dr[:, parameters[i]],
            color=color,
            linewidth=settings.width,
            label=label,
        )
    plt.xlim(0, settings.xlim)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Density response", fontsize=settings.labelsz, fontweight="bold")
    plt.legend(fontsize=settings.ticksz, loc="upper right")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def plot_ssf(plt: plt, scheme: qstls.Solver, settings: PlotSettings):
    results = scheme.results
    plt.subplot(2, 2, 4)
    plt.plot(results.wvg, results.ssf, color=settings.color, linewidth=settings.width)
    plt.xlim(0, settings.xlim)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Static structure factor", fontsize=settings.labelsz, fontweight="bold")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def plot_error(plt: plt, error: list[float], settings: PlotSettings):
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


if __name__ == "__main__":
    main()
