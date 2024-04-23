import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
from scour.scour import start, getInOut, parse_args
import qupled.classic as qpc
import qupled.quantum as qpq

def main():
    darkmode = True
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
    <svg width="864pt" height="576pt" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
     {}
    </svg>
    """
    image_template = """
    <g visibility="hidden">
     {}
    </g>
    """
    image_duration = 0.18 # in seconds
    svg_image = ""
    animation_file = "qupled_animation_light.svg"
    if (darkmode):
        animation_file = "qupled_animation_dark.svg"
    for i in range(len(svg_files)):
        svg_file = svg_files[i]
        begin_visible = i * image_duration
        begin_hidden = begin_visible + image_duration
        with open(svg_file, "r") as f:
            svg_content = f.read()
            svg_content = svg_content[svg_content.index("<svg"):]
            svg_content = add_animation(svg_content, begin_visible, begin_hidden)
            svg_image += image_template.format(svg_content)
        os.remove(svg_file)
    with open(animation_file, "w") as fw:
         fw.write(svg_template.format(svg_image))
    optimise_svg(animation_file)


def add_animation(svg_content, begin, end):
    animation_xml = f"""
     <animate attributeName="visibility" values="hidden;visible" begin="{begin}s"/>
     <animate attributeName="visibility" values="visible;hidden" begin="{end}s"/>
    """
    index = svg_content.index("</svg>")
    modified_svg_content = svg_content[:index] + animation_xml + svg_content[index:]
    return modified_svg_content

         
def create_one_svg_file(i, errorList, darkmode):
    # Solve scheme
    plot_data = solve_qstls(i)
    # Get plot settings
    settings = PlotSettings(darkmode)
    plt.figure(figsize=settings.figure_size)
    plt.style.use(settings.theme)
    # Clip plot data
    plot_data.clip(settings)
    # Plot quantities of interest
    plot_density_response(plt, plot_data, settings)
    plot_ssf(plt, plot_data, settings)
    plot_error(plt, i, errorList, plot_data.error, settings)
    # Combine plots
    plt.tight_layout()
    # Save figure
    file_name = f"plot{i:03}.svg"
    plt.savefig(file_name)
    plt.close()
    # Optimise svg file
    optimise_svg(file_name)
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
    return QStlsData(qstls.scheme.wvg,
                     qstls.scheme.adr,
                     qstls.scheme.idr,
                     qstls.scheme.ssf,
                     qstls.scheme.error)


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
        if (i == 0) :
            label = r"$\omega = 0$"
        else:
            label = r"$\omega = {}\pi/\beta\hbar$".format(parameters[i]*2)
        color = settings.colormap(1.0 - 1.0*i/numParameters)
        plt.plot(plot_data.wvg, dr[:,parameters[i]], color=color,
                 linewidth=settings.width, label=label)
    plt.xlim(0, settings.xlim)
    plt.xlabel("Wave-vector", fontsize=settings.labelsz)
    plt.title("Density response", fontsize=settings.labelsz,
              fontweight="bold")
    plt.legend(fontsize=settings.ticksz, loc="lower right")
    plt.xticks(fontsize=settings.ticksz)
    plt.yticks(fontsize=settings.ticksz)


def plot_ssf(plt, plot_data, settings):
    plt.subplot(2, 2, 4)
    plt.plot(plot_data.wvg, plot_data.ssf, color=settings.color,
             linewidth=settings.width)
    plt.xlim(0, settings.xlim)
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


class QStlsData():

    def __init__(self, wvg, adr, idr,
                 ssf, error):
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
    
class PlotSettings():

    def __init__(self, darkmode):
        self.labelsz = 16
        self.ticksz = 14
        self.width = 2.0
        self.theme = "ggplot"
        self.colormap = cm["viridis"].reversed()
        self.xlim = 6
        if darkmode:
            self.theme = "dark_background"
            self.colormap = cm["plasma"]
        self.color = self.colormap(1.0)
        self.figure_size = (12,8)


if __name__ == "__main__":
    main()
