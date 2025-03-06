import numpy as np
from matplotlib import pyplot as plt
from qupled.util import Plot


def test_plot1D(mocker):
    mock_plot_show = mocker.patch.object(plt, plt.show.__name__)
    x = np.arange(0, 10, 0.1)
    y = np.zeros(len(x))
    Plot.plot_1D(x, y, "x", "y")
    assert mock_plot_show.call_count == 1


def test_plot1D_parametric(mocker):
    mock_plot_show = mocker.patch.object(plt, plt.show.__name__)
    x = np.arange(0, 10, 0.1)
    y = np.zeros((len(x), 2))
    parameter = np.array([0, 1])
    Plot.plot_1D_parametric(x, y, "x", "y", parameter)
    assert mock_plot_show.call_count == 1
    try:
        parameter = np.array([1, 2])
        Plot.plot_1D_parametric(x, y, "x", "y", parameter)
        assert False
    except:
        assert mock_plot_show.call_count == 1
