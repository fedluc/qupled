import os
import pytest
import numpy as np
from qupled.util import Plot


@pytest.fixture
def plot_instance():
    return Plot


def test_plot1D(plot_instance, mocker):
    mockPlotShow = mocker.patch("matplotlib.pyplot.show")
    x = np.arange(0, 10, 0.1)
    y = np.zeros(len(x))
    plot_instance.plot1D(x, y, "x", "y")
    assert mockPlotShow.call_count == 1


def test_plot1DParametric(plot_instance, mocker):
    mockPlotShow = mocker.patch("matplotlib.pyplot.show")
    x = np.arange(0, 10, 0.1)
    y = np.zeros((len(x), 2))
    parameter = np.array([0, 1])
    plot_instance.plot1DParametric(x, y, "x", "y", parameter)
    assert mockPlotShow.call_count == 1
    try:
        parameter = np.array([1, 2])
        plot_instance.plot1DParametric(x, y, "x", "y", parameter)
        assert False
    except:
        assert mockPlotShow.call_count == 1
