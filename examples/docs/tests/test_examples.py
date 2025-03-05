import os
import sys
import glob
import importlib
import pytest


def clean_example_files():
    for file_extension in ["h5", "bin", "zip"]:
        file_names = glob.glob("*." + file_extension)
        for file_name in file_names:
            os.remove(file_name)


def run_example(example_name, mocker):
    examples_dir = os.path.abspath("examples/docs")
    try:
        if examples_dir not in sys.path:
            sys.path.insert(0, examples_dir)
        mocker.patch("matplotlib.pyplot.show")
        importlib.import_module(example_name)
    finally:
        clean_example_files()
        if examples_dir in sys.path:
            sys.path.remove(examples_dir)


def run_example_with_error(example_name, mocker, expected_error_message):
    try:
        mocker.patch("matplotlib.pyplot.show")
        with pytest.raises(SystemExit) as excinfo:
            importlib.import_module(example_name)
        assert excinfo.value.code == expected_error_message
    finally:
        clean_example_files()


def test_fixed_adr_qstls(mocker):
    run_example_with_error(
        "fixedAdrQstls", mocker, "Error while solving the dielectric theory"
    )


def test_fixed_adr_qstls_iet(mocker):
    run_example("fixedAdrQstlsIet", mocker)


def test_fixed_adr_qvs_stls(mocker):
    run_example("fixedAdrQVSStls", mocker)


def test_initial_guess_qstls(mocker):
    run_example("initialGuessQstls", mocker)


def test_initial_guess_qstls_iet(mocker):
    run_example("initialGuessQstlsIet", mocker)


def test_initial_guess_stls(mocker):
    run_example("initialGuessStls", mocker)


def test_solve_quantum_schemes(mocker):
    run_example("solveQuantumSchemes", mocker)


def test_solve_qvs_stls(mocker):
    run_example("solveQVSStls", mocker)


def test_solve_rpa_and_esa(mocker):
    run_example("solveRpaAndESA", mocker)


def test_solve_stls(mocker):
    run_example("solveStls", mocker)


def test_solve_stls_iet(mocker):
    run_example("solveStlsIet", mocker)


def test_solve_vs_stls(mocker):
    run_example("solveVSStls", mocker)
