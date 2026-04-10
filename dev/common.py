from pathlib import Path


def get_wheel_file():
    wheel_file = list(Path().rglob("qupled*.whl"))
    if not wheel_file:
        print("No .whl files found. Ensure the package is built first.")
        return None
    return str(wheel_file[0])
