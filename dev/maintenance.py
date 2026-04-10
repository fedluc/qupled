import shutil
from pathlib import Path


def clean():
    folders_to_clean = [
        Path("dist"),
        Path("dist-native"),
        Path("dist-native-tests"),
        Path("dist-native-only"),
        Path("src", "qupled.egg-info"),
        Path("docs", "_build"),
        Path("docs", "_generated"),
    ]
    for folder in folders_to_clean:
        if folder.exists():
            print(f"Removing folder: {folder}")
            shutil.rmtree(folder)


def update_version(build_version):
    pyproject_file = Path("pyproject.toml")
    if not pyproject_file.exists():
        return
    with pyproject_file.open("r") as file:
        content = file.readlines()
    with pyproject_file.open("w") as file:
        for line in content:
            if line.startswith("version = "):
                file.write(f'version = "{build_version}"')
                file.write("\n")
            else:
                file.write(line)
