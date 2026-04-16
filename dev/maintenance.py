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
