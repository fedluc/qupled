import platform

system = platform.system()
if system == "Linux":
    import qupled.Linux.qupled as native
elif system == "Darwin":
    import qupled.Darwin.qupled as native
else:
    raise ImportError(f"Platform {system} is not supported")
