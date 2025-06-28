import time


def timer(func):
    """
    A decorator that measures and prints the execution time of the decorated function.

    The elapsed time is displayed in hours, minutes, and seconds as appropriate.

    Args:
        func (callable): The function whose execution time is to be measured.

    Returns:
        callable: A wrapper function that executes the original function and prints the elapsed time.
    """

    def wrapper(*args, **kwargs):
        tic = time.perf_counter()
        result = func(*args, **kwargs)
        toc = time.perf_counter()
        dt = toc - tic
        hours = int(dt // 3600)
        minutes = int((dt % 3600) // 60)
        seconds = dt % 60
        if hours > 0:
            print(f"Elapsed time: {hours} h, {minutes} m, {seconds:.1f} s.")
        elif minutes > 0:
            print(f"Elapsed time: {minutes} m, {seconds:.1f} s.")
        else:
            print(f"Elapsed time: {seconds:.1f} s.")
        return result

    return wrapper
