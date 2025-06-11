import functools

from qupled import native


class MPI:
    """Class to handle the calls to the MPI API"""

    def rank(self):
        """Get rank of the process"""
        return native.MPI.rank()

    def is_root(self):
        """Check if the current process is root (rank 0)"""
        return native.MPI.is_root()

    def barrier(self):
        """Setup an MPI barrier"""
        native.MPI.barrier()

    def timer(self):
        """Get wall time"""
        return native.MPI.timer()

    @staticmethod
    def run_only_on_root(func):
        """Python decorator for all methods that have to be run only by root"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if native.MPI.is_root():
                return func(*args, **kwargs)

        return wrapper

    @staticmethod
    def synchronize_ranks(func):
        """Python decorator for all methods that need rank synchronization"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            native.MPI.barrier()

        return wrapper

    @staticmethod
    def record_time(func):
        """Python decorator for all methods that have to be timed"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            mpi = native.MPI
            tic = mpi.timer()
            func(*args, **kwargs)
            toc = mpi.timer()
            dt = toc - tic
            hours = dt // 3600
            minutes = (dt % 3600) // 60
            seconds = dt % 60
            if mpi.is_root():
                if hours > 0:
                    print("Elapsed time: %d h, %d m, %d s." % (hours, minutes, seconds))
                elif minutes > 0:
                    print("Elapsed time: %d m, %d s." % (minutes, seconds))
                else:
                    print("Elapsed time: %.1f s." % seconds)

        return wrapper
