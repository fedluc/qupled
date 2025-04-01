from qupled import native


import functools


class MPI:
    """Class to handle the calls to the MPI API"""

    def __init__(self):
        self.qp_mpi = native.MPI()

    def rank(self):
        """Get rank of the process"""
        return self.qp_mpi.rank()

    def is_root(self):
        """Check if the current process is root (rank 0)"""
        return self.qp_mpi.is_root()

    def barrier(self):
        """Setup an MPI barrier"""
        self.qp_mpi.barrier()

    def timer(self):
        """Get wall time"""
        return self.qp_mpi.timer()

    @staticmethod
    def run_only_on_root(func):
        """Python decorator for all methods that have to be run only by root"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if MPI().is_root():
                return func(*args, **kwargs)

        return wrapper

    @staticmethod
    def synchronize_ranks(func):
        """Python decorator for all methods that need rank synchronization"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            MPI().barrier()

        return wrapper

    @staticmethod
    def record_time(func):
        """Python decorator for all methods that have to be timed"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            mpi = MPI()
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
