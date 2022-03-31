"""
parallel_helper.py -- Helper for Multiprocessing

This includes a very basic job pool and some helper classes
"""

import gc
import multiprocessing
import os
import subprocess
import threading
import time
from abc import abstractmethod
from typing import Union, List

from commonutils.importer.tqdm_importer import tqdm

_JOB_TYPE = Union[multiprocessing.Process, threading.Thread]
_PROCESS_TYPE = Union[multiprocessing.Process, subprocess.Popen]


class BaseParallelJobExecutor(threading.Thread):
    """
    Base class of all parallel job executors.

    This class is an abstraction of a sized pool.
    """
    _pool_size: Union[int, float]
    """
    How many jobs is allowed to be executed in one time.

    Use ``math.inf`` to set unlimited or ``0`` to auto determine.
    """

    _pool_name: str
    """
    name of pool to be showed on progress bar, etc.
    """

    _refresh_interval: float
    """
    Interval for probing job status.
    """

    _is_terminated: bool
    """
    Whether termination signal was sent to this thread
    """

    _is_appendable: bool
    """
    Whether this queue is appendable.
    """

    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: Union[int, float] = 0,
                 **kwargs):
        super().__init__()
        self._pool_size = pool_size
        if self._pool_size == 0:
            self._pool_size = multiprocessing.cpu_count()
        self._pool_name = pool_name
        self._is_terminated = False
        self._is_appendable = True

    @abstractmethod
    def append(self, mp_instance: _JOB_TYPE):
        pass

    def stop(self):
        """
        Send termination signal.
        This will stop the job queue from adding more jobs.
        """
        self._is_terminated = True
        # The job queue may receive this signal before being started.
        self._is_appendable = False

    def run(self):
        while not self._is_terminated:
            pass


class ParallelJobExecutor(BaseParallelJobExecutor):
    """
    This is a parallel job executor,
    for jobs in a format of :py:class:`multiprocessing.Process` or :py:class:`threading.Thread`.

    This queue is designed for "batch" jobs.
    That is, the user should append all jobs before they start the executor.

    This executor is designed for non-stated jobs.
    That is, the executor will NOT save the state of any job.
    """

    _pending_job_queue: List[_JOB_TYPE]
    """
    Job waiting to be executed
    """

    _running_job_queue: List[_JOB_TYPE]

    _n_jobs: int
    """
    Number of jobs to be executed
    """

    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: Union[int, float] = 0,
                 refresh_interval: float = 0.01):
        super().__init__(pool_name=pool_name, pool_size=pool_size)
        self._pending_job_queue = []
        self._n_jobs = 0
        self._running_job_queue = []
        self._refresh_interval = refresh_interval

    def run(self):
        """
        Run the queue.
        """
        self._is_appendable = False

        def _scan_through_process():
            """
            Scan through all processes and terminate the exited process.
            """
            for process in self._running_job_queue:
                if not process.is_alive():
                    process.join()
                    if isinstance(process, multiprocessing.Process):
                        process.close()
                    self._running_job_queue.remove(process)
                    del process
                    gc.collect()
                    pbar.update(1)

        pbar = tqdm(desc=self._pool_name, total=self._n_jobs)
        while len(self._pending_job_queue) > 0 and not self._is_terminated:
            while len(self._pending_job_queue) > 0 and len(self._running_job_queue) < self._pool_size:
                new_processs = self._pending_job_queue.pop(0)
                self._running_job_queue.append(new_processs)
                new_processs.start()
            _scan_through_process()
            time.sleep(self._refresh_interval)
        while len(self._running_job_queue) > 0 and not self._is_terminated:
            _scan_through_process()
            time.sleep(self._refresh_interval)
        self._is_terminated = True

    def append(self, mp_instance: _JOB_TYPE):
        """
        Commit a new job to the queue
        """
        if self._is_appendable:
            self._pending_job_queue.append(mp_instance)
            self._n_jobs += 1
        else:
            raise ValueError("Job queue not appendable!")


class TimeOutKiller(threading.Thread):
    """
    A timer that kills a process if time out is reached.

    A process can be either represented using :py:class:`multiprocessing.Process` or :py:class`subprocess.Popen`,
    or by its PID.

    After reaching the timeout, the killer will firstly send SIGTERM (15).
    If the process is alive after 3 seconds, it will send SIGKILL (9).

    TODO: Check whether the PIDs are in the same round.
    """

    _pid: int
    """
    Monitored process ID
    """

    timeout: float
    """
    The timeout in seconds, default 30.0
    """

    def __init__(self, process_or_pid: Union[_PROCESS_TYPE, int], timeout: float = 30.0):
        """
        .. warning :: Initialize the object after starting the monitored process!
        """
        super().__init__()
        if isinstance(process_or_pid, int):
            self._pid = process_or_pid
        else:
            self._pid = process_or_pid.pid
        self.timeout = timeout

    def run(self):
        time.sleep(self.timeout)
        try:
            os.kill(self._pid, 15)
        except (ProcessLookupError, PermissionError):
            return
        time.sleep(3)
        try:
            os.kill(self._pid, 9)
        except (ProcessLookupError, PermissionError):
            pass
