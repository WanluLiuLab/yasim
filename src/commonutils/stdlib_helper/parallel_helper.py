"""
parallel_helper.py -- Helper for Multiprocessing

This includes a very basic job pool.
"""

import gc
import multiprocessing
import os
import subprocess
import threading
import time
from typing import Union, List

from commonutils.importer.tqdm_importer import tqdm

_JOB_TYPE = Union[multiprocessing.Process, threading.Thread]
_PROCESS_TYPE = Union[multiprocessing.Process, subprocess.Popen]


class ParallelJobQueue(threading.Thread):
    """
    This is a parallel job queue,
    for jobs in a format of :py:class:`multiprocessing.Process` or :py:class:`threading.Thread`.

    This queue is designed for "batch" jobs.
    That is, the user should append all jobs before they start the queue.


    """

    pool_size: Union[int, float]
    """
    How many jobs is allowed to be executed in one time.
    
    Use ``math.inf`` to set unlimited or ``0`` to auto determine.
    """

    pool_name: str
    """
    name of pool to be showed on progress bar, etc.
    """

    refresh_interval: float
    """
    Interval for probing job status.
    """

    _pending_job_queue: List[_JOB_TYPE]
    """
    Job waiting to be executed
    """

    _is_terminated: bool
    """
    Whether termination signal was sent to this thread
    """

    _is_appendable: bool
    """
    Whether this queue is appendable.
    """

    _max_queue_len: int
    """
    Number of jobs to be executed
    """

    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: Union[int, float] = 0,
                 refresh_interval: float = 0.01):
        super().__init__()
        self.pool_size = pool_size
        if self.pool_size == 0:
            self.pool_size = multiprocessing.cpu_count()
        self._pending_job_queue = []
        self._is_terminated = False
        self._max_queue_len = 0
        self.pool_name = pool_name
        self._running_job_queue = []
        self.refresh_interval = refresh_interval
        self._is_appendable = True

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

        pbar = tqdm(desc=self.pool_name, total=self._max_queue_len)
        while len(self._pending_job_queue) > 0 and not self._is_terminated:
            while len(self._pending_job_queue) > 0 and len(self._running_job_queue) < self.pool_size:
                new_processs = self._pending_job_queue.pop(0)
                self._running_job_queue.append(new_processs)
                new_processs.start()
            _scan_through_process()
            time.sleep(self.refresh_interval)
        while len(self._running_job_queue) > 0 and not self._is_terminated:
            _scan_through_process()
            time.sleep(self.refresh_interval)
        self._is_terminated = True

    def stop(self):
        """
        Send termination signal.
        This will stop the job queue from adding more jobs.
        """
        self._is_terminated = True
        # The job queue may receive this signal before being started.
        self._is_appendable = False

    def append(self, mp_instance: _JOB_TYPE):
        """
        Commit a new job to the queue
        """
        if self._is_appendable:
            self._pending_job_queue.append(mp_instance)
            self._max_queue_len += 1
        else:
            raise ValueError("Job queue not appendable!")


class TimeOutKiller(threading.Thread):
    _pid:int
    timeout: int

    def __init__(self, process_or_pid:Union[_PROCESS_TYPE, int], timeout: int = 30):
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
        except ProcessLookupError:
            return
        time.sleep(3)
        try:
            os.kill(self._pid, 9)
        except ProcessLookupError:
            pass
