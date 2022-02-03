"""Helper for multiprocessing"""
import gc
import multiprocessing
import threading
import time
from typing import Union

from commonutils.tqdm_importer import tqdm


class ParallelJobQueue(threading.Thread):
    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: int = 0,
                 refresh_interval: float = 0.01):
        super().__init__()
        self.pool_size = pool_size
        if self.pool_size == 0:
            self.pool_size = multiprocessing.cpu_count()
        self.pending_job_queue = []
        self._is_terminated = False
        self._max_queue_len = 0
        self.pool_name = pool_name
        self.running_job_queue = []
        self.refresh_interval = refresh_interval

    def run(self):
        def _scan_through_process():
            """
            Scan through all processes and terminate the exited process.
            """
            for process in self.running_job_queue:
                if not process.is_alive():
                    process.join()
                    if isinstance(process, multiprocessing.Process):
                        process.close()
                    self.running_job_queue.remove(process)
                    del process
                    gc.collect()
                    pbar.update(1)

        pbar = tqdm(desc=self.pool_name, total=self._max_queue_len)
        while len(self.pending_job_queue) > 0 and not self._is_terminated:
            while len(self.pending_job_queue) > 0 and len(self.running_job_queue) < self.pool_size:
                new_processs = self.pending_job_queue[0]
                self.pending_job_queue.remove(new_processs)
                self.running_job_queue.append(new_processs)
                new_processs.start()
            _scan_through_process()
            time.sleep(self.refresh_interval)
        while len(self.running_job_queue) > 0 and not self._is_terminated:
            _scan_through_process()
            time.sleep(self.refresh_interval)
        pbar.close()
        self._is_terminated = True

    def stop(self):
        self._is_terminated = True

    def append(self, mp_instance: Union[multiprocessing.Process, threading.Thread]):
        self.pending_job_queue.append(mp_instance)
        self._max_queue_len += 1
