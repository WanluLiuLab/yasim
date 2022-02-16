# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: logger.py -- System-wide logger.
#
#  VERSION HISTORY:
#  2021-08-10 0.1  : Purposed and added by YU Zhejian, uses logger class.
#  2021-08-11 0.1  : Rewritten using global logger.
#  2021-08-13 0.1  : RegisteredLoggerHandler and chronolog (time_recorder) added.
#  2021-08-15 0.1  : Logger class and RegisteredLoggerHandler deprecated.
#  2021-08-15 0.1  : Trace level added.
#  2021-09-07 0.1  : Logger injection removed.
#
# ==============================================================================

"""
logger.py -- System-wide logger.

Features
--------

This logger defines a "trace" level (TRACE = 8) and a logging decorator.
It can also reads environment variable named 'LOG_LEVEL' and fallback to DEBUG (10) by default.

Usage
-----

See :py:func:`chronolog` and :py:func:`trace`.

There are also :py:func:`set_level` and :py:func:`get_logger`,
which is only snake-case wrappers for those contents inside :py:mod:`logging` standard module.
"""

import logging
import os

__all__ = ['__version__', 'TRACE', 'chronolog', 'set_level', 'get_logger']
__version__ = 0.1

from typing import Union

_SB = os.environ.get('SPHINX_BUILD')

# The following contents adds a new log level called trace.
TRACE = 8


def trace(self, msg, *args, **kwargs):
    """
    Log 'msg % args' with severity 'TRACE'.
    """
    if self.isEnabledFor(TRACE):
        self._log(TRACE, msg, args, **kwargs)


logging.addLevelName(TRACE, "TRACE")
logging.Logger.trace = trace
logging.trace = trace

_lh = logging.getLogger()


def _get_formatter(level: int) -> logging.Formatter:
    if level > logging.DEBUG:
        log_format = '%(asctime)s\t[%(levelname)s] %(message)s'
    else:
        log_format = '%(asctime)s %(name)s:%(lineno)d::%(funcName)s\t[%(levelname)s]\t%(message)s'
    return logging.Formatter(log_format)


def set_level(level: Union[str, int], quiet: bool = True) -> int:
    """
    Set the global logging level, and update the format.
    The log will be more verbose if the level is below debug.

    # FIXME: File set DEBUG.
    """
    this_level = _lh.getEffectiveLevel()
    try:
        _lh.setLevel(level)
        this_level = _lh.getEffectiveLevel()
    except ValueError:
        pass
    if not quiet:
        _lh.info(f'Resetting log level: {logging.getLevelName(this_level)}')

    logging.basicConfig(
        # level=this_level,
        handlers=[
            logging.StreamHandler()
        ]
    )
    file_handler = logging.FileHandler(filename="log.log")
    file_handler.setLevel(logging.DEBUG)
    # logging.root.setLevel(logging.DEBUG)
    for handler in logging.root.handlers:
        handler.formatter = _get_formatter(this_level)
    logging.root.addHandler(file_handler)
    return this_level


if not "_global_level" in locals() and not "_global_level" in globals():
    # set the global log-level.
    # Will read from LOG_LEVEL environment variable.
    # And fall to DEBUG if fails.
    _global_level = os.environ.get('LOG_LEVEL')
    if _global_level is None:
        _global_level = logging.INFO
    _global_level = set_level(_global_level)


def chronolog(display_time: bool = False, log_error: bool = False):
    """
    The logging decorator, will inject a logger variable named lh to the code.
    From <https://stackoverflow.com/questions/17862185/how-to-inject-variable-into-scope-with-a-decorator>

    .. note::
        The :py:func:`error` (or :py:func:`exception`, :py:func:`critical`, :py:func:`fatal`
        functions DO NOT exit the program! You have to exit the program by yourself!

    .. warning::
        Call this function, do NOT call functions inside this function!

    :param display_time: Whether to display calling time, arguments and return value in [TRACE] log level.
    :param log_error: Whether add error captured
    """

    def msg_decorator(f):
        if _SB == '1':
            return f  # To make Sphinx get the right result.

        def inner_dec(*args, **kwargs):
            """
            Decorator which performs the logging and do the work.

            :param args: Unnamed arguments of the decorated function call.
            :param kwargs: Named arguments of the decorated function call.
            :return: The return value of the decorated function call.
            :raise: The return value of the decorated function call.
            """
            try:
                _ = f.__globals__
            except AttributeError:
                return f(*args, **kwargs)
            lh = logging.getLogger(f.__module__)
            if display_time:
                args_repr = [repr(a) for a in args]
                kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
                signature = ", ".join(args_repr + kwargs_repr)
                # FIXME: Function name incorrect!
                lh.trace(f"{f.__name__}({signature})", )
            res = None
            try:
                res = f(*args, **kwargs)
            except Exception as e:
                if log_error:
                    lh.exception(f"Exception inside func: {e}", stack_info=True, exc_info=True)
                raise e
            finally:
                if display_time:
                    lh.trace(f"{f.__name__} -> {res}", )
            return res

        return inner_dec

    return msg_decorator


def get_logger(name: str):
    """
    A Simple logging.getLogger() wrapper.
    """
    return logging.getLogger(name)
