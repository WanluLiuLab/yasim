from typing import Callable


def copy_doc(copy_func: Callable) -> Callable:
    """
    The following piece of code is from
    https://stackoverflow.com/questions/68901049/copying-the-docstring-of-function-onto-another-function-by-name
    by Iced Chai at Aug 24, 2021 at 2:56

    This wrapper copies docstring from one function to another.

    Use Example: copy_doc(self.copy_func)(self.func) or used as deco

    >>> class Test:
    ...     def foo(self) -> None:
    ...         \"\"\"Woa\"\"\"
    ...         ...
    ...
    ...     @copy_doc(foo)
    ...     def this(self) -> None:
    ...         pass
    >>> Test.this.__doc__
    'Woa'
    """

    def wrapper(func: Callable) -> Callable:
        func.__doc__ = copy_func.__doc__
        return func

    return wrapper
