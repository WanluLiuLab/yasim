# TODO: to be refactored

from typing import Type, Callable, TypeVar, Any

_HookedType = TypeVar("_HookedType")

try:
    _HookType = Callable[[_HookedType, ...], Any]
except TypeError:
    _HookType = Callable[..., Any]
_NewAttrType = TypeVar("_NewAttrType")


def register_new_attribute(
        cls: Type[_HookedType],
        attribute_name: str,
        attribute_type: Type[_NewAttrType] = Any
):
    setattr(cls, attribute_name, None)


def register_new_hook(cls: Type[_HookedType], hook: _HookType):
    setattr(cls, hook.__name__, hook)


def register_new_attribute_decorator(cls: Type[_HookedType]):
    def wrapper(attribute_name: str, attribute_type: Type[_NewAttrType] = Any):
        register_new_attribute(cls, attribute_name, attribute_type)

    return wrapper


def register_new_hook_decorator(cls: Type[_HookedType]):
    def wrapper(hook: _HookType):
        register_new_hook(cls, hook)

    return wrapper


def hookable_decorator(cls: Type[_HookedType]) -> Type[_HookedType]:
    cls.register_new_attribute = register_new_attribute_decorator(cls)
    cls.register_new_hook = register_new_hook_decorator(cls)
    return cls


if __name__ == "__main__":
    @hookable_decorator
    class A: pass


    def a_func(_a: A): return id(_a)


    A.register_new_hook(a_func)
    A.register_new_attribute("attr", int)
    a = A()
    print(a.a_func())
    assert a.attr is None
