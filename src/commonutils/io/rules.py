import io
from abc import abstractmethod
from typing import IO, Union, List, Callable

from commonutils.typing import PathType

_IOType = Union[IO, io.IOBase]

# This goes wrong in Python 3.7
try:
    _RuleType = Callable[[PathType, ...], bool]
except TypeError:
    _RuleType = Callable[..., bool]
try:
    _OpenerType = Callable[[PathType, ...], Union[IO, io.IOBase]]
except TypeError:
    _OpenerType = Callable[..., Union[IO, io.IOBase]]


class FileRuleType:
    _rule: _RuleType
    _opener: _OpenerType
    rule_name: str

    @abstractmethod
    def apply_rule(self, path: PathType, *args, **kwargs) -> bool:
        pass

    @abstractmethod
    def get_opener(self, path: PathType, *args, **kwargs) -> _IOType:
        pass


class BaseFileRule(FileRuleType):

    def apply_rule(self, path: PathType, *args, **kwargs) -> bool:
        return self._rule(path, *args, **kwargs)

    def get_opener(self, path: PathType, *args, **kwargs) -> _IOType:
        return self._opener(path, *args, **kwargs)

    def __init__(self, rule: _RuleType, opener: _OpenerType, rule_name: str = "unnames_rule"):
        self._rule = rule
        self._opener = opener
        self.rule_name = rule_name

    @classmethod
    def from_extension(cls, *extensions: str, opener: _OpenerType, rule_name: str = "unnames_rule"):
        def extension_rule(path: PathType, *args, **kwargs):
            for extension in extensions:
                if path.endswith(extension):
                    return True
            else:
                return False

        new_instance = cls(extension_rule, opener=opener, rule_name=rule_name)
        return new_instance


class FileRuleRing:
    _registered_rules: List[FileRuleType] = []

    @staticmethod
    def register(rule: FileRuleType):
        FileRuleRing._registered_rules.append(rule)

    @staticmethod
    def open(path: str, *args, **kwargs) -> _IOType:
        for rule in FileRuleRing._registered_rules:
            if rule.apply_rule(path, *args, **kwargs):
                return rule.get_opener(path, *args, **kwargs)
        return io.open(path, *args, **kwargs)
