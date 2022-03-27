import copy
from io import StringIO
from collections.abc import Iterable
import abc
from typing import Union, Dict, Any, Type, TypeVar

import numpy as np
import yaml

from ._typehints import FileHandle
from . import Rotation
from . import util

MyType = TypeVar('MyType', bound='Config')

class NiceDumper(yaml.SafeDumper):
    """Make YAML readable for humans."""

    def write_line_break(self,
                         data: str = None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()

    def increase_indent(self,
                        flow: bool = False,
                        indentless: bool = False):
        return super().increase_indent(flow, False)

    def represent_data(self,
                       data: Any):
        """Cast Config objects and its subclasses to dict."""
        if isinstance(data, dict) and type(data) != dict:
            return self.represent_data(dict(data))
        if isinstance(data, np.ndarray):
            return self.represent_data(data.tolist())
        if isinstance(data, Rotation):
            return self.represent_data(data.quaternion.tolist())
        else:
            return super().represent_data(data)

    def ignore_aliases(self,
                       data: Any) -> bool:
        """Do not use references to existing objects."""
        return True

class Config(dict):
    """YAML-based configuration."""

    def __init__(self,
                 yml: Union[str, Dict[str, Any]] = None,
                 **kwargs):
        """Initialize from YAML, dict, or key=value pairs."""
        if isinstance(yml,str):
            kwargs.update(yaml.safe_load(yml))
        elif isinstance(yml,dict):
            kwargs.update(yml)

        super().__init__(**kwargs)

    def __repr__(self) -> str:
        """Show as in file."""
        output = StringIO()
        self.save(output)
        output.seek(0)
        return ''.join(output.readlines())


    def __copy__(self: MyType) -> MyType:
        """Create deep copy."""
        return copy.deepcopy(self)

    copy = __copy__


    def __or__(self: MyType,
               other) -> MyType:
        """
        Update configuration with contents of other.

        Parameters
        ----------
        other : damask.Config or dict
            Key-value pairs that update self.

        Returns
        -------
        updated : damask.Config
            Updated configuration.

        Note
        ----
        This functionality is a backport for Python 3.8

        """
        duplicate = self.copy()
        duplicate.update(other)
        return duplicate


    def __ior__(self: MyType,
                other) -> MyType:
        """Update configuration with contents of other."""
        return self.__or__(other)


    def delete(self: MyType,
               keys: Union[Iterable, str]) -> MyType:
        """
        Remove configuration keys.

        Parameters
        ----------
        keys : iterable or scalar
            Label of the key(s) to remove.

        Returns
        -------
        updated : damask.Config
            Updated configuration.

        """
        duplicate = self.copy()
        for k in keys if isinstance(keys, Iterable) and not isinstance(keys, str) else [keys]:
            del duplicate[k]
        return duplicate


    @classmethod
    def load(cls: Type[MyType],
             fname: FileHandle) -> MyType:
        """
        Load from yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.

        Returns
        -------
        loaded : damask.Config
            Configuration from file.

        """
        return cls(yaml.safe_load(util.open_text(fname)))

    def save(self,
             fname: FileHandle,
             **kwargs):
        """
        Save to yaml file.

        Parameters
        ----------
        fname : file, str, or pathlib.Path
            Filename or file for writing.
        **kwargs : dict
            Keyword arguments parsed to yaml.dump.

        """
        if 'width' not in kwargs:
            kwargs['width'] = 256
        if 'default_flow_style' not in kwargs:
            kwargs['default_flow_style'] = None
        if 'sort_keys' not in kwargs:
            kwargs['sort_keys'] = False

        fhandle = util.open_text(fname,'w')
        try:
            fhandle.write(yaml.dump(self,Dumper=NiceDumper,**kwargs))
        except TypeError:                                                                           # compatibility with old pyyaml
            del kwargs['sort_keys']
            fhandle.write(yaml.dump(self,Dumper=NiceDumper,**kwargs))


    @property
    @abc.abstractmethod
    def is_complete(self):
        """Check for completeness."""
        raise NotImplementedError


    @property
    @abc.abstractmethod
    def is_valid(self):
        """Check for valid file layout."""
        raise NotImplementedError
