# RA, 2020-12-04

import types
import inspect
import contextlib
import pathlib


def import_from_source(name: str, file_path: str) -> types.ModuleType:
    """
    Usage:
        mod = import_from_source("main", str(path / "main.py"))
        mod.some_function()

    From
    https://stackoverflow.com/questions/33469246
    """
    from importlib.machinery import ModuleSpec, SourceFileLoader
    from importlib.util import spec_from_loader, module_from_spec
    loader: SourceFileLoader = SourceFileLoader(name, file_path)
    spec: ModuleSpec = spec_from_loader(loader.name, loader)
    module: types.ModuleType = module_from_spec(spec)
    loader.exec_module(module)
    return module


def whatsmyname() -> str:
    return inspect.currentframe().f_back.f_code.co_name


@contextlib.contextmanager
def redirect(stdout: pathlib.Path = None, stderr: pathlib.Path = None):
    from contextlib import redirect_stdout, redirect_stderr

    with contextlib.ExitStack() as stack:
        if stdout is not None:
            stdout = stack.enter_context(stdout.open(mode="w"))
            stack.enter_context(redirect_stdout(stdout))
        if stderr is not None:
            stderr = stack.enter_context(stderr.open(mode="w"))
            stack.enter_context(redirect_stderr(stderr))
        yield
