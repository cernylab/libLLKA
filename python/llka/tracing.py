"""Debug tracing of the libLLKA internals (``llka_tracing.h``)."""

from __future__ import annotations

import ctypes
from ctypes import byref
from typing import NamedTuple

from . import _capi
from ._capi import lib

__all__ = ['TracepointInfo', 'toggle_all_tracepoints', 'toggle_tracepoint', 'trace', 'tracepoint_info', 'tracepoint_state']


class TracepointInfo(NamedTuple):
    """A tracepoint and the description of what it logs."""

    tpid: int
    description: str


def toggle_all_tracepoints(enabled: bool) -> None:
    """Enables or disables all tracepoints."""
    lib.LLKA_toggleAllTracepoints(bool(enabled))


def toggle_tracepoint(tpid: int, enabled: bool) -> None:
    """Enables or disables one tracepoint."""
    lib.LLKA_toggleTracepoint(tpid, bool(enabled))


def trace(dont_clear: bool = False) -> str:
    """Returns the trace logged so far. The trace log is cleared unless ``dont_clear`` is true.

    Returns an empty string if libLLKA was built without tracing.
    """
    ptr = lib.LLKA_trace(bool(dont_clear))
    if not ptr:
        return ''
    try:
        return ctypes.string_at(ptr).decode('utf-8', errors='replace')
    finally:
        lib.LLKA_destroyTrace(ptr)


def tracepoint_info() -> list[TracepointInfo]:
    """Returns all available tracepoints."""
    infos = lib.LLKA_tracepointInfo()
    try:
        return [
            TracepointInfo(infos.infos[idx].TPID, _capi.decode(infos.infos[idx].description))
            for idx in range(infos.nInfos)
        ]
    finally:
        if infos.infos:
            lib.LLKA_destroyTracepointInfo(byref(infos))


def tracepoint_state(tpid: int) -> bool:
    """Returns whether a tracepoint is enabled."""
    return bool(lib.LLKA_tracepointState(tpid))
