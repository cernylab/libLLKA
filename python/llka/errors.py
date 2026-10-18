"""Exceptions raised when a libLLKA function reports an error.

Every exception derives from :class:`LLKAError` and carries the original
:class:`~llka.enums.RetCode` in its ``retcode`` attribute. Where a standard
Python exception describes the problem well, the LLKA exception also derives
from it so that callers may catch e.g. ``ValueError`` or ``FileNotFoundError``.
"""

from __future__ import annotations

from .enums import RetCode

__all__ = [
    'BadAverageNuAnglesError',
    'BadClassificationClustersError',
    'BadClassificationLimitsError',
    'BadConfalsError',
    'BadDataError',
    'BadGoldenStepsError',
    'CannotReadFileError',
    'CifParseError',
    'ClassificationContextError',
    'InvalidArgumentError',
    'LLKAError',
    'LLKANotImplementedError',
    'MismatchingDataError',
    'MismatchingSizesError',
    'MissingAtomsError',
    'MultipleAltIdsError',
    'NoDataError',
    'NoFileError',
    'NothingToClassifyError',
]


class LLKAError(Exception):
    """Base class of all errors reported by libLLKA."""

    retcode: RetCode = RetCode.OK

    def __init__(self, message: str, retcode: RetCode | None = None) -> None:
        super().__init__(message)
        if retcode is not None:
            self.retcode = retcode


class InvalidArgumentError(LLKAError, ValueError):
    """An invalid argument was passed to a function."""

    retcode = RetCode.E_INVALID_ARGUMENT


class MismatchingSizesError(LLKAError, ValueError):
    """Arguments that were supposed to have the same size did not."""

    retcode = RetCode.E_MISMATCHING_SIZES


class LLKANotImplementedError(LLKAError, NotImplementedError):
    """The function is not implemented in this build of libLLKA."""

    retcode = RetCode.E_NOT_IMPLEMENTED


class MultipleAltIdsError(LLKAError, ValueError):
    """A structure contains atoms from more than one alternate position."""

    retcode = RetCode.E_MULTIPLE_ALT_IDS


class MissingAtomsError(LLKAError, ValueError):
    """A structure does not contain atoms that are required."""

    retcode = RetCode.E_MISSING_ATOMS


class MismatchingDataError(LLKAError, ValueError):
    """Data that were supposed to match did not."""

    retcode = RetCode.E_MISMATCHING_DATA


class ClassificationContextError(LLKAError, ValueError):
    """Base class of errors raised when a classification context cannot be initialized."""


class BadClassificationClustersError(ClassificationContextError):
    """Classification cluster definitions are invalid."""

    retcode = RetCode.E_BAD_CLASSIFICATION_CLUSTERS


class BadGoldenStepsError(ClassificationContextError):
    """Golden step definitions are invalid."""

    retcode = RetCode.E_BAD_GOLDEN_STEPS


class BadConfalsError(ClassificationContextError):
    """Confal definitions are invalid."""

    retcode = RetCode.E_BAD_CONFALS


class BadAverageNuAnglesError(ClassificationContextError):
    """Average nu angles definitions are invalid."""

    retcode = RetCode.E_BAD_AVERAGE_NU_ANGLES


class BadClassificationLimitsError(ClassificationContextError):
    """Classification limits are invalid."""

    retcode = RetCode.E_BAD_CLASSIFICATION_LIMITS


class NoFileError(LLKAError, FileNotFoundError):
    """A file does not exist."""

    retcode = RetCode.E_NO_FILE


class CannotReadFileError(LLKAError, OSError):
    """A file exists but cannot be read."""

    retcode = RetCode.E_CANNOT_READ_FILE


class BadDataError(LLKAError, ValueError):
    """Data contains wrong or unexpected values."""

    retcode = RetCode.E_BAD_DATA


class CifParseError(BadDataError):
    """Cif data could not be parsed or processed. The message describes the problem."""


class NoDataError(LLKAError, ValueError):
    """A data block is empty when it should not be."""

    retcode = RetCode.E_NO_DATA


class NothingToClassifyError(LLKAError, ValueError):
    """There was nothing to classify, e.g. a structure contains no nucleotides."""

    retcode = RetCode.E_NOTHING_TO_CLASSIFY


_ERRORS_BY_RETCODE: dict[RetCode, type[LLKAError]] = {
    cls.retcode: cls
    for cls in (
        InvalidArgumentError,
        MismatchingSizesError,
        LLKANotImplementedError,
        MultipleAltIdsError,
        MissingAtomsError,
        MismatchingDataError,
        BadClassificationClustersError,
        BadGoldenStepsError,
        BadConfalsError,
        BadAverageNuAnglesError,
        BadClassificationLimitsError,
        NoFileError,
        CannotReadFileError,
        BadDataError,
        NoDataError,
        NothingToClassifyError,
    )
}


def error_for_retcode(code: int, message: str) -> LLKAError:
    """Returns an exception instance that corresponds to a libLLKA return code."""
    try:
        retcode = RetCode(code)
    except ValueError:
        return LLKAError(f'{message} (unknown return code {code})')

    cls = _ERRORS_BY_RETCODE.get(retcode, LLKAError)
    return cls(message, retcode)
