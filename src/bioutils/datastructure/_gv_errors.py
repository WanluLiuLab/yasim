_all = [
    'GVPError',
    'ShouldOperateThroughGeneViewError',
    'ExonInATranscriptOnDifferentStrandError',
    'ExonInATranscriptOnDifferentChromosomeError',
    'DuplicatedExonError',
    'TranscriptInAGeneOnDifferentStrandError',
    'DuplicatedTranscriptIDError',
    'DuplicatedTranscriptError',
    'TranscriptInAGeneOnDifferentChromosomeError'
]

__all__ = _all

__all__.append('_all')


class GVPError(ValueError):
    pass


class ShouldOperateThroughGeneViewError(GVPError):
    pass


class ExonInATranscriptOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedExonError(GVPError):
    pass


class ExonInATranscriptOnDifferentStrandError(GVPError):
    pass


class TranscriptInAGeneOnDifferentChromosomeError(GVPError):
    pass


class DuplicatedTranscriptError(GVPError):
    pass


class DuplicatedTranscriptIDError(GVPError):
    pass


class TranscriptInAGeneOnDifferentStrandError(GVPError):
    pass
