
__all__ = (
    'GVPError',
    'ShouldOperateThroughGeneViewError',
    'ExonInATranscriptOnDifferentStrandError',
    'ExonInATranscriptOnDifferentChromosomeError',
    'DuplicatedExonError',
    'TranscriptInAGeneOnDifferentStrandError',
    'DuplicatedTranscriptIDError',
    'DuplicatedTranscriptError',
    'TranscriptInAGeneOnDifferentChromosomeError'
)

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
