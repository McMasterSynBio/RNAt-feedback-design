from data.kozak_pattern import KOZAK

def locate_kozak(seq: str, kozak_pattern: str = KOZAK) -> tuple[int, int] | None:
    """
    Locate the start and end positions of the KOZAK sequence of a given RNA sequence.
    """
    start = seq.rfind(kozak_pattern) # Find from rightmost position
    if start == -1: return None
    return start, start + len(kozak_pattern)
    