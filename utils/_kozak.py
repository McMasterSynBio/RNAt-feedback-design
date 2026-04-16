def locate_kozak(seq: str) -> tuple[int, int] | None:
    """
    Locate the start and end positions of the KOZAK sequence of a given RNA sequence.
    """
    kozak_pattern = "AUGG"
    start = seq.find(kozak_pattern)
    if start == -1: return None
    return start, start + len(kozak_pattern)
    