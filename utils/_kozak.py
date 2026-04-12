def locate_kozak(seq: str) -> tuple[int, int] | None:
    """
    Locate the start and end positions of the KOZAK sequence of a given RNA sequence.
    """
    kozak_pattern = "GCCRCCAUGG"  # R = A or G
    for i in range(len(seq) - len(kozak_pattern) + 1):
        match = True
        for j, p in enumerate(kozak_pattern):
            if p == "R":
                if seq[i + j] not in ("A", "G"):
                    match = False
                    break
            else:
                if seq[i + j] != p:
                    match = False
                    break
        if match:
            return i, i + len(kozak_pattern)
    return None
    