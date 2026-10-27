import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    n_sites = 0
    n_diff = 0

    # count the gap and the differents sites.
    for a, b in zip(reference_sequence, distant_sequence):
        a, b = a.upper(), b.upper()
        # skip gap
        if a == "-" or b == "-":
            continue
        n_sites += 1
        if a != b:
            n_diff += 1

    if n_sites == 0:
        return 0.0
    # caculate the P;
    p = n_diff / n_sites
    if p == 0:
        return 0.0
    if p >= 0.75:
        return float("inf")
    # JC
    d_site = -0.75 * math.log(1 - 4 * p / 3)

    return d_site * n_sites


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    n_sites = 0
    transitions = 0
    transversions = 0

    for a, b in zip(reference_sequence, distant_sequence):
        a, b = a.upper(), b.upper()

        # skip the gaps
        if a not in "AGTC" or b not in "AGTC":
            continue

        n_sites += 1
        if a == b:
            continue

        # transitions: A<->G, C<->T
        if(a, b) in {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}:
            transitions += 1
        else:
            transversions += 1

    if n_sites == 0:
        return 0.0

    p = transitions / n_sites
    q = transversions / n_sites

    # saturations
    if 1 - 2 * p - q <= 0 or 1 - 2 * q <= 0:
        return float("inf")

    d_site = -0.5 * math.log(1 - 2 * p - q) - 0.25 * math.log(1 - 2 * q)

    return d_site * n_sites

    # raise NotImplementedError()
