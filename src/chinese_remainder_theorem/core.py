from __future__ import annotations

from typing import List, Sequence, Tuple


class CRTInputError(ValueError):
    """Raised when the inputs to :func:`crt` are structurally invalid.

    A dedicated subclass (rather than a bare :class:`ValueError`) lets callers
    filter CRT-specific validation failures out of unrelated value errors that
    might surface from elsewhere in their own code.
    """


def _extended_gcd(a: int, b: int) -> Tuple[int, int, int]:
    """Return ``(g, x, y)`` such that ``a*x + b*y == g == gcd(a, b)``.

    Uses the iterative form of the extended Euclidean algorithm. The
    recursive form is equally correct but risks ``RecursionError`` on large
    coprime moduli where the stack depth matches the number of Euclidean
    steps; the loop cannot overflow.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    # old_r is the gcd (non-negative by construction for non-negative inputs).
    return old_r, old_s, old_t


def crt(residues: Sequence[int], moduli: Sequence[int]) -> int:
    """Solve a system of simultaneous congruences via the Chinese Remainder Theorem.

    Given residues ``r_i`` and pairwise-coprime moduli ``m_i``, returns the
    unique ``x`` with ``0 <= x < prod(m_i)`` such that ``x == r_i (mod m_i)``
    for every ``i``.

    Implementation:
        Combines congruences two at a time. To merge ``x == a (mod m)`` and
        ``x == b (mod n)`` with ``gcd(m, n) == 1``, solve via the extended
        Euclidean algorithm: find ``s, t`` with ``s*m + t*n == 1``; then
        ``x == a + m * ((b - a) * s mod n) (mod m*n)``. Iterating merges all
        congruences and the running modulus is exactly the product of the
        merged moduli, so the final result lies in ``[0, prod(m_i))``.

    Args:
        residues: Remainders for each congruence. May be negative or larger
            than the corresponding modulus; each is normalised internally.
        moduli: Positive integer moduli. Must be non-empty and pairwise
            coprime (``gcd(m_i, m_j) == 1`` for ``i != j``).

    Returns:
        The unique solution in ``[0, prod(moduli))``.

    Raises:
        CRTInputError: If ``residues`` and ``moduli`` have differing or zero
            length, if any modulus is not a positive integer, or if the
            moduli are not pairwise coprime.
    """
    if len(residues) != len(moduli):
        raise CRTInputError(
            f"residues and moduli must have the same length "
            f"(got {len(residues)} and {len(moduli)})"
        )
    if len(moduli) == 0:
        raise CRTInputError("at least one congruence is required")

    for m in moduli:
        if not isinstance(m, int) or isinstance(m, bool):
            raise CRTInputError(f"modulus must be an int, got {m!r}")
        if m <= 0:
            raise CRTInputError(f"modulus must be positive, got {m}")

    for r in residues:
        if not isinstance(r, int) or isinstance(r, bool):
            raise CRTInputError(f"residue must be an int, got {r!r}")

    # Check pairwise coprimality once, up front, so we can raise a clear,
    # specific error before doing any work rather than a confusing one mid-merge.
    n = len(moduli)
    for i in range(n):
        for j in range(i + 1, n):
            if _extended_gcd(moduli[i], moduli[j])[0] != 1:
                raise CRTInputError(
                    f"moduli must be pairwise coprime; "
                    f"gcd({moduli[i]}, {moduli[j]}) != 1"
                )

    # Seed with the first congruence, normalising the residue into [0, m).
    modulus = moduli[0]
    solution = residues[0] % modulus

    for k in range(1, n):
        next_modulus = moduli[k]
        next_residue = residues[k] % next_modulus

        # gcd(m, n) == 1 is guaranteed by the pairwise check above, so the
        # Bézout coefficients always yield the inverse we need.
        g, s, t = _extended_gcd(modulus, next_modulus)
        # g is 1 here; kept in the tuple only because _extended_gcd returns it.
        _ = g

        # x = a + m * ((b - a) * s mod n), with the whole thing taken mod m*n.
        diff = next_residue - solution
        factor = (diff * s) % next_modulus
        solution = solution + modulus * factor
        modulus = modulus * next_modulus
        solution = solution % modulus

    return solution
