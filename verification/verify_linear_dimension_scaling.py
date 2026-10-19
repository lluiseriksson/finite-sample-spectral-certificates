"""Replay the arithmetic behind the sharpened d=N scaling theorem.

This standard-library check is not a substitute for the manuscript proof.
It verifies the exact node counts and separation, the right-inverse constant,
and the two numerical inequalities at N=5.  The logarithmic derivatives
recorded below are then negative for every real N>=5, so those inequalities
remain valid for all larger integer N.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction as F


def acosh(value: Decimal) -> Decimal:
    return (value + (value * value - 1).sqrt()).ln()


def exact_combinatorics(max_degree: int = 64) -> None:
    for degree in range(5, max_degree + 1):
        dimension = degree
        nodes = [F(1) + F(r, 10 * degree) for r in range(1, degree + 1)]
        nodes += [F(2) + F(r, 10 * degree) for r in range(1, degree + 1)]
        assert len(nodes) == 2 * degree == dimension + degree
        assert len(set(nodes)) == len(nodes)
        assert min(nodes) >= 1 and max(nodes) <= F(7, 2)
        separation = min(
            abs(left - right)
            for index, left in enumerate(nodes)
            for right in nodes[index + 1 :]
        )
        assert separation >= F(1, 10 * degree)

        # Two nodes per coordinate; an off-diagonal entry sees four.
        assert 2 + 2 <= degree + 1
        # Full spark leaves at least M-d+1=N+1 nonorthogonal signatures.
        assert 2 * degree - dimension + 1 == degree + 1

        right_inverse_constant = 4 * (45 * degree) ** 3 * dimension
        assert right_inverse_constant == 364_500 * degree**4


def analytic_tail_check() -> tuple[Decimal, Decimal]:
    with localcontext() as context:
        context.prec = 80
        n = Decimal(5)
        eta0 = acosh(Decimal(3))
        eta1 = acosh(Decimal(8))
        delta_eta = eta1 - eta0

        neumann = Decimal(364_500) * n ** Decimal("4.5") * (-n**3).exp()
        correction = (
            Decimal(8_019_000)
            * n ** Decimal("4.5")
            * (-n**3 + n * delta_eta).exp()
        )
        base = Decimal(52) / 5 * eta0.exp() * (-eta0 * n).exp()
        assert neumann <= Decimal("0.5")
        assert correction <= base

        # d/dN log(neumann) = 9/(2N)-3N^2.
        assert Decimal(9) / (2 * n) - 3 * n**2 < 0
        # d/dN log(correction/base) = 9/(2N)-3N^2+eta1.
        assert eta1 < 3
        assert Decimal(9) / (2 * n) - 3 * n**2 + eta1 < 0
        return neumann, correction / base


def main() -> None:
    exact_combinatorics()
    neumann, correction_ratio = analytic_tail_check()
    print("linear-dimension scaling arithmetic: PASS")
    print("dimension=d=N, signatures=M=2N, checked exact node data for 5<=N<=64")
    print(f"N=5 Neumann product = {neumann:.8E}")
    print(f"N=5 correction/base = {correction_ratio:.8E}")
    print("both logarithmic derivatives are negative for every N>=5")


if __name__ == "__main__":
    main()
