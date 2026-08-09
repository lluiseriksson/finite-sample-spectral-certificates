"""Arithmetic replay for the fixed three-channel separation theorem.

The manuscript is the proof.  This standard-library check independently
replays its rational node geometry, dimension counts, exponential rate, and
the conservative perturbation envelopes for a broad finite range.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import ceil


def acosh(value: Decimal) -> Decimal:
    return (value + (value * value - 1).sqrt()).ln()


def group_sizes(total: int) -> list[int]:
    quotient, remainder = divmod(total, 3)
    return [quotient + (index < remainder) for index in range(3)]


def exact_geometry(max_degree: int = 128) -> None:
    for degree in range(7, max_degree + 1):
        total = degree + 3
        sizes = group_sizes(total)
        largest = max(sizes)
        assert largest == ceil(total / 3)
        assert min(sizes) >= largest - 1
        assert 2 * largest <= degree + 1

        groups: list[list[F]] = []
        for group, size in enumerate(sizes):
            groups.append(
                [
                    F(1) + F(5 * (3 * index + group), 6 * largest)
                    for index in range(size)
                ]
            )
        nodes = [node for group in groups for node in group]
        assert len(nodes) == total and len(set(nodes)) == total
        assert min(nodes) >= 1 and max(nodes) < F(7, 2)
        global_separation = min(
            abs(left - right)
            for index, left in enumerate(nodes)
            for right in nodes[index + 1 :]
        )
        assert global_separation >= F(5, 6 * largest)
        for group in groups:
            if len(group) > 1:
                assert {group[i + 1] - group[i] for i in range(len(group) - 1)} == {
                    F(5, 2 * largest)
                }


def envelope_check(max_degree: int = 128) -> tuple[Decimal, Decimal, int]:
    with localcontext() as context:
        context.prec = 80
        eta0 = acosh(Decimal(3))
        eta1 = acosh(Decimal(8))
        delta_eta = eta1 - eta0
        e = Decimal(1).exp()
        lagrange_constant = Decimal(36) * e / 5
        rate = (2 * eta0 - lagrange_constant.ln()) / 3
        assert rate > 0
        base_constant = 2 * ((2 * (lagrange_constant.ln() + eta0) / 3).exp())
        total_constant = base_constant + 1
        noncommuting_threshold = ceil(float(total_constant.ln() / rate))

        worst_neumann = Decimal(0)
        worst_correction_ratio = Decimal(0)
        for degree in range(7, max_degree + 1):
            n = Decimal(degree)
            largest = Decimal(ceil((degree + 3) / 3))
            right_inverse = (
                6
                * largest
                * (Decimal(27) * largest / 5) ** (2 * largest - 1)
            )
            epsilon = (-n**3).exp()
            neumann = Decimal(3).sqrt() * epsilon * right_inverse
            pass_bound = (
                2
                * (4 * e) ** (largest - 1)
                * (n * delta_eta).exp()
            )
            correction = 2 * Decimal(3).sqrt() * right_inverse * epsilon * (
                1 + pass_bound
            )
            correction_ratio = correction * (rate * n).exp()
            assert neumann < Decimal("0.5")
            assert correction_ratio < 1
            worst_neumann = max(worst_neumann, neumann)
            worst_correction_ratio = max(worst_correction_ratio, correction_ratio)

        assert total_constant * (-rate * Decimal(noncommuting_threshold)).exp() < 1
        return worst_neumann, worst_correction_ratio, noncommuting_threshold


def main() -> None:
    exact_geometry()
    neumann, correction_ratio, threshold = envelope_check()
    print("constant-channel scaling arithmetic: PASS")
    print("dimension d=3, signatures M=N+3, exact node data checked for 7<=N<=128")
    print(f"largest Neumann product = {neumann:.8E}")
    print(f"largest correction/base-scale ratio = {correction_ratio:.8E}")
    print(f"conservative noncommuting threshold from C exp(-cN)<1: N={threshold}")


if __name__ == "__main__":
    main()
