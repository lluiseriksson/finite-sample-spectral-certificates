"""Arithmetic replay for the robust fixed-three-channel separation theorem.

The manuscript is the proof.  This standard-library check independently
replays its rational node geometry, dimension counts, exponential rate, and
the conservative perturbation, full-spark, and noise-scale envelopes.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations, permutations
from math import ceil, comb


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


def polynomial_product(left: list[F], right: list[F]) -> list[F]:
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def determinant_polynomial(groups: list[int], indices: list[int], total: int) -> list[F]:
    """Return det(e_g + epsilon*(1,t,t^2)) as coefficients in epsilon."""

    columns: list[list[list[F]]] = []
    for group, index in zip(groups, indices, strict=True):
        t = F(index, total + 1)
        moment = [F(1), t, t * t]
        columns.append(
            [[F(row == group), moment[row]] for row in range(3)]
        )

    coefficients = [F(0)] * 4
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        sign = -1 if inversions % 2 else 1
        term = [F(1)]
        for column, row in enumerate(permutation):
            term = polynomial_product(term, columns[column][row])
        for power, value in enumerate(term):
            coefficients[power] += sign * value
    return coefficients


def quantitative_full_spark_check(max_degree: int = 24) -> F:
    """Exact finite replay of the determinant lower-envelope argument."""

    smallest_ratio: F | None = None
    for degree in range(7, max_degree + 1):
        total = degree + 3
        sizes = group_sizes(total)
        group_of = [group for group, size in enumerate(sizes) for _ in range(size)]
        epsilon = F(1, 2 ** (13 * degree))
        assert epsilon <= F(1, 108 * (total + 1) ** 6)
        lower = epsilon**3 / (2 * (total + 1) ** 6)

        for triple in combinations(range(1, total + 1), 3):
            groups = [group_of[index - 1] for index in triple]
            coefficients = determinant_polynomial(groups, list(triple), total)
            assert coefficients[3] != 0
            nonzero = [abs(value) for value in coefficients if value]
            assert min(nonzero) >= F(1, (total + 1) ** 6)
            assert max(nonzero) <= 18
            determinant = sum(
                value * epsilon**power
                for power, value in enumerate(coefficients)
            )
            assert abs(determinant) >= lower
            ratio = abs(determinant) / lower
            smallest_ratio = ratio if smallest_ratio is None else min(smallest_ratio, ratio)

    assert smallest_ratio is not None
    return smallest_ratio


def envelope_check(max_degree: int = 512) -> tuple[Decimal, Decimal, int, Decimal]:
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
        robust_separation_threshold = ceil(
            float((total_constant.ln() - Decimal("0.75").ln()) / rate)
        )

        right_inverse_base = Decimal(108) * e / 5
        worst_neumann = Decimal(0)
        worst_correction_ratio = Decimal(0)
        smallest_robust_log10 = Decimal("Infinity")
        for degree in range(7, max_degree + 1):
            n = Decimal(degree)
            largest = Decimal(ceil((degree + 3) / 3))
            right_inverse = 6 * right_inverse_base ** (2 * largest - 1)
            epsilon = Decimal(2) ** (-13 * degree)
            neumann = Decimal(3).sqrt() * epsilon * right_inverse
            pass_bound = (
                2
                * (4 * e) ** (largest - 1)
                * (n * delta_eta).exp()
            )
            assert right_inverse <= Decimal(2) ** (4 * degree + 17)
            assert pass_bound <= Decimal(2) ** (Decimal(10 * degree + 11) / 3)
            correction = 2 * Decimal(3).sqrt() * right_inverse * epsilon * (
                1 + pass_bound
            )
            correction_ratio = correction * (rate * n).exp()
            assert neumann < Decimal("0.5")
            assert correction_ratio < 1
            assert neumann <= Decimal(2) ** (18 - 9 * degree)
            assert correction <= Decimal(2) ** (
                Decimal(71 - 17 * degree) / 3
            )
            worst_neumann = max(worst_neumann, neumann)
            worst_correction_ratio = max(worst_correction_ratio, correction_ratio)

            full_spark_floor = epsilon**3 / (48 * (n + 4) ** 6)
            lagrange_at_zero = (Decimal(24) * e / 5) ** degree
            robust_delta = full_spark_floor / (
                4 * Decimal(3).sqrt() * lagrange_at_zero
            )
            smallest_robust_log10 = min(smallest_robust_log10, robust_delta.log10())

        assert (
            total_constant
            * (-rate * Decimal(robust_separation_threshold)).exp()
            < Decimal("0.75")
        )
        return (
            worst_neumann,
            worst_correction_ratio,
            robust_separation_threshold,
            smallest_robust_log10,
        )


def binomial_tail(degree: int, probability: F, threshold: int) -> F:
    return sum(
        F(comb(degree, successes))
        * probability**successes
        * (1 - probability) ** (degree - successes)
        for successes in range(threshold, degree + 1)
    )


def scalar_barrier_check(max_degree: int = 128) -> Decimal:
    """Check the endpoint binomial tails against the Hoeffding envelope."""

    worst_ratio = Decimal(0)
    with localcontext() as context:
        context.prec = 80
        for degree in range(7, max_degree + 1):
            threshold = ceil(degree / 3)
            stop_tail = binomial_tail(degree, F(2, 9), threshold)
            pass_error = 1 - binomial_tail(degree, F(4, 9), threshold)
            bound = (-Decimal(2 * degree) / 81).exp()
            for exact_value in (stop_tail, pass_error):
                value = Decimal(exact_value.numerator) / Decimal(exact_value.denominator)
                assert value <= bound
                worst_ratio = max(worst_ratio, value / bound)
    return worst_ratio


def main() -> None:
    exact_geometry()
    determinant_ratio = quantitative_full_spark_check()
    neumann, correction_ratio, threshold, robust_log10 = envelope_check()
    barrier_ratio = scalar_barrier_check()
    print("robust constant-channel scaling arithmetic: PASS")
    print("dimension d=3, signatures M=N+3, exact node data checked for 7<=N<=128")
    print("rational perturbation epsilon_N = 2^(-13N)")
    print(f"largest Neumann product = {neumann:.8E}")
    print(f"largest correction/base-scale ratio = {correction_ratio:.8E}")
    print(f"smallest exact determinant/lower-envelope ratio = {float(determinant_ratio):.6g}")
    print(f"smallest certified log10 calibration radius through N=512 = {robust_log10:.6f}")
    print(f"largest exact binomial-tail/Hoeffding ratio = {barrier_ratio:.8E}")
    print(f"conservative strict robust-separation threshold: N={threshold}")


if __name__ == "__main__":
    main()
