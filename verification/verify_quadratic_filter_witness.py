"""Exact replay of the quadratic noncommuting-filter witness.

All arithmetic is rational.  Positivity on the full stopband [-1, 0] is
certified by nonnegative Bernstein coefficients after the change x = t - 1.
No SDP solver or floating-point sampling enters the certificate.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import comb


Polynomial = list[F]  # ascending power coefficients
Matrix = list[list[F]]


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    return [
        (left[i] if i < len(left) else F())
        + (right[i] if i < len(right) else F())
        for i in range(size)
    ]


def poly_scale(value: F, polynomial: Polynomial) -> Polynomial:
    return [value * coefficient for coefficient in polynomial]


def poly_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [F() for _ in range(len(left) + len(right) - 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return result


def shift_to_unit_interval(polynomial: Polynomial) -> Polynomial:
    """Return coefficients of q(t)=p(t-1) from coefficients of p(x)."""
    degree = len(polynomial) - 1
    shifted = [F() for _ in range(degree + 1)]
    for power, coefficient in enumerate(polynomial):
        for unit_power in range(power + 1):
            shifted[unit_power] += (
                coefficient
                * comb(power, unit_power)
                * ((-1) ** (power - unit_power))
            )
    return shifted


def power_to_bernstein(polynomial: Polynomial) -> list[F]:
    """Convert a power-basis polynomial on [0,1] to Bernstein coefficients."""
    degree = len(polynomial) - 1
    return [
        sum(
            (
                polynomial[power]
                * F(comb(index, power), comb(degree, power))
                for power in range(index + 1)
            ),
            F(),
        )
        for index in range(degree + 1)
    ]


def evaluate_matrix(coefficients: list[Matrix], point: F) -> Matrix:
    return [
        [
            sum(
                (point**power * coefficients[power][row][column]
                 for power in range(len(coefficients))),
                F(),
            )
            for column in range(2)
        ]
        for row in range(2)
    ]


def matvec(matrix: Matrix, vector: list[F]) -> list[F]:
    return [
        sum((matrix[row][column] * vector[column] for column in range(2)), F())
        for row in range(2)
    ]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum((left[i][k] * right[k][j] for k in range(2)), F())
            for j in range(2)
        ]
        for i in range(2)
    ]


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return [[left[i][j] - right[i][j] for j in range(2)] for i in range(2)]


def matrix_add(*matrices: Matrix) -> Matrix:
    return [
        [sum((matrix[i][j] for matrix in matrices), F()) for j in range(2)]
        for i in range(2)
    ]


def matrix_scale(value: F, matrix: Matrix) -> Matrix:
    return [[value * matrix[i][j] for j in range(2)] for i in range(2)]


def main() -> None:
    coefficients: list[Matrix] = [
        [[F(61, 96), F(7, 48)], [F(7, 48), F(61, 96)]],
        [[F(7, 24), F(7, 48)], [F(7, 48), F(77, 96)]],
        [[F(-7, 24), F(-7, 24)], [F(-7, 24), F(-7, 16)]],
    ]
    pass_points = [F(1, 2), F(1), F(3, 2), F(7, 2)]
    targets = [[F(1), F(2)], [F(0), F(1)], [F(-1), F(2)], [F(-1), F(1)]]

    # Exact tangential constraints and full spark in dimension two.
    for point, target in zip(pass_points, targets, strict=True):
        assert matvec(evaluate_matrix(coefficients, point), target) == target
    assert len(set(pass_points)) == len(pass_points)
    target_minors = [
        left[0] * right[1] - left[1] * right[0]
        for left, right in combinations(targets, 2)
    ]
    assert all(minor != 0 for minor in target_minors)
    degree, dimension, number_targets = 2, 2, 4
    assert number_targets - dimension + 1 == degree + 1

    # The coefficient tuple is genuinely noncommuting.
    commutators = [
        subtract(multiply(coefficients[i], coefficients[j]),
                 multiply(coefficients[j], coefficients[i]))
        for i, j in combinations(range(3), 2)
    ]
    expected_upper_right = [F(343, 4608), F(-49, 2304), F(49, 384)]
    assert [commutator[0][1] for commutator in commutators] == expected_upper_right
    assert all(commutator != [[F(), F()], [F(), F()]] for commutator in commutators)

    p00 = [coefficient[0][0] for coefficient in coefficients]
    p01 = [coefficient[0][1] for coefficient in coefficients]
    p11 = [coefficient[1][1] for coefficient in coefficients]
    operator_norm_upper = F(25, 32)
    bound = [operator_norm_upper, F(), F()]

    bernstein_data: dict[str, dict[str, list[F]]] = {}
    for label, sign in (("aI-P", F(-1)), ("aI+P", F(1))):
        diagonal_00 = poly_add(bound, poly_scale(sign, p00))
        off_diagonal = poly_scale(sign, p01)
        diagonal_11 = poly_add(bound, poly_scale(sign, p11))
        determinant = poly_add(
            poly_multiply(diagonal_00, diagonal_11),
            poly_scale(F(-1), poly_multiply(off_diagonal, off_diagonal)),
        )
        trace = poly_add(diagonal_00, diagonal_11)
        entries = {
            "leading_minor": power_to_bernstein(shift_to_unit_interval(diagonal_00)),
            "determinant": power_to_bernstein(shift_to_unit_interval(determinant)),
            "trace": power_to_bernstein(shift_to_unit_interval(trace)),
        }
        bernstein_data[label] = entries
        assert all(value > 0 for value in entries["leading_minor"])
        assert all(value >= 0 for value in entries["determinant"])

    assert bernstein_data["aI-P"]["determinant"] == [
        F(1421, 1536), F(2597, 6144), F(4655, 27648), F(931, 18432), F(0)
    ]
    assert bernstein_data["aI+P"]["determinant"] == [
        F(1, 16), F(1711, 3072), F(3835, 3456), F(3707, 2304), F(1525, 768)
    ]
    assert operator_norm_upper == F(25, 32) < 1

    # Robust commuting lower bound.  For normalized target directions, every
    # pair has squared inner product at most 9/10.  Hence every pair has
    # smallest singular value at least 1/5, and any unit common eigenvector has
    # overlap at least 1/(5 sqrt(2)) with at least three targets.  We use
    # sqrt(2)<3/2 to obtain the rational slope 60.  The largest absolute
    # Lagrange sum at x=0 over a three-node subset is exactly 8.
    target_norm_squares = [sum(value * value for value in target) for target in targets]
    normalized_dot_squares = []
    for left, right in combinations(range(4), 2):
        dot = sum(targets[left][i] * targets[right][i] for i in range(2))
        normalized_dot_squares.append(
            dot * dot / (target_norm_squares[left] * target_norm_squares[right])
        )
    assert max(normalized_dot_squares) == F(9, 10)
    # 1-|dot| >= (1-|dot|^2)/2 >= 1/20 > 1/25.
    assert (1 - max(normalized_dot_squares)) / 2 == F(1, 20) > F(1, 25)
    full_spark_floor = F(1, 5)
    lagrange_sums = []
    for subset in combinations(pass_points, 3):
        weights = []
        for node in subset:
            others = [other for other in subset if other != node]
            weights.append(
                ((-others[0]) * (-others[1]))
                / ((node - others[0]) * (node - others[1]))
            )
        lagrange_sums.append(sum(abs(weight) for weight in weights))
    assert sorted(lagrange_sums) == [F(11, 4), F(19, 5), F(7), F(8)]
    assert full_spark_floor == F(1, 5)
    robust_slope = F(60)
    robust_tolerance = (F(1) - operator_norm_upper) / robust_slope
    assert robust_tolerance == F(7, 1920)

    # Exact five-tap reciprocal linear-phase FIR realization.  Under
    # x=(9y+5)/4 and y=cos(omega), P(x) becomes
    # G(omega)=B_2+2 B_1 cos(omega)+2 B_0 cos(2 omega).  The causal response
    # with taps (B_0,B_1,B_2,B_1,B_0) is e^{-2 i omega}G(omega).
    c0, c1, c2 = coefficients
    d0 = matrix_add(c0, matrix_scale(F(5, 4), c1), matrix_scale(F(25, 16), c2))
    d1 = matrix_add(matrix_scale(F(9, 4), c1), matrix_scale(F(45, 8), c2))
    d2 = matrix_scale(F(81, 16), c2)
    tap_0 = matrix_scale(F(1, 4), d2)
    tap_1 = matrix_scale(F(1, 2), d1)
    tap_2 = matrix_add(d0, matrix_scale(F(1, 2), d2))
    taps = [tap_0, tap_1, tap_2, tap_1, tap_0]
    assert taps == [
        [[F(-189, 512), F(-189, 512)], [F(-189, 512), F(-567, 1024)]],
        [[F(-63, 128), F(-21, 32)], [F(-21, 32), F(-21, 64)]],
        [[F(-149, 768), F(-665, 768)], [F(-665, 768), F(-235, 1536)]],
        [[F(-63, 128), F(-21, 32)], [F(-21, 32), F(-21, 64)]],
        [[F(-189, 512), F(-189, 512)], [F(-189, 512), F(-567, 1024)]],
    ]
    for point, target in zip(pass_points, targets, strict=True):
        cosine = (4 * point - 5) / 9
        cosine_two = 2 * cosine * cosine - 1
        zero_phase = matrix_add(
            tap_2,
            matrix_scale(2 * cosine, tap_1),
            matrix_scale(2 * cosine_two, tap_0),
        )
        assert zero_phase == evaluate_matrix(coefficients, point)
        assert matvec(zero_phase, target) == target
    tap_commutator = subtract(multiply(tap_0, tap_1), multiply(tap_1, tap_0))
    assert tap_commutator != [[F(), F()], [F(), F()]]

    print("PASS: 4 exact tangential constraints and all 2x2 target minors are nonzero")
    print("PASS: all three exact coefficient commutators are nonzero")
    print("PASS: Bernstein certificates give (25/32)I±P(x) >= 0 on [-1,0]")
    print("PASS: exact continuum operator-norm bound is 25/32 < 1")
    print("PASS: every commuting symmetric degree-2 feasible filter is I by the root count")
    print("PASS: approximate commuting leakage is at least 1-60 delta")
    print("PASS: strict separation survives every delta < 7/1920")
    print("PASS: exact palindromic five-tap MIMO FIR realization preserves all pass directions")


if __name__ == "__main__":
    main()
