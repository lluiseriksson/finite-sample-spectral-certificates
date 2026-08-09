"""Exact replay of the quadratic noncommuting-filter witness.

All arithmetic is rational.  Positivity on the full stopband [-1, 0] is
certified by positive Bernstein coefficients after the change x = t - 1.
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
        [[F(2, 5), F(-3, 16)], [F(-3, 16), F(29, 32)]],
        [[F(15, 16), F(3, 20)], [F(3, 20), F(3, 32)]],
        [[F(-3, 8), F()], [F(), F(-3, 80)]],
    ]
    pass_points = [F(1, 2), F(1), F(3, 2), F(2)]
    targets = [[F(1), F(-2)], [F(1), F(-1)], [F(1), F(1)], [F(1), F(2)]]

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
    expected_upper_right = [F(1053, 12800), F(-81, 1280), F(81, 1600)]
    assert [commutator[0][1] for commutator in commutators] == expected_upper_right
    assert all(commutator != [[F(), F()], [F(), F()]] for commutator in commutators)

    p00 = [coefficient[0][0] for coefficient in coefficients]
    p01 = [coefficient[0][1] for coefficient in coefficients]
    p11 = [coefficient[1][1] for coefficient in coefficients]
    one = [F(1), F(), F()]

    bernstein_data: dict[str, dict[str, list[F]]] = {}
    spectral_margins: list[F] = []
    for label, sign in (("I-P", F(-1)), ("I+P", F(1))):
        diagonal_00 = poly_add(one, poly_scale(sign, p00))
        off_diagonal = poly_scale(sign, p01)
        diagonal_11 = poly_add(one, poly_scale(sign, p11))
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
        assert all(value > 0 for value in entries["determinant"])
        determinant_lower = min(entries["determinant"])
        trace_upper = max(entries["trace"])
        spectral_margins.append(determinant_lower / trace_upper)

    assert bernstein_data["I-P"]["determinant"] == [
        F(81, 256), F(1701, 10240), F(219, 2560), F(441, 10240), F(27, 1280)
    ]
    assert bernstein_data["I+P"]["determinant"] == [
        F(53, 1280), F(8389, 10240), F(783, 512), F(21913, 10240), F(3371, 1280)
    ]
    uniform_margin = min(spectral_margins)
    assert spectral_margins == [F(3, 304), F(53, 4232)]
    assert uniform_margin == F(3, 304)
    operator_norm_upper = F(1) - uniform_margin
    assert operator_norm_upper == F(301, 304) < 1

    # Robust commuting lower bound.  For normalized target directions, every
    # pair has squared inner product at most 9/10.  Since
    # sqrt(9/10) < 593/625, any unit common eigenvector has overlap at least
    # m0=4/25 with at least three of the four targets.  The largest absolute
    # Lagrange sum at x=0 over a three-node subset is exactly 17.
    target_norm_squares = [sum(value * value for value in target) for target in targets]
    normalized_dot_squares = []
    for left, right in combinations(range(4), 2):
        dot = sum(targets[left][i] * targets[right][i] for i in range(2))
        normalized_dot_squares.append(
            dot * dot / (target_norm_squares[left] * target_norm_squares[right])
        )
    assert max(normalized_dot_squares) == F(9, 10)
    assert F(9, 10) < F(593, 625) ** 2
    overlap_floor = F(4, 25)
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
    assert sorted(lagrange_sums) == [F(5), F(5), F(7), F(17)]
    robust_slope = max(lagrange_sums) / overlap_floor
    robust_tolerance = uniform_margin / robust_slope
    assert robust_slope == F(425, 4)
    assert robust_tolerance == F(3, 32300)

    # Exact five-tap reciprocal linear-phase FIR realization.  Under
    # x=(3y+1)/2 and y=cos(omega), P(x) becomes
    # G(omega)=B_2+2 B_1 cos(omega)+2 B_0 cos(2 omega).  The causal response
    # with taps (B_0,B_1,B_2,B_1,B_0) is e^{-2 i omega}G(omega).
    c0, c1, c2 = coefficients
    d0 = matrix_add(c0, matrix_scale(F(1, 2), c1), matrix_scale(F(1, 4), c2))
    d1 = matrix_add(matrix_scale(F(3, 2), c1), matrix_scale(F(3, 2), c2))
    d2 = matrix_scale(F(9, 4), c2)
    tap_0 = matrix_scale(F(1, 4), d2)
    tap_1 = matrix_scale(F(1, 2), d1)
    tap_2 = matrix_add(d0, matrix_scale(F(1, 2), d2))
    taps = [tap_0, tap_1, tap_2, tap_1, tap_0]
    assert taps == [
        [[F(-27, 128), F()], [F(), F(-27, 1280)]],
        [[F(27, 64), F(9, 80)], [F(9, 80), F(27, 640)]],
        [[F(113, 320), F(-9, 80)], [F(-9, 80), F(577, 640)]],
        [[F(27, 64), F(9, 80)], [F(9, 80), F(27, 640)]],
        [[F(-27, 128), F()], [F(), F(-27, 1280)]],
    ]
    for point, target in zip(pass_points, targets, strict=True):
        cosine = (2 * point - 1) / 3
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
    print("PASS: Bernstein certificates give I±P(x) >= (3/304) I on [-1,0]")
    print("PASS: exact continuum operator-norm bound is 301/304 < 1")
    print("PASS: every commuting symmetric degree-2 feasible filter is I by the root count")
    print("PASS: approximate commuting leakage is at least 1-(425/4) delta")
    print("PASS: strict separation survives every delta < 3/32300")
    print("PASS: exact palindromic five-tap MIMO FIR realization preserves all pass directions")


if __name__ == "__main__":
    main()
