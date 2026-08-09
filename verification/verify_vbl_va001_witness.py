"""Exact replay of the data-grounded VBL-VA001 quadratic witness.

The calibration nodes and signatures are six-decimal rationalizations of five
dominant cospectral directions extracted from six normal-condition records.
No floating-point solver is used in the certificate below.
"""

from fractions import Fraction as F
from itertools import combinations
from math import sqrt


NODES = [F("0.968522"), F("0.876070"), F("-0.296151"), F("-0.524590"), F("-0.720003")]
SIGNATURES = [
    [F("-0.054531"), F("0.996129"), F("-0.068944")],
    [F("0.980644"), F("0.167125"), F("0.102012")],
    [F("0.152303"), F("0.986844"), F("-0.054250")],
    [F("0.869001"), F("-0.350335"), F("0.349432")],
    [F("0.606651"), F("0.716228"), F("-0.344951")],
]

# Symmetric-entry order within each coefficient: 00, 01, 02, 11, 12, 22.
PAIRS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
FIXED = {15: F("-0.316171"), 16: F("-0.537470"), 17: F("-0.909040")}


def det2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def det3(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def solve_square(matrix, rhs):
    """Exact Gauss-Jordan solution of a nonsingular square system."""
    augmented = [list(row) + [value] for row, value in zip(matrix, rhs, strict=True)]
    size = len(augmented)
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            if scale:
                augmented[row] = [
                    value - scale * pivot_value
                    for value, pivot_value in zip(augmented[row], augmented[column], strict=True)
                ]
    return [row[-1] for row in augmented]


def equation_coefficient(variable, degree, output, signature):
    pair = PAIRS[variable % 6]
    if variable // 6 != degree:
        return F(0)
    row, column = pair
    if row == column:
        return signature[column] if output == row else F(0)
    if output == row:
        return signature[column]
    if output == column:
        return signature[row]
    return F(0)


def build_coefficients():
    free_indices = [index for index in range(18) if index not in FIXED]
    matrix = []
    rhs = []
    for node, signature in zip(NODES, SIGNATURES, strict=True):
        for output in range(3):
            row = []
            target = signature[output]
            for variable in free_indices:
                value = sum(
                    node**degree * equation_coefficient(variable, degree, output, signature)
                    for degree in range(3)
                )
                row.append(value)
            for variable, fixed_value in FIXED.items():
                target -= fixed_value * sum(
                    node**degree * equation_coefficient(variable, degree, output, signature)
                    for degree in range(3)
                )
            matrix.append(row)
            rhs.append(target)
    solution = solve_square(matrix, rhs)
    entries = [F(0)] * 18
    for index, value in FIXED.items():
        entries[index] = value
    for index, value in zip(free_indices, solution, strict=True):
        entries[index] = value
    coefficients = []
    for degree in range(3):
        coefficient = [[F(0) for _ in range(3)] for _ in range(3)]
        for offset, (row, column) in enumerate(PAIRS):
            value = entries[6 * degree + offset]
            coefficient[row][column] = value
            coefficient[column][row] = value
        coefficients.append(coefficient)
    return coefficients


def matmul(left, right):
    return [
        [sum(left[row][middle] * right[middle][column] for middle in range(3)) for column in range(3)]
        for row in range(3)
    ]


def matsub(left, right):
    return [[left[row][column] - right[row][column] for column in range(3)] for row in range(3)]


def polynomial(coefficients, node):
    return [
        [sum(node**degree * coefficients[degree][row][column] for degree in range(3)) for column in range(3)]
        for row in range(3)
    ]


def apply(matrix, vector):
    return [sum(matrix[row][column] * vector[column] for column in range(3)) for row in range(3)]


def commutator(left, right):
    return matsub(matmul(left, right), matmul(right, left))


def max_row_sum(matrix):
    return max(sum(abs(value) for value in row) for row in matrix)


def normalized_determinant(rows):
    determinant = abs(det3([[rows[column][row] for column in range(3)] for row in range(3)]))
    norms = [sqrt(float(sum(value * value for value in vector))) for vector in rows]
    return float(determinant) / (norms[0] * norms[1] * norms[2])


def main():
    coefficients = build_coefficients()

    for node, signature in zip(NODES, SIGNATURES, strict=True):
        assert apply(polynomial(coefficients, node), signature) == signature

    normalized_minors = []
    for subset in combinations(range(5), 3):
        rows = [SIGNATURES[index] for index in subset]
        determinant = det3([[rows[column][row] for column in range(3)] for row in range(3)])
        assert determinant != 0
        normalized_minors.append(normalized_determinant(rows))

    commutators = [commutator(coefficients[left], coefficients[right]) for left, right in combinations(range(3), 2)]
    assert any(any(value for row in matrix for value in row) for matrix in commutators)

    # Grid-point certificate: ||P(x_i)|| < 191/200 because
    # (191/200)^2 I - P(x_i)^2 is positive definite by Sylvester's criterion.
    left = F(-1)
    right = F(-929, 1000)
    cell_count = 100
    step = (right - left) / cell_count
    grid_radius = F(191, 200)
    identity = [[F(int(row == column)) for column in range(3)] for row in range(3)]
    for index in range(cell_count + 1):
        value = polynomial(coefficients, left + index * step)
        square = matmul(value, value)
        slack = [
            [grid_radius * grid_radius * identity[row][column] - square[row][column] for column in range(3)]
            for row in range(3)
        ]
        assert slack[0][0] > 0
        assert det2([row[:2] for row in slack[:2]]) > 0
        assert det3(slack) > 0

    # For symmetric matrices ||A||_2 <= ||A||_infinity.  Thus on |x|<=1,
    # ||P'(x)||_2 <= ||C1||_infinity + 2 ||C2||_infinity.
    derivative_bound = max_row_sum(coefficients[1]) + 2 * max_row_sum(coefficients[2])
    continuum_bound = grid_radius + derivative_bound * step / 2
    assert continuum_bound < F(24, 25)

    print("VBL-VA001 exact witness: PASS")
    print(f"minimum normalized 3x3 target minor: {min(normalized_minors):.12f}")
    print(f"derivative operator-norm upper bound: {float(derivative_bound):.12f}")
    print(f"certified continuum upper bound: {float(continuum_bound):.12f} < 24/25")
    print("commuting exact-calibration lower bound: 1 (full-spark obstruction)")


if __name__ == "__main__":
    main()
