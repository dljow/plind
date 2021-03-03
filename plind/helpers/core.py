import numpy as np

# from https://gist.github.com/quantumfx/3cd2c3c4f673baddf76667975ea6589a
def unordered_pairing(a, b):
    return ((a+1) * (b+1)).astype(int) + ((np.abs(a - b) - 1)**2 / 4).astype(int)


def untangle(data):
    weights, points = zip(*data)
    return (
        np.concatenate(points),
        np.repeat(weights, [len(grp) for grp in points]),
    )

def get_all_exponents(dim, max_degree):
    """Get all exponent combinations of dimension `dim` and maximum degree
    `max_degree`. This method is actually meant for evaluating all polynomials
    with these exponents.

    This function is used in the Grundman-Moeller integration scheme and is from

    Source code: quadpy
    Author: Nico Schlomer
    Code version: 0.16.6
    Availability: https://github.com/nschloe/quadpy
    """

    def augment(exponents):
        """This function takes the values and exponents of a given monomial
        level, e.g., [(1,0,0), (0,1,0), (0,0,1)], and augments them by one
        level, i.e., [(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)].
        The methods works for all dimension and is based on the observation
        that the augmentation can happen by

          1. adding 1 to all exponents from the previous level (i.e.,
             multiplication with x[0]),
          2. adding 1 to the second exponent of all exponent tuples where
             ex[0]==0 (i.e., multiplication with x[1]),
          3. adding 1 to the third exponent of all exponent tuples where
             ex[0]==0, ex[1]=0 (i.e., multiplication with x[2]),

        etc. The function call is recursive.
        """

        if len(exponents) == 0 or len(exponents[0]) == 0:
            return []

        idx_leading_zero = [k for k in range(len(exponents)) if exponents[k][0] == 0]
        exponents_with_leading_zero = [exponents[k][1:] for k in idx_leading_zero]
        # val1 = vals[idx_leading_zero]
        # x1 = x[1:]
        out1 = augment(exponents_with_leading_zero)

        # increment leading exponent by 1
        out = [[e[0] + 1] + e[1:] for e in exponents]
        # vals0 = vals * x[0]

        out += [[0] + e for e in out1]
        # out_vals = np.concatenate([vals0, vals1])
        return out

    # dim = x.shape[0]

    # Initialization, level 0
    exponents = [dim * [0]]
    # vals = np.array(np.ones(x.shape[1:]))

    # all_vals = []
    all_exponents = []

    # all_vals.append(vals)
    all_exponents.append(exponents)
    for _ in range(max_degree):
        exponents = augment(exponents)
        # all_vals.append(vals)
        all_exponents.append(exponents)

    return all_exponents
