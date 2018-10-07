import numpy as nm


class DGBasis:

    def __init__(self, degree):
        self.degree = degree
        self.pdegree = degree + 1 # practical degree used in calculations

    def values(self, points):

        if isinstance(points, (int, float)):
            n = 1
        else:
            n = len(points)
        values = nm.ones((self.pdegree, n))
        for i in range(1, self.pdegree):
            values[i] = points * values[i-1]
        return values.T


if __name__ == '__main__':
    bs = DGBasis(5)
    vals = bs.values(nm.linspace(0, 1))
    from matplotlib import pylab as plt

    plt.plot(nm.linspace(0, 1), vals)
    plt.show()

