import numpy as np

"""
python version of ROOT's TMatrixDSymEigen
"""
class TMatrixDSymEigen:
    def __init__(self, m):
        self.m = m
        self.eigenvalues = None
        self.eigenvectors = None

    def diagonalize(self):
        if self.eigenvalues is not None and self.eigenvectors is not None:
            return
        w, v = np.linalg.eigh(self.m)
        self.eigenvalues = w
        self.eigenvectors = v.T

    def get_eigen_values(self):
        if self.eigenvalues is None:
            self.diagonalize()
        return self.eigenvalues

    def get_eigen_vectors(self):
        if self.eigenvectors is None:
            self.diagonalize()
        return self.eigenvectors
