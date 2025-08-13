import numpy as np
import time
import math
import matplotlib.pyplot as plt

class Complex_operations():
    def __init__(self, re, im):
        self.re = re
        self.im = im

    def __str__(self):
        return f"{self.re:.4f} + {self.im:.4f}i"

    def add(c1, c2):
        return Complex_operations(c1.re + c2.re, c1.im + c2.im)

    def sub(c1, c2):
        return Complex_operations(c1.re - c2.re, c1.im - c2.im)

    def mul(c1, c2):
        re = c1.re*c2.re - c1.im*c2.im
        im = c1.re*c2.im + c1.im*c2.re
        return Complex_operations(re, im)

    def div(c1, c2):
        denom = c2.re**2 + c2.im**2
        re = (c1.re*c2.re + c1.im*c2.im) / denom
        im = (c1.im*c2.re - c1.re*c2.im) / denom
        return Complex_operations(re, im)

    def conjugate(c):
        return Complex_operations(c.re, -c.im)

    def modulus(c):
        return math.sqrt(c.re**2 + c.im**2)

    def argument(c):
        return math.atan2(c.im, c.re)

    def to_polar(c):
        r = Complex_operations.modulus(c)
        theta =Complex_operations.argument(c)
        return (r, theta)

    def from_polar(r, theta):
        return Complex_operations(r * math.cos(theta), r * math.sin(theta))
    
def read_vectors(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        n = int(lines[0])
        A = list(map(float, lines[1].split()))
        B = list(map(float, lines[2].split()))
    return A, B

def read_matrix(filename):
    with open(filename,'r') as f:
        matrix=[]
        for line in f:
            row=[float(num)  for num in line.strip().split()]
            matrix.append(row)
    return matrix

class RNG():
    def __init__(self, seed=1, a=1664525, c=1013904223, m=2**32):
        self.a = a
        self.c = c
        self.m = m
        self.state = seed

    def random(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state / self.m  # Normalize to [0,1)
class VectorAlgebra:
    @staticmethod
    def dot_product(a, b):
        """Return dot product of vectors a and b."""
        return sum(a[i] * b[i] for i in range(len(a)))

    @staticmethod
    def cross_product(a, b):
        """Return cross product of 3D vectors a and b."""
        return [
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        ]

    @staticmethod
    def magnitude(a):
        """Return magnitude of vector a."""
        return math.sqrt(sum(x**2 for x in a))

    @staticmethod
    def normalize(a):
        """Return unit vector in direction of a."""
        mag = VectorAlgebra.magnitude(a)
        return [x/mag for x in a]

    @staticmethod
    def angle_between(a, b):
        """Return angle in radians between vectors a and b."""
        dot = VectorAlgebra.dot_product(a, b)
        mag_a = VectorAlgebra.magnitude(a)
        mag_b = VectorAlgebra.magnitude(b)
        # Clamp value to avoid domain errors in acos
        cos_theta = max(min(dot / (mag_a * mag_b), 1), -1)
        return math.acos(cos_theta)

    @staticmethod
    def projection_of_a_on_b(a, b):
        """Return projection vector of a onto b."""
        mag_b = VectorAlgebra.magnitude(b)
        unit_b = [x/mag_b for x in b]
        scalar_proj = VectorAlgebra.dot_product(a, unit_b)
        return [scalar_proj * x for x in unit_b]
class matrix_operation:
    @staticmethod
    def add(A, B):
        """Return A + B."""
        rows = len(A)
        cols = len(A[0])
        return [[A[i][j] + B[i][j] for j in range(cols)] for i in range(rows)]

    @staticmethod
    def subtract(A, B):
        """Return A - B."""
        rows = len(A)
        cols = len(A[0])
        return [[A[i][j] - B[i][j] for j in range(cols)] for i in range(rows)]

    @staticmethod
    def multiply(A, B):
        """Return A * B (matrix multiplication)."""
        rows_A = len(A)
        cols_A = len(A[0])
        cols_B = len(B[0])
        result = [[0 for _ in range(cols_B)] for _ in range(rows_A)]
        for i in range(rows_A):
            for j in range(cols_B):
                for k in range(cols_A):
                    result[i][j] += A[i][k] * B[k][j]
        return result

    @staticmethod
    def transpose(A):
        """Return transpose of A."""
        rows = len(A)
        cols = len(A[0])
        return [[A[j][i] for j in range(rows)] for i in range(cols)]

    @staticmethod
    def identity(n):
        """Return n x n identity matrix."""
        I = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            I[i][i] = 1
        return I

    @staticmethod
    def determinant(A):
        """Recursive computation of determinant."""
        n = len(A)
        if n == 1:
            return A[0][0]
        if n == 2:
            return A[0][0]*A[1][1] - A[0][1]*A[1][0]
        
        det = 0
        for c in range(n):
            sub_matrix = [row[:c] + row[c+1:] for row in A[1:]]
            det += ((-1)**c) * A[0][c] * matrix_operation.determinant(sub_matrix)
        return det

    @staticmethod
    def inverse_gauss_jordan(A):
        """Return inverse of A using Gauss–Jordan elimination."""
        n = len(A)
        AM = [row[:] for row in A]  # Copy
        I = matrix_operation.identity(n)

        for fd in range(n):
            fdScaler = 1.0 / AM[fd][fd]
            for j in range(n):
                AM[fd][j] *= fdScaler
                I[fd][j] *= fdScaler
            for i in list(range(n)):
                if i != fd:
                    crScaler = AM[i][fd]
                    for j in range(n):
                        AM[i][j] -= crScaler * AM[fd][j]
                        I[i][j] -= crScaler * I[fd][j]
        return I