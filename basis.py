from math import sin, cos
import numpy

class Basis:
    def __init__(self, origin=numpy.array([0, 0, 0]),
                 array=numpy.array([[1, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]])):
        self.origin = origin
        self.array = array
        self.transpose = numpy.transpose(self.array)

    def translate(self, shift):
        return Basis(self.origin + shift, self.array)

    def rotate_x(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[1, 0,  0],
                                [0, c, -s],
                                [0, s,  c]])
        return Basis(self.origin, numpy.dot(rotation, self.array))

    def rotate_y(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[ c, 0, s],
                                [ 0, 1, 0],
                                [-s, 0, c]])
        return Basis(self.origin, numpy.dot(rotation, self.array))

    def rotate_z(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[c, -s, 0],
                                [s,  c, 0],
                                [0,  0, 1]])
        return Basis(self.origin, numpy.dot(rotation, self.array))

    def convert(self, v):
        return numpy.dot(self.array, v - self.origin) + self.origin

    def deconvert(self, v):
        return numpy.dot(self.transpose, v - self.origin) + self.origin
