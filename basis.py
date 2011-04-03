from math import sin, cos
import numpy

class Basis:
    def __init__(self):
        self.origin = numpy.array([0, 0, 0])
        self.array = numpy.array([[1, 0, 0],
                                  [0, 1, 0],
                                  [0, 0, 1]])
        self.transpose = numpy.transpose(self.array)

    def translate(self, shift):
        self.origin += shift

    def rotate_x(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[1, 0,  0],
                                [0, c, -s],
                                [0, s,  c]])
        self.array = numpy.dot(rotation, self.array)
        self.transpose = numpy.transpose(self.array)

    def rotate_y(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[ c, 0, s],
                                [ 0, 1, 0],
                                [-s, 0, c]])
        self.array = numpy.dot(rotation, self.array)
        self.transpose = numpy.transpose(self.array)

    def rotate_z(self, radians):
        c = cos(radians)
        s = sin(radians)
        rotation = numpy.array([[c, -s, 0],
                                [s,  c, 0],
                                [0,  0, 1]])
        self.array = numpy.dot(rotation, self.array)
        self.transpose = numpy.transpose(self.array)

    def convert(self, v):
        return numpy.dot(self.array, v - self.origin) + self.origin

    def deconvert(self, v):
        return numpy.dot(self.transpose, v - self.origin) + self.origin
