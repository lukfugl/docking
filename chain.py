from Bio.PDB import *
import numpy
from basis import Basis
from math import sqrt

class Chain:
    def __init__(self, pdbChain, charges):
        self.pdbChain = pdbChain
        self.charges = charges
        self.basis = Basis()
        self.atoms = Selection.unfold_entities(pdbChain, 'A')
        self.searcher = NeighborSearch(self.atoms)
        centroid = numpy.array((0, 0, 0))
        for atom in self.atoms:
            centroid += atom.get_coord()
        centroid /= len(self.atoms)
        self.translate(centroid)

    def translate(self, shift, recalculate_positions=True):
        self.basis.translate(shift)
        if recalculate_positions:
            self.recalculate_positions()

    def rotate_x(self, radians, recalculate_positions=True):
        self.basis.rotate_x(radians)
        if recalculate_positions:
            self.recalculate_positions()

    def rotate_y(self, radians, recalculate_positions=True):
        self.basis.rotate_y(radians)
        if recalculate_positions:
            self.recalculate_positions()

    def rotate_z(self, radians, recalculate_positions=True):
        self.basis.rotate_z(radians)
        if recalculate_positions:
            self.recalculate_positions()

    def recalculate_positions(self):
        self.positions = dict()
        for atom in self.atoms:
            self.positions[atom] = self.basis.convert(atom.get_coord())

    def neighborhood(self, center, radius=12):
        return self.searcher.search(self.basis.deconvert(center), radius)

    def position(self, atom):
        return self.positions[atom]

    def charge(self, atom):
        return self.charges[atom]

    def score_atom(self, other, chain, radius=12):
        p1 = chain.position(other)
        q1 = chain.charge(other)
        neighborhood = self.neighborhood(p1, radius)
        score = 0
        for atom in neighborhood:
            p2 = self.position(atom)
            q2 = self.charge(atom)
            delta = p1 - p2
            r = sqrt(numpy.dot(delta, delta))
            if r > 0:
                score += q2 / r
        return score * q1 / (4 * 3.141592)

    def score_chain(self, chain, radius=12):
        score = 0
        for atom in chain.atoms:
            score += self.score_atom(atom, chain, radius)
        return score
