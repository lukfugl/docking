from Bio.PDB import *
import numpy
from basis import Basis
from math import sqrt, pi, sin, cos
import string
from pymol import cmd
import random

class Chain:
    def __init__(self, pdbChain, charges):
        self.pdbChain = pdbChain
        self.charges = charges
        self.atoms = Selection.unfold_entities(pdbChain, 'A')
        self.searcher = NeighborSearch(self.atoms)
        self.reset_basis()

    def reset_basis(self):
        centroid = numpy.array((0, 0, 0))
        for atom in self.atoms:
            centroid += atom.get_coord()
        centroid /= len(self.atoms)
        self.set_basis(Basis().translate(centroid))
        self.recalculate_positions()

    def set_basis(self, basis):
        self.basis = basis
        self.recalculate_positions()

    def calculate_positions(self, basis=None):
        if not basis:
            basis = self.basis
        positions = dict()
        for atom in self.atoms:
            positions[atom] = basis.convert(atom.get_coord())
        return positions

    def recalculate_positions(self):
        self.positions = self.calculate_positions()

    def neighborhood(self, center, radius=12, basis=None):
        if not basis:
            basis = self.basis
        return self.searcher.search(basis.deconvert(center), radius)

    def position(self, atom, positions=None):
        if not positions:
            positions = self.positions
        return positions[atom]

    def charge(self, atom):
        return self.charges[atom]

    def atompdb(self, atom):
        residue = atom.get_parent()
        chain = self.pdbChain
        record = 'HETATM'
        coords = self.position(atom)
        return "%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f\n" \
            % (record, atom.get_serial_number(), atom.get_name(),
               residue.get_resname(), chain.get_id()[0], residue.get_id()[1],
               coords[0], coords[1], coords[2])

    def pdbstr(self):
        pdb = []
        for atom in self.atoms:
            pdb.append(self.atompdb(atom))
        return string.join(pdb, '')

    def push2pymol(self, name, state=1):
        cmd.read_pdbstr(self.pdbstr(), name, state)
        
    def score_atom(self, other, chain, radius=12, positions=None, basis=None):
        if not positions:
            positions = self.positions
        if not basis:
            basis = self.basis
        p1 = chain.position(other)
        q1 = chain.charge(other)
        neighborhood = self.neighborhood(p1, radius, basis=basis)
        score = 0
        for atom in neighborhood:
            p2 = self.position(atom, positions=positions)
            q2 = self.charge(atom)
            delta = p1 - p2
            r = sqrt(numpy.dot(delta, delta))
            if r > 0:
                invr = 1 / r
                # coulombic force
                score += q1 * q2 * invr / (4 * pi)
                # vanderwaal force
                score += invr ** 12 - 2 * invr ** 6
        return score

    def score_chain(self, chain, radius=12, basis=None):
        if not basis:
            basis = self.basis
        if basis == self.basis:
            positions = self.positions
        else:
            positions = self.calculate_positions(basis)
        score = 0
        for atom in chain.atoms:
            score += self.score_atom(atom, chain, radius, positions, basis)
        return score

    def optimize_versus(self, other, level=0, max_level=8, trials=100, coverage=1, radius=12, basis=None):
        if not basis:
            basis = self.basis
        for i in xrange(0, coverage):
            best_basis = basis
            best_score = self.score_chain(other, radius=radius, basis=basis)
            for i in xrange(0, trials):
                new_basis = self.random_perturbation(level, basis=basis)
                score = self.score_chain(other, radius=radius, basis=new_basis)
                if score < best_score:
                    best_basis = new_basis
                    best_score = score
            basis = best_basis
        if level < max_level:
            self.optimize_versus(other, level=level+1, max_level=max_level, trials=trials, coverage=coverage, radius=radius, basis=basis)
        else:
            self.set_basis(basis)

    def random_perturbation(self, level=0, basis=None):
        if not basis:
            basis = self.basis
        # 20% of perturbations are translations, 80% are rotations
        if random.random() < 0.2:
            return self.random_translation(level, basis=basis)
        else:
            return self.random_rotation(level, basis=basis)

    def random_translation(self, level=0, basis=None):
        if not basis:
            basis = self.basis
        d = 0.5 ** level
        shift = self.random_from_ball() * d
        return basis.translate(shift)
        
    def random_rotation(self, level=0, basis=None):
        if not basis:
            basis = self.basis
        phi, theta = self.random_from_cap(0.5 ** level)
        return basis.rotate_x(phi).rotate_y(theta)
        
    # uniform distributions of points within the unit circle        
    def random_from_disc_polar(self):
        r = random.random()
        t = random.random()
        if t > r:
            t, r = r, t
        if r == 0:
            return 0.0, 0.0
        else:
            return r, 2 * pi * t / r

    def random_from_disc_cartesian(self):
        r, theta = self.random_from_disc_polar()
        if r == 0:
            return 0.0, 0.0
        else:
            x = r * cos(theta)
            y = r * sin(theta)
            return x, y

    # *non-uniform* distribution of spherical coordinates on the unit sphere
    # with phi <= t. smaller t => closer to uniform; t >= PI/2 just plain
    # doesn't work
    def random_from_cap(self, t):
        phi, theta = self.random_from_disc_polar()
        phi *= t
        return phi, theta
        
    # uniform distribution of points within the unit ball
    def random_from_ball(self):
        x1, y1 = self.random_from_disc_cartesian()
        x2, y2 = self.random_from_disc_cartesian()
        return numpy.array([x1, y1 * x2, y1 * y2])
