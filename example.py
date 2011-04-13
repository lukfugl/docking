from xbgf import XBGFParser
from pymol import cmd

parser = XBGFParser()
chains = parser.parse('2PWS.xbgf', '2PWS.xbgf')
protein = chains[0]
chains = parser.parse('ligand_out.xbgf', 'ligand_out.xbgf')
ligand = chains[0]

protein.push2pymol('protein')
ligand.push2pymol('ligand-before')
print ligand.score_chain(protein)
ligand.optimize_versus(protein)
ligand.push2pymol('ligand-after1')
print ligand.score_chain(protein)
ligand.reset_basis()
ligand.optimize_versus(protein)
ligand.push2pymol('ligand-after2')
print ligand.score_chain(protein)
ligand.reset_basis()
ligand.optimize_versus(protein)
ligand.push2pymol('ligand-after3')
print ligand.score_chain(protein)

cmd.zoom('protein')
cmd.mplay()
