from xbgf import XBGFParser

parser = XBGFParser()
chains = parser.parse('2PWS.xbgf', '2PWS.xbgf')
protein = chains[0]
ligand = chains[1]
print protein.score_chain(ligand)
ligand.rotate_y(1)
print protein.score_chain(ligand)
