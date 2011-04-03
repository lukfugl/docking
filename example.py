from xbgf import XBGFParser

parser = XBGFParser()
chains = parser.parse('2PWS.xbgf', '2PWS.xbgf')
chain = chains[0]
print chain.score_chain(chain)
chain.rotate_y(1)
print chain.score_chain(chain)
