import warnings

import numpy

from Bio.PDB.PDBExceptions import \
        PDBConstructionException, PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import Selection

from chain import Chain
  
ELEMENTS = {  1: 'H',     2: 'HE',    3: 'LI',    4: 'BE',    5: 'B',
              6: 'C',     7: 'N',     8: 'O',     9: 'F',    10: 'NE',
             11: 'NA',   12: 'MG',   13: 'AL',   14: 'SI',   15: 'P',
             16: 'S',    17: 'CL',   18: 'AR',   19: 'K',    20: 'CA',
             21: 'SC',   22: 'TI',   23: 'V',    24: 'CR',   25: 'MN',
             26: 'FE',   27: 'CO',   28: 'NI',   29: 'CU',   30: 'ZN',
             31: 'GA',   32: 'GE',   33: 'AS',   34: 'SE',   35: 'BR',
             36: 'KR',   37: 'RB',   38: 'SR',   39: 'Y',    40: 'ZR',
             41: 'NB',   42: 'MO',   43: 'TC',   44: 'RU',   45: 'RH',
             46: 'PD',   47: 'AG',   48: 'CD',   49: 'IN',   50: 'SN',
             51: 'SB',   52: 'TE',   53: 'I',    54: 'XE',   55: 'CS',
             56: 'BA',   57: 'LA',   58: 'CE',   59: 'PR',   60: 'ND',
             61: 'PM',   62: 'SM',   63: 'EU',   64: 'GD',   65: 'TB',
             66: 'DY',   67: 'HO',   68: 'ER',   69: 'TM',   70: 'YB',
             71: 'LU',   72: 'HF',   73: 'TA',   74: 'W',    75: 'RE',
             76: 'OS',   77: 'IR',   78: 'PT',   79: 'AU',   80: 'HG',
             81: 'TL',   82: 'PB',   83: 'BI',   84: 'PO',   85: 'AT',
             86: 'RN',   87: 'FR',   88: 'RA',   89: 'AC',   90: 'TH',
             91: 'PA',   92: 'U',    93: 'NP',   94: 'PU',   95: 'AM',
             96: 'CM',   97: 'BK',   98: 'CF',   99: 'ES',  100: 'FM',
            101: 'MD',  102: 'NO',  103: 'LR',  104: 'RF',  105: 'DB',
            106: 'SG',  107: 'BH',  108: 'HS',  109: 'MT',  110: 'DS',
            111: 'RG',  112: 'CN',  113: 'UUT', 114: 'UUQ', 115: 'UUP',
            116: 'UUH', 117: 'UUS', 118: 'UUO' }

class XBGFParser:
    def __init__(self, PERMISSIVE=1, structure_builder=None):
        if structure_builder != None:
            self.structure_builder = structure_builder
        else:
            self.structure_builder = StructureBuilder()
        self.PERMISSIVE = PERMISSIVE

    # public interface

    def parse(self, id, file):
        self.structure_builder.init_structure(id)
        if isinstance(file, basestring):
            file=open(file)
        self.charges = dict()
        self._parse(file.readlines())
        self.structure = self.structure_builder.get_structure()
        return self._process_structure()

    # private methods

    def _parse(self, lines):
        self.structure_builder.init_model(0)
        self.structure_builder.init_seg("")
        self.current_chain_id = None
        self.current_residue_id = None
        self.current_resname = None
        for i in range(0, len(lines)):
            self.line_counter = i + 1
            self.structure_builder.set_line_counter(self.line_counter)
            line = lines[i]
            if line[0:6] == 'ATOM  ':
                self._parse_atom(line)

    def _parse_atom(self, line):
        self._update_chain(line)
        self._update_residue(line)
        self._update_atom(line)

    def _update_chain(self, line):
        chain_id = self._extract_chain(line)
        if self.current_chain_id != chain_id:
            try:
                self.structure_builder.init_chain(chain_id)
                self.current_chain_id = chain_id
                self.current_residue_id = None
                self.current_resname = None
            except PDBConstructionException, message:
                self._handle_PDB_exception(message) 

    def _update_residue(self, line):
        resseq, resname, hetero_flag = self._extract_residue(line)
        residue_id = (hetero_flag, resseq, " ")
        if self.current_residue_id != residue_id or self.current_resname != resname:
            try:
                self.structure_builder.init_residue(resname, hetero_flag, resseq, " ")
                self.current_residue_id = residue_id
                self.current_resname = resname
            except PDBConstructionException, message:
                self._handle_PDB_exception(message) 

    def _update_atom(self, line):
        fullname, name = self._extract_name(line)
        serial_number = self._extract_serial_number(line)
        coord = self._extract_coord(line)
        bfactor = self._extract_bfactor(line)
        occupancy = self._extract_occupancy(line)
        element = self._extract_element(line)
        try:
            self.structure_builder.init_atom(name, coord, bfactor, occupancy, " ",
                                        fullname, serial_number, element)
        except PDBConstructionException, message:
            self._handle_PDB_exception(message)
            return
        self.charges[self.structure_builder.atom] = self._extract_charge(line)

    def _extract_chain(self, line):
        return line[25:26]

    def _extract_residue(self, line):
        resseq = int(line[27:32].split()[0])
        resname = line[20:24]
        if is_aa(resname):
            hetero_flag = " "
        elif resname == "HOH" or resname == "WAT":
            hetero_flag = "W"
        else:
            hetero_flag = "H"
        return resseq, resname, hetero_flag

    def _extract_name(self, line):
        fullname = line[14:19]
        split_list = fullname.split()
        if len(split_list) != 1:
            # atom name has internal spaces, e.g. " N B ", so we do not strip
            # spaces
            name = fullname
        else:
            # atom name is like " CA ", so we can strip spaces
            name = split_list[0]
        return fullname, name

    def _extract_serial_number(self, line):
        try:
            return int(line[10:13])
        except:
            return 0

    def _extract_coord(self, line):
        try:
            x = float(line[32:42]) 
            y = float(line[42:52]) 
            z = float(line[52:62])
            return numpy.array((x, y, z), 'f')
        except:
            self._handle_PDB_exception("Invalid or missing coordinate(s)")
            return numpy.array((0, 0, 0), 'f')

    def _extract_bfactor(self, line):
        try:
            return float(line[83:90])
        except:
            self._handle_PDB_exception("Invalid or missing bfactor")
            return 0

    def _extract_occupancy(self, line):
        try:
            return float(line[90:97])
        except:
            self._handle_PDB_exception("Invalid or missing occupancy")
            return 0

    def _extract_element(self, line):
        try:
            return ELEMENTS[int(line[97:101])]
        except:
            self._handle_PDB_exception("Invalid or missing element")
            return ""

    def _extract_charge(self, line):
        try:
            return float(line[74:82])
        except:
            self._handle_PDB_exception("Invalid or missing charge")
            return 0

    def _handle_PDB_exception(self, message):
        """
        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE==1), or raises it again, this time adding the 
        PDB line number to the error message.
        """
        message="%s at line %i." % (message, self.line_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn("PDBConstructionException: %s\n"
                          "Exception ignored.\n"
                          "Some atoms or residues may be missing in the data structure."
                          % message, PDBConstructionWarning)
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)

    def _process_structure(self):
        return [Chain(chain, self.charges) for chain in
                Selection.unfold_entities(self.structure, 'C')]
