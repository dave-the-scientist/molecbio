"""
Simple objects for comparing protein alignment scores.

    These are implemented as dictionaries, allowing to residues to be compared
by calling B['A','C'], which would return the value. The objects also implement
an 'alphabet' attribute, listing all letters they contain, as well as the
'min' and 'max' values. This module currently defines the following objects:
    blosum30
    blosum45
    blosum50
    blosum62
    blosum80
    matrices -- A list containing the names of the above.
"""
import os, math
import matrices


class Blosum(dict):
    def __init__(self, matrixStr):
        dict.__init__(self)
        self.alphabet = []
        self.min = 0
        self.max = 0
        self.lambda_value = 0.0
        self.entropy = 0.0
        self.expected = 0.0
        self.__initFromNcbiStr(matrixStr)

    def __initFromNcbiStr(self, matStr):
        minVal, maxVal = 100, None
        it = iter(matStr.splitlines())
        for line in it:
            if line.strip().lower().endswith('bit units'):
                bitFraction = line.rpartition('/')[2].split()[0]
                self.lambda_value = math.log(2, math.e)/float(bitFraction)
                continue
            if 'entropy' in line.lower() and 'expected' in line.lower():
                ent, _, exp = line.partition(',')
                try:
                    self.entropy = float(ent[ent.index('=')+1:])
                    self.expected = float(exp[exp.index('=')+1:])
                except ValueError: pass
                continue
            if line.startswith('#') or line == '\n': continue
            headerLine = line.upper().split()
            self.alphabet = headerLine
            break
        for line in it:
            if not line.strip(): continue
            line = line.split()
            c = line[0].upper()
            c_low = c.lower()
            for c2, num in zip(headerLine, line[1:]):
                c2_low = c2.lower()
                num = int(num)
                if num < minVal: minVal = num
                if num > maxVal: maxVal = num
                self[c, c2] = num
                self[c2, c] = num
                self[c_low, c2_low] = num
                self[c2_low, c_low] = num
        self.min = minVal
        self.max = maxVal

# # # # #  Constants  # # # # #
blosum30 = Blosum(matrices.BLOSUM30)
blosum45 = Blosum(matrices.BLOSUM45)
blosum50 = Blosum(matrices.BLOSUM50)
blosum62 = Blosum(matrices.BLOSUM62)
blosum80 = Blosum(matrices.BLOSUM80)

matrices = ['blosum30', 'blosum45', 'blosum50', 'blosum62', 'blosum80']
