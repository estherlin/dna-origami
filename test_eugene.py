"""
from mfold_library import Strand, Mfold, EnergyMatrix

strands = [Strand('AAAAA', 'a1b3c1'), Strand('TGGGT', 'A1D3C1')]
mfold = Mfold('/home/h/henryx/eugene_enph479/mfold_test', '/home/h/henryx/.local/bin/mfold')
energy_matrix = EnergyMatrix(mfold, strands)
energy_matrix.create()
print(energy_matrix.matrix)
"""

from genetic import GeneticAlgorithm
from mfold_library import Region

GeneticAlgorithm([[Region('b', 4), Region('c', 4)], [Region('C', 4)]], population_size=3, iterations=4).run()
