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

#GeneticAlgorithm([[Region('b', 4), Region('c', 4)], [Region('C', 4)]], population_size=3, iterations=4).run()
GeneticAlgorithm(
	[[Region('a', 25), Region('B', 25)], [Region('b', 25), Region('C', 25)], [Region('c', 25), Region('D', 25)], [Region('d', 25), Region('A', 25)]],
	population_size=20,
	iterations=20,
	mutation_rate=500
).run()

