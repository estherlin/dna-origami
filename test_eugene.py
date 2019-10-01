from mfold_library import Strand, Mfold, EnergyMatrix

strands = [Strand('AAAAA', 'a1b3c1'), Strand('TGGGT', 'A1D3C1')]
mfold = Mfold('/home/h/henryx/eugene_enph479/mfold_test', '/home/h/henryx/.local/bin/mfold')
energy_matrix = EnergyMatrix(mfold, strands)
energy_matrix.create()
print(energy_matrix.matrix)