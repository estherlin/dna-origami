import numpy as np
from mfold_library import *

#bases = ['AGGCACAGCTATAATAACGCAATCCTCTCCGGCCTCAAACTACTTTACCT', 'AGGTAAAGTAGTTTGAGGCCGGAGACCGAATGGAGTCTGTTCTCGACGCT', 'AGCGTCGAGAACAGACTCCATTCGGACAATTACGAACCAACTTAGGACCT', 'AGGTCCTAAGTTGGTTCGTAATTGTGGATTGCGTTATTATAGCTGTGCCT']
bases = ['AGGAGCCTATCGGGTAGATCGAAGACGTACAGGTGTGACTTGAATTTGCT', 'AGCAAATTCAAGTCACACCTGTACGAGTGTTAGAATACAACAAGCGACCT', 'AGGTCGCTTGTTGTATTCTAACACTGCATCTCATACGGCAGTATCCGCCT', 'AGGCGGATACTGCCGTATGAGATGCTCTTCGATCTACCCGATAGGCTCCT']
regions = [[Region('A', 25), Region('b', 25)], [Region('B', 25), Region('c', 25)], [Region('C', 25), Region('d', 25)], [Region('D', 25), Region('a', 25)]]
mfold=Mfold(output_folder='./', mfold_command='/home/ubuntu/dna-origami/mfold_quik')
em=EnergyMatrix(mfold, [Strand(bases[i], regions[i]) for i in range(4)])
em.create()
print(em.matrix)
print(np.linalg.norm(em.matrix))
