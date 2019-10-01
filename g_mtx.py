import os
import string

def get_restrictions(struct1, struct2):
        regions1 = {}
        regions2 = {}
        index = 0
        for i in range(0, len(struct1), 2):
                regions1[struct1[i]] = (index, ord(struct1[i+1]) - ord('0'))
                index += ord(struct1[i+1]) - ord('0')
        index = len(struct1) + 3
        for i in range(0, len(struct2), 2):
                regions2[struct2[i]] = (index, ord(struct2[i+1]) - ord('0'))
                index += ord(struct2[i+1]) - ord('0')

        restrictions = []

        # a very very bad assumption: there are 26 or fewer pairing regions
        for i in range(26):
                lower_seg = False
                upper_seg = False
                if string.ascii_lowercase[i] in regions1:
                        lower_seg = regions1[string.ascii_lowercase[i]]
                if string.ascii_lowercase[i] in regions2:
                        lower_seg = regions2[string.ascii_lowercase[i]]
                if string.ascii_uppercase[i] in regions1:
                        upper_seg = regions1[string.ascii_uppercase[i]]
                if string.ascii_uppercase[i] in regions2:
                        upper_seg = regions2[string.ascii_uppercase[i]]

                if lower_seg and upper_seg:
                        restrictions.append('P {0} {1} {2} {3}'.format(lower_seg[0], lower_seg[0] + lower_seg[1], upper_seg[0], upper_seg[0] + upper_seg[1]))
        return restrictions

def run_mfold(strand1, strand2):
        seq_file = open('a.seq', 'w+')
        seq_file.write(strand1[0] + 'LLL' + strand2[0])
        seq_file.close()
        settings_file = open('a.aux', 'w+')
        print(strand1, strand2)
        for restriction in get_restrictions(strand1[1], strand2[1]):
                print(restriction)
                settings_file.write(restriction)
        settings_file.close()
        os.system('mfold SEQ=a.seq AUX=a.aux')
        # get the energy output

def make_energy_matrix(strands):
        matrix = [[run_mfold(i, j) for j in strands] for i in strands]

strands = [('AAAAA', 'a1b3c1'),('TGGGT', 'A1D3C1')]
print(make_energy_matrix(strands))
