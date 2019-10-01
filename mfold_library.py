import os
import string
import subprocess


class Strand:
    allowed_bases = set('ATCG')
    allowed_constraints = set(string.ascii_letters + string.digits)
    

    def __init__(self, bases, constraints):
        self.bases = bases.upper()
        self.constraints = constraints
        if set(self.bases) > Strand.allowed_bases:
            raise TypeError('The selected bases contain letters '
                          + 'other than A, T, C, and G: ' + bases)
        if set(self.constraints) > Strand.allowed_constraints:
            raise TypeError('The selected constraints contain '
                          + 'non-alphanumeric characters: ' + constraints)



class Region:
    def __init__(self, char, index, length):
        # the character that represents the region, e.g. 'A3' -> 'A'
        self.char = char
        # the index that the region starts at in the sequence
        self.start = index
        # one past the index that the region ends at in the sequence
        self.end = index + length



class Mfold:
    energy_string = ' Initial dG = '
    linker_sequence = 'LLL'
    output_suffixes = ['.aux', '.cmd', '.con', '.log', '.pnt', '.sav', '.seq'
            '-local.pnt', '-local.seq', '.ann', '.ct', '.ps', '.det', '.out',
            '.h-num', '.plot', '.pdf', '.ss-count', '-temp.det', '-temp.out']


    def __init__(self, output_folder='', mfold_command=''):
        self.folder = output_folder
        self.command = mfold_command


    def run(self, strand1, strand2, sequence_file='a.seq', settings_file='a.aux'):
        seq_path = os.path.join(self.folder, sequence_file)
        set_path = os.path.join(self.folder, settings_file)

        with open(seq_path, 'w') as seqfile:
            seqfile.write(strand1.bases + Mfold.linker_sequence + strand2.bases)
        with open(set_path, 'w') as setfile:
            for constraint in Mfold.get_constraints(strand1, strand2):
                setfile.write(constraint)

        subprocess.run([self.command, f'SEQ={seq_path}', f'AUX={set_path}'],
                cwd=self.folder)


    def clean(self, file_prefix):
        for suffix in Mfold.output_suffixes:
            file_path = os.path.join(self.folder, f'{file_prefix}{suffix}')
            if os.path.exists(file_path):
                os.remove(file_path)


    def get_energy(self, details_file='a.det'):
        details_path = os.path.join(self.folder, details_file)
        if os.path.exists(details_path):
            with open(details_path, 'r') as detfile:
                for line in detfile:
                    if line.startswith(Mfold.energy_string):
                        return float(line[len(Mfold.energy_string):])
        return None
            

    def get_constraints(strand1, strand2):
        regions1 = Mfold.index_constraints(strand1, 0)
        regions2 = Mfold.index_constraints(
                strand2, len(strand1.constraints) + 3)
        constraints = []

        for lower_char in string.ascii_lowercase:
            upper_char = lower_char.upper()
            if lower_char in regions1:
                lower_segment = regions1[lower_char]
            elif lower_char in regions2:
                lower_segment = regions2[lower_char]
            else:
                lower_segment = None

            if upper_char in regions1:
                upper_segment = regions1[upper_char]
            elif upper_char in regions2:
                upper_segment = regions2[upper_char]
            else:
                upper_segment = None

            if lower_segment and upper_segment:
                constraints.append(
                        f'P {lower_segment.start} {lower_segment.end} '
                        + f'{upper_segment.start} {upper_segment.end}')
        return constraints


    # TODO: can only handle 26 or fewer pairing regions
    def index_constraints(strand, starting_index=0):
        # index is the index of the constraint in the final strand
        index = starting_index
        regions = dict()

        for i in range(len(strand.constraints)):
            char = strand.constraints[i]
            if char.isalpha():
                j = i+1
                while (j < len(strand.constraints) 
                        and strand.constraints[j].isdigit()):
                    j += 1
                length = int(strand.constraints[i+1:j])
                regions[char] = Region(char, index, length)
                index += length

        assert(index - starting_index == len(strand.bases))
        return regions



class EnergyMatrix:
    def __init__(self, mfold, strands):
        self.mfold = mfold
        self.strands = strands
        self.matrix = [[None for strand1 in strands] for strand2 in strands]

    def create(self):
        for i, strand1 in enumerate(self.strands):
            for j, strand2 in enumerate(self.strands):
                self.mfold.clean(f'{i}_{j}')
                self.mfold.run(strand1, strand2, f'{i}_{j}.seq', f'{i}_{j}.aux')
        self.matrix = [[self.mfold.get_energy(f'{i}_{j}.det')
                for i in range(len(self.strands))] 
                    for j in range(len(self.strands))]
