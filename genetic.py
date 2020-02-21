from mfold_library import Strand, Region, Mfold, EnergyMatrix
from math import exp, sqrt
import numpy as np
from random import randrange, random, sample, choice
import json

TEMPERATURE = 310.15
BOLTZMANN = 1.38064852e-23
AVOGADRO = 6.0221409e23

class Sequence:
	"""
	A class used to represent a Sequence (which is made up of several Strands). It is defined by a strand
	structure and a dictionary defining the bases of the Regions with lowercased names.
	"""
	def __init__(self, region_definitions, strand_structures):
		"""
		Args:
			region_definitions: A map from lowercased region names to string of bases
			strand_structures: A list of strand structures. Each strand structure is represented by a list of Region.
		"""
		self.region_definitions = region_definitions
		self.strand_structures = strand_structures

	@staticmethod
	def random_sequence(sequence_structure):
		"""
		Generates a random Sequence object with a given structure.
		Args:
			sequence_structure: A list of Regions for each strand.
		Returns:
			A Sequence object with randomized bases and the given structure.
		"""
		region_defs = {}
		for strand_structure in sequence_structure:
			for region in strand_structure:
				if not region.name.lower() in region_defs:
					region_defs[region.name.lower()] = "".join([choice(list(Strand.allowed_bases))
														for i in range(0, region.length)])
		return Sequence(region_defs, sequence_structure)

	def mutate(self, mutation_rate):
		"""
		Mutates the sequence.
		Args:
			mutation_rate: Mutate 1 out of every mutation_rate bases
		"""
		for region in self.region_definitions:
			bases = list(self.region_definitions[region])
			for i in range(len(bases)):
				if randrange(mutation_rate) == 0:
					bases[i] = sample(Strand.allowed_bases, 1)[0]
			self.region_definitions[region] = "".join(bases)

	@staticmethod
	def _mate_bases(bases1, bases2):
		"""
		Creates a new string of bases by mating 2 bases together.
		Args:
			bases1: A string of bases representing the first parent.
			bases2: A string of bases representing the second parent.
		Returns:
			A string of bases, with each base randomly chosen from one of the two parent strings.
		"""
		child = list(bases1)
		for i in range(len(bases1)):
			if randrange(1) == 0:
				child[i] = bases2[i]
		return "".join(child)

	@staticmethod
	def _mate_bases_crossover(bases1, bases2):
		crosspoint = np.random.randint(0, high=len(bases1))
		if random() > 0.5:
			return bases1[:crosspoint] + bases2[crosspoint:]
		else:
			return bases2[:crosspoint] + bases1[crosspoint:]

	@staticmethod
	def mate(sequence1, sequence2):
		"""
		Create a new Sequence object by mating 2 sequences together.
		Args:
			sequence1: The first parent to mate.
			sequence2: The second parent to mate.
		Returns:
			A new Sequence object, created by mating each Region definition.
		"""
		if sequence1.strand_structures != sequence2.strand_structures:
			raise ValueError('The sequences being mated have different structures')

		child_regions = {}
		for region in sequence1.region_definitions:
			child_regions[region] = Sequence._mate_bases_crossover(sequence1.region_definitions[region], sequence2.region_definitions[region])

		return Sequence(child_regions, sequence1.strand_structures)

	def build_strand(self, strand_structure):
		"""
		Given a strand structure, generate the Strand object.
		Args:
			strand_structure: A list of Regions
		Returns:
			A Strand object with bases from the region_definitions.
		"""
		bases = ""
		for region in strand_structure:
			if region.name.islower():
				bases += self.region_definitions[region.name]
			else:
				bases += Strand.complement(self.region_definitions[region.name.lower()])[::-1] # Reversed, because strands only bind in the opposite direction

		return Strand(bases, strand_structure)

	def fitness(self, mfold, cache):
		"""
		Calculate the fitness of the sequence.
		Args:
			mfold: The mfold object to run calculations with.
		Returns:
			A number representing the fitness of the sequence.
		"""
		region_hash = json.dumps(self.region_definitions, sort_keys=True)
		if not region_hash in cache:
			strands = [self.build_strand(strand_structure) for strand_structure in self.strand_structures]
			base_content = [strand.base_content() for strand in strands]
			at = sum([bc[0] for bc in base_content])
			gc = sum([bc[1] for bc in base_content])
			maxrun = max([bc[2] for bc in base_content])
			x = at/(at + gc)
			penalty = ((8.0/13 * x + 1.0)/(4.0/13 + 1.0))**4
			if maxrun > 3:
				penalty *= maxrun / 3
			penalty /=  len(strands)
			energy_matrix = EnergyMatrix(mfold, strands, penalty)
			energy_matrix.create()
			cache[region_hash] = energy_matrix.matrix
		return np.linalg.norm(cache[region_hash], ord=1)

	def print(self):
		"""
		Prints out the strands in the sequence.
		"""
		print("SEQUENCE:")
		for strand_struct in self.strand_structures:
			built_strand = self.build_strand(strand_struct)
			print(built_strand.bases)

class GeneticAlgorithm:
	"""
	Implementation of the genetic algorithm
	"""
	def __init__(self, structure, mfold_command, population_size=50, mutation_rate=100, iterations=100, boltzmann_factor=1, initial_sequences=[]):
		"""
		Args:
			structure: A list of strand structures
			population_size: The number of sequences in a population
			mutation_rate: Reciprocal of the rate of mutation
			initial_sequences: A list of user defined sequences to include in the initial population
		Attributes:
			population: A list of sequences
		"""
		self.iterations = iterations
		self.population_size = population_size
		self.mutation_rate = mutation_rate
		self.boltzmann_factor = boltzmann_factor * 1000/(TEMPERATURE * AVOGADRO * BOLTZMANN)
		self.population = initial_sequences + [Sequence.random_sequence(structure) for i in range(population_size - len(initial_sequences))]
		self.mfold = Mfold(output_folder='./', mfold_command=mfold_command)
		self.cache = {}
		self.fitness_history = []
		self.diversity_history = []
		self.best_child = None
		
	def iterate(self):
		"""
		Do one iteration of the genetic algorithm.
		"""
		# Find the fitness of each sequence in the population
		fitnesses = [sequence.fitness(self.mfold, self.cache) for sequence in self.population]
		power_sum = sum([exp(-fitness*self.boltzmann_factor) for fitness in fitnesses])
		weighted_fitnesses = [exp(-fitness*self.boltzmann_factor)/power_sum for fitness in fitnesses]
		self.fitness_history += [fitnesses]

		# Save the best child
		best_child = self.population[np.argmax(weighted_fitnesses)]
		self.best_child = best_child

		# Mate strands at random, weighted by fitness level
		#midpoint = sum(weighted_fitnesses[:int(self.population_size/2)])
		self.population = [self.generate_child(weighted_fitnesses) for i in range(self.population_size - 1)]

		# Mutate the strands
		for sequence in self.population:
			sequence.mutate(self.mutation_rate)

		self.population.append(best_child)

	def _round_up(self, weights, number):
		"""
		Given a list of weights, find out which member of the population a number refers to.
		Example: ([0.5, 0.5], 0.75) => 1
				 ([0.5, 0.5], 0.25) => 0
		Args:
			weights: The probability of each member of the population to be chosen as a parent.
			number: The number that we want to find the Sequence of.
		Returns:
			The Sequence in the population.
		"""
		curr = 0
		for i in range(len(weights)):
			curr += weights[i]
			if curr > number:
				return self.population[i]

		return self.population[-1]

	def generate_child(self, weighted_fitnesses):
		"""
		Generate a child from the population.
		Args:
			weighted_fitnesses: The probability of each member of the population to be chosen as a parent.
		Returns:
			A child Sequence.
		"""
		parent1 = self._round_up(weighted_fitnesses, random())
		parent2 = self._round_up(weighted_fitnesses, random())
		return Sequence.mate(parent1, parent2)

	def generate_child_segregated(self, weighted_fitnesses, divide):
		"""
		Generate a child from the population where the population is segregated at the divide point.
		"""
		if random() > 0.5:
			# Pick from below divide
			parent1 = self._round_up(weighted_fitnesses, divide * random())
			parent2 = self._round_up(weighted_fitnesses, divide * random())
			return Sequence.mate(parent1, parent2)
		else:
			parent1 = self._round_up(weighted_fitnesses, 1 - (1 - divide) * random())
			parent2 = self._round_up(weighted_fitnesses, 1 - (1 - divide) * random())
			return Sequence.mate(parent1, parent2)

	def run(self):
		"""
		Run the genetic algorithm.
		"""
		for i in range(self.iterations):
			print("ITERATION", i)
			self.diversity_history.append(self.diversity())
			self.iterate()

	def diversity(self):
		"""
		Computes the diversity of the Sequences in the population.
		"""
		first_region_label = self.population[0].strand_structures[0][0].name.lower()
		region_length = len(self.population[0].region_definitions[first_region_label])
		base_counts = [{'A': 0, 'C': 0 , 'T': 0, 'G': 0} for x in range(region_length)]
		for seq in self.population:
			bases = seq.region_definitions[first_region_label]
			for i in range(region_length):
				base_counts[i][bases[i]] += 1
		avg = self.population_size/4.0
		base_diversity = [sqrt((base['A'] - avg)**2 + (base['C'] - avg)**2 + (base['T'] - avg)**2 + (base['G'] - avg)**2)/sqrt(3) for base in base_counts]
		return sum(base_diversity)/len(base_diversity)

	def print_population(self):
		"""
		Prints all the sequences in the population.
		"""
		for sequence in self.population:
			sequence.print()
