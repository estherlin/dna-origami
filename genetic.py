from mfold_library import Strand, Region, Mfold, EnergyMatrix
from math import exp
import numpy as np
from random import randrange, random, sample, choice
import json

class Sequence:
	"""
	Args:
		region_definitions: A map from lowercased region names to string of bases
		strand_structures: A list of strand structures. Each strand structure is represented by a list of Region.
	"""
	def __init__(self, region_definitions, strand_structures):
		self.region_definitions = region_definitions
		self.strand_structures = strand_structures

	"""
	Generate a random Sequence object with a given structure.
	Args:
		sequence_structure: A list of Regions for each strand.
	Returns:
		A Sequence object with randomized bases and the given structure.
	"""
	@staticmethod
	def random_sequence(sequence_structure):
		region_defs = {}
		for strand_structure in sequence_structure:
			for region in strand_structure:
				if not region.name.lower() in region_defs:
					region_defs[region.name] = "".join([choice(list(Strand.allowed_bases))
														for i in range(0, region.length)])
		return Sequence(region_defs, sequence_structure)

	"""
	Mutate the sequence.
	Args:
		mutation_rate: Mutate 1 out of every mutation_rate bases
	"""
	def mutate(self, mutation_rate):
		for region in self.region_definitions:
			bases = list(self.region_definitions[region])
			for i in range(len(bases)):
				if randrange(mutation_rate) == 0:
					bases[i] = sample(Strand.allowed_bases, 1)[0]
			self.region_definitions[region] = "".join(bases)

	"""
	Create a new string of bases by mating 2 bases together.
	Args:
		bases1: A string of bases representing the first parent.
		bases2: A string of bases representing the second parent.
	Returns:
		A string of bases, with each base randomly chosen from one of the two parent strings.
	"""
	@staticmethod
	def _mate_bases(bases1, bases2):
		child = list(bases1)
		for i in range(len(bases1)):
			if randrange(1) == 0:
				child[i] = bases2[i]
		return "".join(child)

	"""
	Create a new Sequence object by mating 2 sequences together.
	Args:
		sequence1: The first parent to mate.
		sequence2: The second parent to mate.
	Returns:
		A new Sequence object, created by mating each Region definition.
	"""
	@staticmethod
	def mate(sequence1, sequence2):
		if sequence1.strand_structures != sequence2.strand_structures:
			raise ValueError('The sequences being mated have different structures')

		child_regions = {}
		for region in sequence1.region_definitions:
			child_regions[region] = Sequence._mate_bases(sequence1.region_definitions[region], sequence2.region_definitions[region])

		return Sequence(child_regions, sequence1.strand_structures)

	"""
	Given a strand structure, generate the Strand object.
	Args:
		strand_structure: A list of Regions
	Returns:
		A Strand object with bases from the region_definitions.
	"""
	def build_strand(self, strand_structure):
		bases = ""
		for region in strand_structure:
			if region.name.islower():
				bases += self.region_definitions[region.name]
			else:
				bases += Strand.complement(self.region_definitions[region.name.lower()])[::-1] # Reversed, because strands only bind in the opposite direction

		return Strand(bases, strand_structure)

	"""
	Calculate the fitness of the sequence.
	Args:
		mfold: The mfold object to run calculations with.
	Returns:
		A number representing the fitness of the sequence.
	"""
	def fitness(self, mfold, cache):
		region_hash = json.dumps(self.region_definitions, sort_keys=True)
		if not region_hash in cache:
			strands = [self.build_strand(strand_structure) for strand_structure in self.strand_structures]
			energy_matrix = EnergyMatrix(mfold, strands)
			energy_matrix.create()
			cache[region_hash] = energy_matrix.matrix
		return exp(-np.linalg.norm(cache[region_hash]))

	"""
	Prints out the strands in the sequence.
	"""
	def print(self):
		print("SEQUENCE:")
		for strand_struct in self.strand_structures:
			built_strand = self.build_strand(strand_struct)
			print(built_strand.bases)

class GeneticAlgorithm:
	"""
	Args:
		structure: A list of strand structures
		population_size: The number of sequences in a population
		mutation_rate: Reciprocal of the rate of mutation
		initial_sequences: A list of user defined sequences to include in the initial population
	Attributes:
		population: A list of sequences
	"""
	def __init__(self, structure, population_size=50, mutation_rate=100, iterations=100, initial_sequences=[]):
		self.iterations = iterations
		self.population_size = population_size
		self.mutation_rate = mutation_rate
		self.population = initial_sequences + [Sequence.random_sequence(structure) for i in range(population_size - len(initial_sequences))]
		self.mfold = Mfold(output_folder='./', mfold_command='mfold')
		self.cache = {}

	"""
	Do one iteration of the genetic algorithm.
	"""
	def iterate(self):
		# Find the fitness of each sequence in the population
		fitnesses = [sequence.fitness(self.mfold, self.cache) for sequence in self.population]
		weighted_fitnesses = [fitness/sum(fitnesses) for fitness in fitnesses]

		# Mate strands at random, weighted by fitness level
		self.population = [self.generate_child(weighted_fitnesses) for i in range(self.population_size)]

		# Mutate the strands
		for sequence in self.population:
			sequence.mutate(self.mutation_rate)

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
	def _round_up(self, weights, number):
		curr = 0
		for i in range(len(weights)):
			curr += weights[i]
			if curr > number:
				return self.population[i]

		return self.population[-1]

	"""
	Generate a child from the population.
	Args:
		weighted_fitnesses: The probability of each member of the population to be chosen as a parent.
	Returns:
		A child Sequence.
	"""
	def generate_child(self, weighted_fitnesses):
		parent1 = self._round_up(weighted_fitnesses, random())
		parent2 = self._round_up(weighted_fitnesses, random())
		return Sequence.mate(parent1, parent2)

	"""
	Run the genetic algorithm.
	"""
	def run(self):
		#self.mfold.clean_all()
		for i in range(self.iterations):
			self.iterate()
		for sequence in self.population:
			sequence.print()
			seq_hash = json.dumps(sequence.region_definitions, sort_keys=True)
			if seq_hash in self.cache:
				print("INTERACIONS:")
				print(self.cache[seq_hash])
