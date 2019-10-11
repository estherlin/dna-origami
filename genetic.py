from mfold_library import Strand, Region, Mfold, EnergyMatrix
from math import exp
import numpy as np
from random import randrange, random, sample

class Sequence:
	"""
	Args:
		region_definitions: A map from lowercased region names to string of bases
		strand_structures: A list of strand structures. Each strand structure is represented by a list of Region.
	"""
	def __init__(self, region_definitions, strand_structures):
		self.region_definitions = region_definitions
		self.strand_structures = strand_structures

	def random_sequence(sequence_structure):
		region_defs = {}
		for strand_structure in sequence_structure:
			for region in strand_structure:
				if not region.name.lower() in region_defs:
					region_defs[region.name] = "".join(sample(Strand.allowed_bases, region.length))
		return Sequence(region_defs, sequence_structure)

	def mutate(self, mutation_rate):
		for region in self.region_definitions:
			bases = list(self.region_definitions[region])
			for i in range(len(bases)):
				if randrange(mutation_rate) == 0:
					bases[i] = sample(Strand.allowed_bases, 1)[0]
			self.region_definitions[region] = "".join(bases)

	def _mate_bases(bases1, bases2):
		child = list(bases1)
		for i in range(len(bases1)):
			if randrange(1) == 0:
				child[i] = bases2[i]
		return "".join(child)

	def mate(sequence1, sequence2):
		if sequence1.strand_structures != sequence2.strand_structures:
			raise ValueError('The sequences being mated have different structures')

		child_regions = {}
		for region in sequence1.region_definitions:
			child_regions[region] = Sequence._mate_bases(sequence1.region_definitions[region], sequence2.region_definitions[region])

		return Sequence(child_regions, sequence1.strand_structures)

	def build_strand(self, strand_structure):
		bases = ""
		for region in strand_structure:
			if region.name.islower():
				bases += self.region_definitions[region.name]
			else:
				bases += Strand.complement(self.region_definitions[region.name.lower()])

		return Strand(bases, strand_structure)

	def fitness(self, mfold):
		strands = [self.build_strand(strand_structure) for strand_structure in self.strand_structures]
		energy_matrix = EnergyMatrix(mfold, strands)
		energy_matrix.create()
		return exp(-np.linalg.norm(energy_matrix.matrix))


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
		self.mfold = Mfold(mfold_command='mfold')

	def iterate(self):
		# Find the fitness of each sequence in the population
		fitnesses = [sequence.fitness(self.mfold) for sequence in self.population]
		weighted_fitnesses = [fitness/sum(fitnesses) for fitness in fitnesses]

		# Mate strands at random, weighted by fitness level
		self.population = [self.generate_child(weighted_fitnesses) for i in range(self.population_size)]

		# Mutate the strands
		for sequence in self.population:
			sequence.mutate(self.mutation_rate)

	def _round_up(self, weights, number):
		curr = 0
		for i in range(len(weights)):
			curr += weights[i]
			if curr > number:
				return self.population[i]

		return self.population[-1]

	def generate_child(self, weighted_fitnesses):
		parent1 = self._round_up(weighted_fitnesses, random())
		parent2 = self._round_up(weighted_fitnesses, random())
		return Sequence.mate(parent1, parent2)

	def run(self):
		for i in range(self.iterations):
			self.iterate()



