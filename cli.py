from genetic import GeneticAlgorithm, Sequence
from mfold_library import Region
import matplotlib.pyplot as plt
import statistics
import sys
import ast
import re

def parse_raw_structure(raw_structure):
	return [[Region(re.findall('\D+', region)[0], int(re.findall('\d+', region)[0])) for region in strand] for strand in [strand.strip().split() for strand in raw_structure.split(',')]]

def parse_raw_sequences(raw_sequences, structure):
	if isinstance(raw_sequences, str) and len(raw_sequences) <= 2:
		return []
	sequences = []

	r_s = raw_sequences
	if isinstance(raw_sequences, str):
		r_s = ast.literal_eval(raw_sequences)

	for raw_sequence in r_s:
		defs = {}
		raw_strands = raw_sequence.split(',')
		for i in range(len(raw_strands)):
			curr_ind = 0
			for region in structure[i]:
				if region.name.islower():
					defs[region.name] = raw_strands[i][curr_ind:(curr_ind + region.length)]
				curr_ind += region.length
		sequences += [Sequence(defs, structure)]
	return sequences


def consume_input(key, default):
	print(f"Enter the {key}: (default: {default})")
	value = input().strip()
	if not value:
		value = default
	print(f"Given {key}: {value}")
	return value


def save_configuration(params):
	print("Enter the file name to save your input configurations: (default: config.dat)")
	configpath = input()
	if not configpath:
		configpath = "config.dat"
	with open(configpath, "w") as configfile:
		for key, value in params.items():
			configfile.write(f"{key}: {value}\n")
	print(f"Configuration file saved to {configpath}.")
	print(f"You can edit the configuration file directly and run `python3 cli.py {configpath}` next time to skip the manual setup steps.")


def load_configuration(configpath):
	print(f"Automatically using inputs from configuration file {configpath}")
	params = {}
	with open(configpath, "r") as configfile:
		for line in configfile:
			split_line = line.split(':', 1)
			params[split_line[0].strip()] = split_line[1].strip()
	return params



if __name__ == '__main__':
	if len(sys.argv) > 1:
		params = load_configuration(sys.argv[1])

	else:
		params = {}

		print("Enter your desired shape (for example: a25 B25, b25 C25, c25 D25, d25 A25)")
		params["raw_structure"] = input().strip()
		print(f"Given desired shape: {parse_raw_structure(params['raw_structure'])}\n")

		params["mfold_command"] = consume_input('the path to Mfold executable', '~/.local/bin/mfold_quik')
		params["population_size"] = consume_input('population size', '25')
		params["mutation_rate"] = consume_input('mutation rate', '100')
		params["iterations"] = consume_input('number of iterations', '100')
		params["boltzmann_factor"] = consume_input('Boltzmann scaling factor', '1')
		num_init_seq = int(consume_input('number of initial sequences', '0'))
		params["raw_input_sequences"] = []
		if num_init_seq > 0:
			print("Enter one sequence per line (for example: AACG...,CCTG...,GGTA...)")
			for i in range(num_init_seq):
				params["raw_input_sequences"] += [input().strip()]

		print("Enter the file name for the output plot of fitness and diversity history: (default: history.png)")
		params["outfile"] = input()
		if not params["outfile"]:
			params["outfile"] = "history.png"
		print(f"Output plot will be saved to: {params['outfile']}\n")

		save_configuration(params)

	structure = parse_raw_structure(params["raw_structure"])
	gen_alg = GeneticAlgorithm(
		structure,
		mfold_command=params["mfold_command"],
		population_size=int(params["population_size"]),
		iterations=int(params["iterations"]),
		mutation_rate=int(params["mutation_rate"]),
		boltzmann_factor=float(params["boltzmann_factor"]),
		initial_sequences=parse_raw_sequences(params["raw_input_sequences"], structure)
	)
	try:
		gen_alg.run()
	finally:
		print(len(gen_alg.diversity_history) - 1, " iterations completed")
		print("Diversity history: ", gen_alg.diversity_history)
		print("Fitness history: ", gen_alg.fitness_history)
		with open("diversity.dat", "w") as outfile:
			for diversity in gen_alg.diversity_history:
				outfile.write(str(diversity) + '\n')
		with open("fitness.dat", "w") as outfile:
			for fitness in gen_alg.fitness_history:
				outfile.write(str(fitness) + '\n')
		gen_alg.print_population(final=True)

	iterations = range(int(params["iterations"]))
	best = [min(iteration) for iteration in gen_alg.fitness_history]
	worst = [max(iteration) for iteration in gen_alg.fitness_history]
	std = [statistics.stdev(iteration) for iteration in gen_alg.fitness_history]

	plt.rcParams["figure.figsize"] = [10, 15]
	fig, axs = plt.subplots(3, 1)

	axs[0].plot(iterations, best, 'r', label='Best solution')
	axs[0].plot(iterations, worst, 'b', label='Worst solution')
	axs[0].set_xlabel('Iteration')
	axs[0].set_ylabel('Norm')
	axs[0].grid(True)
	axs[0].legend()
	axs[0].set_title('Norms of best and worst solutions per iteration')
	axs[0].set_ylim([0.8*min(best),1.2*max(worst)])
	axs[0].set_xlim([min(iterations),max(iterations)])

	axs[1].plot(iterations, std)
	axs[1].set_xlabel('Iteration')
	axs[1].set_ylabel('Standard deviation')
	axs[1].grid(True)
	axs[1].set_title('Standard deviation of norms in population per iteration')
	axs[1].set_ylim([0, 1.2*max(std)])
	axs[1].set_xlim([min(iterations),max(iterations)])

	axs[2].plot(iterations, gen_alg.diversity_history)
	axs[2].axhline(y=12.5, color='r', linestyle='-')
	axs[2].set_xlabel('Iteration')
	axs[2].set_ylabel('Diversity')
	axs[2].grid(True)
	axs[2].set_title('Diversity of population per iteration')
	axs[2].set_ylim([0,1.2*max(gen_alg.diversity_history)])
	axs[2].set_xlim([min(iterations),max(iterations)])

	fig.tight_layout()
	plt.savefig(params["outfile"])
