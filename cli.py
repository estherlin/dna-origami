from genetic import GeneticAlgorithm
from mfold_library import Region
import matplotlib.pyplot as plt
import statistics
import sys

def parse_raw_structure(raw_structure):
	return [[Region(region[0], int(region[1:])) for region in strand] for strand in [strand.strip().split() for strand in raw_structure.split(',')]]


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

		print("Enter the file name for the output plot of fitness and diversity history: (default: history.png)")
		params["outfile"] = input()
		if not params["outfile"]:
			params["outfile"] = "history.png"
		print(f"Output plot will be saved to: {params['outfile']}\n")

		save_configuration(params)


	gen_alg = GeneticAlgorithm(
		parse_raw_structure(params["raw_structure"]),
		mfold_command=params["mfold_command"],
		population_size=int(params["population_size"]),
		iterations=int(params["iterations"]),
		mutation_rate=int(params["mutation_rate"])
	)
	gen_alg.run()
	print("Diversity history: ", gen_alg.diversity_history)
	print("Fitness history: ", gen_alg.fitness_history)
	with open("diversity.dat", "w") as outfile:
		outfile.write(gen_alg.diversity_history)
	with open("fitness.dat", "w") as outfile:
		outfile.write(gen_alg.fitness_history)

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
	axs[0].set_ylim([0,25])
	axs[0].set_xlim([0,100])

	axs[1].plot(iterations, std)
	axs[1].set_xlabel('Iteration')
	axs[1].set_ylabel('Standard deviation')
	axs[1].grid(True)
	axs[1].set_title('Standard deviation of norms in population per iteration')
	axs[1].set_ylim([0,5])
	axs[1].set_xlim([0,100])

	axs[2].plot(iterations, gen_alg.diversity_history)
	axs[2].axhline(y=12.5, color='r', linestyle='-')
	axs[2].set_xlabel('Iteration')
	axs[2].set_ylabel('Diversity')
	axs[2].grid(True)
	axs[2].set_title('Diversity of population per iteration')
	axs[2].set_ylim([0,15])
	axs[2].set_xlim([0,100])

	fig.tight_layout()
	plt.savefig(params["outfile"])
