from genetic import GeneticAlgorithm
from mfold_library import Region
import matplotlib.pyplot as plt
import statistics

ITERATIONS=100
POPULATION_SIZE=25

gen_alg = GeneticAlgorithm(
	[[Region('a', 25), Region('B', 25)], [Region('b', 25), Region('C', 25)], [Region('c', 25), Region('D', 25)], [Region('d', 25), Region('A', 25)]],
	mfold_command='/home/ubuntu/dna-origami/mfold_quik',
	population_size=POPULATION_SIZE,
	iterations=ITERATIONS,
	mutation_rate=100
)

gen_alg.run()
print(gen_alg.diversity_history)
print(gen_alg.fitness_history)

iterations = range(ITERATIONS)
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

plt.savefig('history.png')
