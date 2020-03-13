# ENPH 479 Capstone Project 1961
## Genetic algorithm pipeline for optimizing DNA origami sequences.

This project is a pipeline for generating DNA origami sequences. Given an input origami structure, the pipeline uses a genetic algorithm to optimize a DNA sequence that will fold into the desired structure, and outputs the sequence for synthesis in a physical wet lab.

## Installation

First, make sure you have Python 3.6+, numpy, and matplotlib somehow installed on your system.

Then, you can simply run the convenient installation script `install.sh` if you use the `apt` package manager.

Otherwise, install Fortran and manually install mfold. Then replace the `mfold_quik` script in your `PATH` with the one in this repository.

## Usage

Run `python3 cli.py` and follow the instructions on the terminal. Your inputs are saved to `config.dat` by default, and you can run the algorithm again with the same parameters again by running `python3 cli.py config.dat`.

The input structure has the following format: `a25 B25, b25 C25, c25 D25, d25 A25`. This defines a sequence, which is composed of comma-separated DNA strands. Each DNA strand is partitioned into segments, and each segment is defined by a single alphabetic letter followed by a number specifying its length. An uppercase segment such as `A25` is complementary to its corresponding lowercase segment `a25`, and vice versa; for example, if `b4` is `AGGT`, then `B4` must be `TCCA`. 

The given example represents a 4-arm star. The sequence is composed of four separate strands, with the latter half of each strand binding to the former half of the next. If you take four strips of paper and glue them together like so, you will have created a star with four arms! In the physical wet lab, the complementary segments among the strands should bind to each other, so that the entire sequence folds into the desired structure. The goal of the genetic algorithm is to minimize the probability of pairs of strands folding into undesirable secondary structures.
