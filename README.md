# CoSMo-**Co**ntracted-**S**trategy-**M**ask-Neg**o**tiation

CoSMo is a tool for computing assumption and strategy templates as contracted strategy masks for both players in a two-objective parity game, presented in the paper [Contract-Based Distributed Logical Controller Synthesis](https://arxiv.org/abs/2307.06212) (HSCC'24).

## Requirements

- <a href='https://gcc.gnu.org/install/'>GCC</a>

## Installation

Run the following commands to build all executable files

```
make build
```
## Usage 

The executable files are generated and stored in the folder `./build/`. Usage of each executable is outlined below.

- `compute` computes an assumption template for the environment and a strategy template for the system in a parity game. It requires stdin input which is the description of a parity game in extended HOA or pgsolver format; and outputs the result to stdout.

- `negotiate` computes an assumption template and a strategy template for each player in a two-objective parity game. It requires stdin input which is the description of a generalized parity game in extended HOA or pgsolver format; and outputs the result to stdout.

- `pg2dpg` converts a parity game to a two-objective parity game by adding two random parity objectives to a game graph. It requires two command-line arguments: (1) filename that contains the description of a parity game in extended HOA or pgsolver format, (2) maximum priority of the parity objectives; and outputs the result game in pgsolver format to stdout.

- `genMaze` generates a maze benchmark as a two-objective parity game. It requires six command-line arguments: (1) number of columns, (2) number of rows, (3) number of walls, (4) number of corridors (optional, default is 0), (5) objective of Player 1 as a parity game in extended HOA format (optional, default is GF ur0), (6) objective of Player 2 as a parity game in extended HOA format (optional, default is GF ul1); and outputs the result game in extended HOA format to stdout. 



