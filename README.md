# CoSMo-**Co**ntracted-**S**trategy-**M**ask-Neg**o**tiation

CoSMo is a tool for computing assumption and strategy templates as contracted strategy masks for both players in a two-objective parity game, presented in the paper [Contract-Based Distributed Logical Controller Synthesis](https://arxiv.org/abs/2307.06212) (HSCC'24).

## Requirements

- <a href='https://gcc.gnu.org/install/'>GCC</a>

## Installation

Run the following commands to build all executable files

```
make
```
## Usage 

The executable files are generated and stored in the folder `./build/`. Usage of each executable is outlined below.

### - cosmo
```
Usage: cosmo [OPTION...]
```

Inputs/Outputs:
- STDIN: 
    - description of a (single-objective) parity game in extended-HOA/pgsolver format OR 
    - description of a two-objective parity game in extended-HOA/pgsolver format OR 
    - concatenated description of two parity games in extended-HOA format
    - description of a multi-objective parity game in pgsolver format (except for the 1st objective, all other objectives are considered for player 1)
- STDOUT: 
    - assumption and strategy templates for each player (in case of a two-objective parity game)
    - assumption on player 1 and strategy template for player 0 in case of a single-objective parity game

The possible OPTIONs are as follows:
- --print-game: print the parity game (same format as input)
- --print-template-size: print size of the templates

Example usage:
```
./build/cosmo --print-game < ./examples/test_dhoa_02.hoa
```

### - genMaze
```
Usage: genMaze [num_cols] [num_rows] [num_walls] [num_corr] [obj1] [obj2]
```

Inputs/Outputs:
- num_cols: number of columns
- num_rows: number of rows
- num_walls: number of walls
- num_corr: number of corridors (optional, default is 0)
- obj1: objective of Player 1 as a parity game in extended-HOA format (optional, default is GF ur0)
- obj2: objective of Player 2 as a parity game in extended-HOA format (optional, default is GF ul1)
- STDOUT: a maze benchmark as a two-objective parity game

Example usage:
```
./build/genMaze 3 3 2 1
```


