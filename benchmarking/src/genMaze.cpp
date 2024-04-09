/*
 * genMaze.cpp
 *
 *  A program to generate a maze
 * */

#include "Maze.hpp"


/* command line inputs:
 *
 *  [1]: #coloumns of maze
 *  [2]: #rows of maze
 *  [3]: #walls in maze (optional, default is max possible)
 *  [4]: max #corridors in maze (optional, default is 0)
 *  [5]: parity game (HOA) for objective 0 (optional, default is GF ur0)
 *  [6]: parity game (HOA) for objective 1 (optional, default is GF ul1)
 */

void printHelp() {
    std::cout << "Usage: genMaze [num_cols] [num_rows] [num_walls] [num_corr] [obj1] [obj2]\n";
    std::cout << "Generate a maze benchmark as a two-objective parity game.\n";
    std::cout << "\nInputs/Outputs:\n";
    std::cout << "- num_cols: number of columns\n";
    std::cout << "- num_rows: number of rows\n";
    std::cout << "- num_walls: number of walls\n";
    std::cout << "- num_corr: number of corridors (optional, default is 0)\n";
    std::cout << "- obj1: objective of Player 1 as a parity game in extended-HOA format (optional, default is GF ur0)\n";
    std::cout << "- obj2: objective of Player 2 as a parity game in extended-HOA format (optional, default is GF ul1)\n";
    std::cout << "- STDOUT: a maze benchmark as a two-objective parity game\n"; 
    std::cout << "\nExample usage:\n";
    std::cout << "genMaze 3 3 2 1\n";
}

int main(int argc, char* argv[]) {
    try {
        /* inputs */
        if (argc < 3){
            printHelp();
            return 0;
        }
        if (std::string(argv[0]) == "--help") {
                printHelp();
                return 0;
        }
        size_t max_x = std::stoi(argv[1]);
        size_t max_y = std::stoi(argv[2]);
        size_t n_walls = max_x*max_y - max_x - max_y + 1;
        if (argc > 3){
            n_walls = std::stoi(argv[3]);
        }
        size_t max_corriders = 0;
        if (argc > 4){
            max_corriders = std::stoi(argv[4]);
        }
        mpa::Maze G;

        /* construct the maze */
        G.dummyMaze(max_x,max_y,n_walls,max_corriders);
        
        if (argc > 5){
            std::string filename1(argv[5]);
            std::string filename2(argv[6]);

            mpa::DistGame G1;
            /* construct the dist game for both LTL formulas */
            G1.product_games(hoa2game(filename2),hoa2game(filename1));
            /* take product of the maze and distgame */
            G.product_maze_game(G1);
        }
    
        distgame2hoa(G);
        return 0;
    }
    /* TODO: print the same output in the output file */
    catch (const std::exception &ex) {
        std::cout << ex.what() << "\n";
        printHelp();
        return -1;
    }
}

