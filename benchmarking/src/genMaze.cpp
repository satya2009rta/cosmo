/*
 * genMaze.cpp
 *
 *  A program to generate a maze
 * */

// #include "FileHandler.hpp"
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

int main(int argc, char* argv[]) {
    /* inputs */
    if (argc < 3){
        std::cerr << "Error: Not enough inputs!\n";
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

