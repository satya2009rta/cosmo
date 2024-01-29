/*
 * pg2dpg.cpp
 *
 *  A program to add two set of random colors to a (pgsolver) game to 
 *  get a dpg (gpg) format distributed parity game */

#include <functional>
#include "FileHandler.hpp"

/* command line inputs:
 *
 *  [1]: input file name
 *  [2]: maximum color
 */

int main(int argc, char* argv[]) {
    try {
        if (argc <= 2){
            std::cerr << "Not enough inputs!\n";
        }
        mpa::Game G1;
        /* read filename */
        std::string filename(argv[1]);
        G1 = file2game(filename);
        
        
        /* construct the dist-game */
        mpa::DistGame G(G1);
        
        int no_games = 2; /* number of new color sets to be added */
        int rand_max_col = std::stoi(argv[2]); /* maximum color that can be generated */
        
        G.all_colors_.clear();
        G.n_games_ = 0;
        G.randDistGame(no_games,rand_max_col,false);
        distgame2gpg(G);
        return 1;
    }
    /* TODO: print the same output in the output file */
    catch (const std::exception &ex) {
        std::cout << ex.what() << "\n";
        return -1;
    }
}

