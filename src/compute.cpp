/*
 * compute.cpp
 *
 *  A program to find assump and strat template for 2-player games 
 * */

#include "FileHandler.hpp"

/* command line inputs:
 *
 *  [1]: input file name (if not given, by default, reads input as a hoa formatted game)
 */

int main(int argc, char* argv[]) {
    try {
        mpa::Game G;
        if (argc <= 1) { /* no filename, then read hoa formatted game from input */
            G = std2game();
        }
        else{
            /* read filename */
            std::string filename(argv[1]);
            G = file2game(filename);
        }

        /* print the game */
        distgame2std(G);
        std::cout << "======================================================\n";

        /* compute the strat template */
        mpa::Template assump;
        mpa::Template strat;
        auto winning_region = G.find_assump_parity(assump, strat);
        G.filter_out_edge_states(winning_region,assump,strat); /* remove edge-states from result (needned for HOA formatted games) */
        
        G.print_set(winning_region.first, "Winning Region");

        std::cout << "\n-- Assumptions --";
        assump.print_template();

        std::cout << "\n-- Strategy Template --";
        strat.print_template();
        std::cout << "*=====================================================\n";
        
        std::cout << "\n#winning_vertices:"<< winning_region.first.size()<< "/"<<winning_region.first.size()+winning_region.second.size()<<"\n\n";
        std::cout << "\n-- Assumptions --";
        assump.print_size();

        std::cout << "\n-- Strategy Template --";
        strat.print_size();
        
        if (winning_region.first.find(G.init_vert_) != winning_region.first.end()){
            std::cout << "REALIZABLE!\n";
            return 1;
        }
        std::cout << "UNREALIZABLE!\n";
        return 0;
    }
    /* TODO: print the same output in the output file */
    catch (const std::exception &ex) {
        std::cout << ex.what() << "\n";
        return 0;
    }
}

