/*
 * negotiate.cpp
 *
 *  A program to compute and compose assump and strat templates for 2-player distributed games */

#include "FileHandler.hpp"

/* command line inputs:
 *
 *  [1]: input file name (if not given, by default, reads input as a hoa formatted games)
 */

int main(int argc, char* argv[]) {
    try {
        mpa::DistGame G;
        if (argc == 2){ /* dist-game given in a single file using gpg format */
            /* read filename */
            std::string filename(argv[1]);
            /* convert the distributed game in gpgsolver format to normal dist-game */
            G = gpg2distgame(filename);
        }
        else if (argc > 2){/* dist-game is product of two games given in two files in HOA format */
                std::string filename(argv[1]);
                std::string filename2(argv[2]);
                
                /* construct the dist game */
                G.product_games(hoa2game(filename),hoa2game(filename2));
                // G = hoa2distgame(filename,filename2);
            }
        else{/* construct games from stdin */
            G = std2distgame();
        }

        /* print the game */
        distgame2std(G);
        
        std::cout << "======================================================\n";
        
        /* compute the strat template */
        std::pair<std::set<size_t>, std::set<size_t>> winning_region;

        std::vector<mpa::Template> assumps;
        std::vector<mpa::Template> strats;
        winning_region = G.find_AG_contract(assumps, strats);
        /* remove edge-states from result (needned for HOA formatted games) */
        G.filter_out_edge_states(winning_region,assumps,strats);

        G.print_set(winning_region.first, "Winning Region");
        for (size_t i = 0; i < 2; i++){
            std::cout<<"\nPlayer "<<i;
            std::cout << "\n-- Assumptions --";
            assumps[i].print_template();
            
            std::cout << "\n-- Strategy Template --";
            strats[i].print_template();
        }
        std::cout << "***===================================================\n";
        
        
        std::cout << "\n#winning_vertices:"<< winning_region.first.size()<< "/"<<winning_region.first.size()+winning_region.second.size()<<"\n";
        for (size_t i = 0; i < 2; i++){
            std::cout<<"\nPlayer "<<i;
            std::cout << "\n-- Assumptions --";
            assumps[i].print_size();

            std::cout << "\n-- Strategy Template --";
            strats[i].print_size();
        }
        std::cout <<"****==================================================\n";
        if (winning_region.first.find(G.init_vert_) != winning_region.first.end()){
            std::cout << "REALIZABLE!\n";
            return 1;
        }
        else{
            std::cout << "UNREALIZABLE!\n";
            return 0;
        }
    }
    /* TODO: print the same output in the output file */
    catch (const std::exception &ex) {
        std::cout << ex.what() << "\n";
        return -1;
    }
}

