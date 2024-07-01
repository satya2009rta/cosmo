/*
 * cosmo.cpp
 *
 *  A program to compute and compose assump and strat templates for 2-player distributed parity games */

#include "FileHandler.hpp"

void printHelpN(const std::string str) {
        std::cout << "cosmo " << str << " is not available.\n";
        std::cout << "Try cosmo --help for more information.\n";
    }
void printHelp() {
    std::cout << "Usage: cosmo [OPTION...]\n";
    std::cout << "Compute assumption and strategy templates as contracted strategy masks for both players in a two-objective parity game.\n";
    std::cout << "\nInputs/Outputs:\n";
    std::cout << "- STDIN: description of a (single-objective) parity game in extended-HOA/pgsolver format\n    or description of a two-objective parity game in extended-HOA/pgsolver format\n    or description of two parity games in extended-HOA format concated in one file\n";
    std::cout << "- STDOUT: assumption and strategy templates for each player (assumption on environment and strategy template for system in case of a single-objective parity game)\n"; 
    std::cout << "\nThe possible OPTIONs are as follows:\n";
    std::cout << "- --help                  Print this help message\n";
    std::cout << "- --print-game            Print the parity game (same format as input)\n";
    std::cout << "- --print-template-size   Print size of the templates\n";
    std::cout << "\nExample usage:\n";
    std::cout << "cosmo --print-template-size < example.gpg\n";
}

int main(int argc, char* argv[]) {
    try {
        bool print_game = false; // Flag to determine if game should be printed (same format as input)
        bool print_template_size = false; // Flag to determine if template size should be printed

        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--print-game") {
                print_game = true;
            } else if (std::string(argv[i]) == "--print-template-size") {
                print_template_size = true;
            } else if (std::string(argv[i]) == "--help") {
                printHelp();
                return 0;
            } else {
                printHelpN(std::string(argv[i]));
                return 0;
            }
        }



        /* construct games from stdin */
        mpa::DistGame G = std2distgame();

        /* print the game if print_game is true */
        if (print_game){
            distgame2std(G);
            std::cout << "===================================================\n";
        }

        /* if it is a single game, just compute assump for env and strat for system */
        if (G.n_games_ == 1){
            mpa::Game G2 = G;
            G2.vert_id_ = G.all_vert_id_[0];
            G2.colors_ = G.all_colors_[0];
            G2.max_color_ = G.all_max_color_[0];
            /* compute the strat template */
            mpa::Template assump;
            mpa::Template strat;
            auto winning_region = G2.find_assump_parity(assump, strat);
            G2.filter_out_edge_states(winning_region,assump,strat); /* remove edge-states from result (needned for HOA formatted games) */
            
            G2.print_set(winning_region.first, "Winning Region");

            std::cout << "\n-- Assumptions (on Envrionment)--";
            assump.print_template();

            std::cout << "\n-- Strategy Template (for System) --";
            strat.print_template();
            std::cout << "*=====================================================\n";
            
            /* print the size of the templates if print_template_size is true */
            if (print_template_size){
                std::cout << "\n#winning_vertices:"<< winning_region.first.size()<< "/"<<winning_region.first.size()+winning_region.second.size()<<"\n\n";
                std::cout << "\n-- Assumptions --";
                assump.print_size();

                std::cout << "\n-- Strategy Template --";
                strat.print_size();
                std::cout << "**===================================================\n";
            }
            
            if (winning_region.first.find(G2.init_vert_) != winning_region.first.end()){
                std::cout << "REALIZABLE!\n";
                return 1;
            }
            std::cout << "UNREALIZABLE!\n";
            return 0;
        }
        
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
        std::cout << "*===================================================\n";
        
        /* print the size of the templates if print_template_size is true */
        if (print_template_size){
            std::cout << "\n#winning_vertices:"<< winning_region.first.size()<< "/"<<winning_region.first.size()+winning_region.second.size()<<"\n";
            for (size_t i = 0; i < 2; i++){
                std::cout<<"\nPlayer "<<i;
                std::cout << "\n-- Assumptions --";
                assumps[i].print_size();

                std::cout << "\n-- Strategy Template --";
                strats[i].print_size();
            }
            std::cout <<"**==================================================\n";
        }

        if (winning_region.first.find(G.init_vert_) != winning_region.first.end()){
            std::cout << "REALIZABLE!\n";
            return 0;
        }
        else{
            std::cout << "UNREALIZABLE!\n";
            return 1;
        }
    }
    /* TODO: print the same output in the output file */
    catch (const std::exception &ex) {
        std::cout << ex.what() << "\n";
        printHelp();
        return -1;
    }
}

