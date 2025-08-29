/*
 * Class: Game
 *
 *  Class formalizing a basic two player (parity) game
 */

#ifndef GAME_HPP_
#define GAME_HPP_

#include "Template.hpp"

#define V0 0 /* vertices belonging to player 0 */
#define V1 1 /* vertices belonging to player 1 */

namespace mpa {
class Game {
public:
    /* number of vertices */
    size_t n_vert_;
    /* number of edges */
    size_t n_edge_;
    /* vertices */
    std::set<size_t> vertices_;
    /* initial vertex */
    size_t init_vert_;
    /* vertex id: V0, when the vertex belongs to player 0 and V1 when it belongs to player 1 */
    std::map<size_t, size_t> vert_id_;
    /* edges as a map from vertices to set of its neighbours */
    std::map<size_t, std::set<size_t>> edges_;
    /* maximum of colors */
    size_t max_color_;
    /* colors of vertices */
    std::map<size_t, size_t> colors_;

    /* variables needed for ehoa formatted games */
    /* pre of edges for original vertices (incoming edges) */
    std::map<size_t, std::set<size_t>> pre_edges_;
    /* ids of atomic proposition */
    std::map<size_t, std::string> ap_id_;
    /* labels of mid-state */
    std::map<size_t, std::vector<size_t>> labels_;
    /* controllable APs */
    std::set<size_t> controllable_ap_;
    /* name of states */
    std::map<size_t, std::string> state_names_;
public:
    /* default constructor */
    Game() {
        n_vert_ = 0;
        n_edge_ = 0;
    }


    ///////////////////////////////////////////////////////////////
    /// Basic game change operations
    ///////////////////////////////////////////////////////////////

    /* convert a parity odd game parity even game by decreasing the colors by 1 */
    int oddToEvenGame(){
        for (auto& pair : colors_){
            colors_[pair.first] = pair.second - 1;
        }
        return 1;
    }
    /* complement vertex ids of vertices */
    int complementPlayers(){
        for (auto& pair : vert_id_){
            if (vert_id_[pair.first] == 1){
                vert_id_[pair.first] = 0;
            }
            else if (vert_id_[pair.first] == 0){
                vert_id_[pair.first] = 1;
            }
        }
        return 1;
    }

    ///////////////////////////////////////////////////////////////
    /// Solving all types of games
    ///////////////////////////////////////////////////////////////

    /* solve reachability game for player i (default is player 0)
     * input: target
     * output: Reach_i(target) */
    std::pair<std::set<size_t>, std::set<size_t>> solve_reachability_game(const std::set<size_t>& target,
                                                    const std::set<size_t>& players = {V0}) const {
        std::set<size_t> losing; /* vertices from which targets might not be reachable */
        std::set<size_t> winning = target; /* vertices from which targets are currently reachable */
        while (true)
        {
            losing = set_complement(winning); /* complement of vertices from which targets are currently reachable */
            std::set<size_t> new_winning; /* new vertices from which targets are currently reachable */
            for (const size_t& v: losing){ /* new_winning =  cpre_i(winning) */
                if(players.find(vert_id_.at(v)) != players.end() && check_set_intersection(edges_.at(v),winning))
                    new_winning.insert(v);
                else if(check_set_inclusion(edges_.at(v),winning) && edges_.find(v) != edges_.end())
                    new_winning.insert(v);
            }
            if (new_winning.empty()){ /* if there is no new winning vertex then return winning */
                losing = set_complement(winning);
                return std::make_pair(winning, losing);
            }
            set_merge(winning, new_winning); /* include the new targets in winning */
        }
    }

    /* solve Buechi game 
     * input: target
     * output: winning region for player 0 */
    std::pair<std::set<size_t>, std::set<size_t>> solve_buchi_game(const std::set<size_t>& target) const {
        Game subgame = *this; /* copy the game */
        std::set<size_t> non_target = set_complement(target); /* vertices which are not target */
        /* in the subgame (copy of the game), target vertices has color 2 
        and non-target vertices have color 1 */
        for (const size_t& v: target)
            subgame.colors_.at(v) = 2;
        for (const size_t& v: non_target)
            subgame.colors_.at(v) = 1;
        return subgame.solve_parity_game(); /* solve the 2-color parity game using parity algo */
    }

    /* solve Co-Buechi game 
     * input: target (Eventually Always [target])
     * output: winning region for player 0 */
    std::pair<std::set<size_t>, std::set<size_t>> solve_cobuchi_game(const std::set<size_t>& target) const {
        Game subgame = *this; /* copy the game */
        std::set<size_t> non_target = set_complement(target); /* vertices which are not target */
        /* in the subgame (copy of the game), target vertices has color 2 
        and non-target vertices have color 1 */
        for (const size_t& v: target)
            subgame.colors_.at(v) = 0;
        for (const size_t& v: non_target)
            subgame.colors_.at(v) = 1;
        return subgame.solve_parity_game(); /* solve the 2-color parity game using parity algo */
    }


    /* solve parity game using Zielonka's algorithm
     * input: game with colors
     * output: winning region for player 0 */
    std::pair<std::set<size_t>, std::set<size_t>> solve_parity_game() const {
        Game game_copy(*this); /* copy of the game as we will change it */
        return game_copy.recursive_zielonka_parity();
    }
    /* recursive zielonka's algorithm */
    std::pair<std::set<size_t>, std::set<size_t>> recursive_zielonka_parity() const {
        if (n_vert_ == 0) /* if current region is empty, nothing to do, return empty */
            return std::make_pair(vertices_, vertices_);
        /* vertices with maximum color */
        std::set<size_t> max_col_vertices = vertex_with_color(max_color_);
        if (max_color_ % 2 == 1){ /* when max_color_ is odd */
            /* vertices from which player 1 can force to visit max_col_vertices (odd) */
            auto regionA = solve_reachability_game(max_col_vertices, {V1});
            /* complement of regionA */
            Game gameCA(subgame(regionA.second)); /* game with regionA removed */
            /* solve the gameCA */
            auto winCA = gameCA.recursive_zielonka_parity();
            if (winCA.first.empty()) /* if winning region is empty in gameCA then return everything empty */
                return std::make_pair(winCA.first, vertices_);
            else {
                /* vertices from which player 0 can force to reach winning region of gameCA */
                auto regionB = solve_reachability_game(winCA.first, {V0});
                Game gameCB(subgame(regionB.second)); /* game with regionB removed */
                /* solve gameCB */
                auto winCB = gameCB.recursive_zielonka_parity();
                return std::make_pair(set_union(winCB.first, regionB.first), winCB.second); /* winning region is regionB \cup winning region of gameCB */
            }
        }
        else{ /* max_color_ is even */
            /* vertices from which player 0 can force to visit max_col_vertices (even) */
            auto regionA = solve_reachability_game(max_col_vertices, {V0});
            /* complement of regionA */
            Game gameCA(subgame(regionA.second)); /* game with regionA removed */
            /* solve the gameCA */
            auto winCA = gameCA.recursive_zielonka_parity();
            if (winCA.second.empty()){ /* if losing region is empty */
                return std::make_pair(vertices_, winCA.second); /* winning region is whole set */
            }
            else {
                /* vertices from which player 1 can force to reach losing region of gameCA */
                auto regionB = solve_reachability_game(winCA.second, {V1});
                Game gameCB(subgame(regionB.second)); /* game with regionB removed */
                /* solve gameCB */
                auto winCB = gameCB.recursive_zielonka_parity();
                return std::make_pair(winCB.first, set_union(winCB.second, regionB.first)); /* winning region is winning region of gameCB */
            }
        }
    }



    ///////////////////////////////////////////////////////////////
    /// Solving all types of Cooperatvive games
    ///////////////////////////////////////////////////////////////

    /* solve cooperative safety game
     * input: target
     * output: winning region */
    std::pair<std::set<size_t>, std::set<size_t>> solve_coop_safety_game(const std::set<size_t>& target) const {
        /* non-targets is the complement of target */
        std::set<size_t> non_target = set_complement(target);
        /* losing region in coop_safety(target) is winning region in forced_reach(non_target) */
        auto losing_region = solve_reachability_game(non_target, {});
        /* return the complement of losing_region */
        return std::make_pair(losing_region.second, losing_region.first);
    }

    /* solve cooperative co-Buchi game
     * input: target (Eventually Always [target])
     * output: winning region */
    std::pair<std::set<size_t>, std::set<size_t>> solve_coop_cobuchi_game(const std::set<size_t>& target) const {
        auto safety_winning_region = solve_coop_safety_game(target); /* winning region of coop_safety(target) */
        /* winning region is the vertices from which players can cooperatively reach safety_winning(target) */
        return solve_reachability_game(safety_winning_region.first, {V0,V1});
    }

    /* solve cooperative Buchi game 
     * input: target
     * output: winning region */
    std::pair<std::set<size_t>, std::set<size_t>> solve_coop_buchi_game(const std::set<size_t>& target) const {
        /* make a copy of the game and make every vertex player 0 vertex */
        Game game_copy = *this;
        for (auto& pair : game_copy.vert_id_){
            pair.second = V0;
        }
        /* solve the normal buchi game on this copy */
        return game_copy.solve_buchi_game(target);
    }

    /* solve cooperative parity game 
     * output: winning region */
    std::pair<std::set<size_t>, std::set<size_t>> solve_coop_parity_game() const {
        /* make a copy of the game and make every vertex player 0 vertex */
        Game game_copy = *this;
        for (auto& pair : game_copy.vert_id_){
            pair.second = V0;
        }
        /* solve the normal buchi game on this copy */
        return game_copy.solve_parity_game();
    }




    ///////////////////////////////////////////////////////////////
    /// assump in Buechi games
    ///////////////////////////////////////////////////////////////
    
    /* compute the assump and strat template for Buechi game */
    std::pair<std::set<size_t>, std::set<size_t>> find_assump_buchi(const std::set<size_t>& target,
                                Template& assump,
                                Template& strat) const {
        /* initialize the template */
        assump.clear();
        strat.clear();

        /* first compute the cooperative winning region */
        auto coop_winning_region = find_live_groups_reach(target, assump, strat);
        /* assump_unsafe edges are the player 1's edges from cooperative winning region to losing region */
        assump.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V1});
        /* strat_unsafe edges are the player 0's edges from cooperative winning region to losing region */
        strat.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V0});
        return coop_winning_region;  /* return the (cooperative) winning region */
    }

    /* find live groups recursively for reaching target */
    std::pair<std::set<size_t>, std::set<size_t>>  find_live_groups_reach(const std::set<size_t>& target,
                                Template& assump,
                                Template& strat) const {
        std::set<size_t> curr_target = target;
        std::pair<std::set<size_t>, std::set<size_t>> curr_winning_region;
        /* keep finding live groups until convergence to the winning region */
        while (1){
            /* keep finding live groups until curr_winning_region = curr_target */
            while (1) {
                /* vertices from which no player can stop reaching cuurent winning region */
                curr_winning_region = solve_reachability_game(curr_target, {});
                curr_target = curr_winning_region.first; /* current target is the winning region of last reachability game */
                /* add strat_live group containing player0 edges from outiside to winning region; and add their sources to curr_target */
                std::map<size_t, std::set<size_t>> new_strat_live_edges = edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, {V0});
                if (curr_winning_region.first.size() == curr_target.size()) {
                    break; /* break the loop if curr_winning_region = curr_target */
                }
                strat.live_groups_.push_back(new_strat_live_edges);
            }
            /* add assump_live group containing player1 edges from outiside to winning region; and add their sources to curr_target */
            std::map<size_t, std::set<size_t>> new_assump_live_edges = edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, {V1});
            if (curr_winning_region.first.size() == curr_target.size()) {
                break; /* break the loop if curr_winning_region = curr_target */
            }
            assump.live_groups_.push_back(new_assump_live_edges);
        }
        return curr_winning_region; /* return current winning region */
    }


    ///////////////////////////////////////////////////////////////
    /// assump in Co-Buechi games
    ///////////////////////////////////////////////////////////////
    
    /* compute the assump and strat template for co-Buechi game */
    std::pair<std::set<size_t>, std::set<size_t>> find_assump_cobuchi(const std::set<size_t>& target,
                                Template& assump,
                                Template& strat) const {
        /* initialize the template */
        assump.clear();
        strat.clear();

        /* first compute the cooperative winning region */
        auto coop_winning_region = find_colive_edges_cobuchi(target, assump, strat);
        /* assump_unsafe edges are the player 1's edges from cooperative winning region to losing region */
        assump.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V1});
        /* strat_unsafe edges are the player 0's edges from cooperative winning region to losing region */
        strat.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V0});
        return coop_winning_region; /* return the (cooperative) winning region */
    }

    /* find colive edges recursively to eventually stay in target */
    std::pair<std::set<size_t>, std::set<size_t>>  find_colive_edges_cobuchi(const std::set<size_t>& target,
                                Template& assump,
                                Template& strat) const {
        auto safety_winning_region = solve_coop_safety_game(target); /* safety winning region for target */
        /* initial colive edges are the edges from safety winning region to losing region */
        edges_between(safety_winning_region.first, safety_winning_region.second, assump.colive_edges_, {V1});
        edges_between(safety_winning_region.first, safety_winning_region.second, strat.colive_edges_, {V0});

        std::set<size_t> curr_target = safety_winning_region.first;
        std::pair<std::set<size_t>, std::set<size_t>> curr_winning_region;
        /* keep finding colive edges until convergence to the winning region */
        while (1){
            /* vertices from which no player can stop reaching cuurent winning region */
            curr_winning_region = solve_reachability_game(curr_target, {});
            curr_target = curr_winning_region.first; /* current target is the winning region of last reachability game */
            /* mark player-i edges from outiside not going to winning region as colive (when there is a way to go to winning region) */
            co_edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, assump.colive_edges_, {V1});
            co_edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, strat.colive_edges_, {V0});
            if (curr_winning_region.first.size() == curr_target.size()) {
                break; /* break the loop if curr_winning_region = curr_target */
            }



            // while (1) {
            //     /* vertices from which no player can stop reaching cuurent winning region */
            //     curr_winning_region = solve_reachability_game(curr_target, {});
            //     curr_target = curr_winning_region.first; /* current target is the winning region of last reachability game */
            //     /* mark player0 edges from outiside to winning region as strat_colive; and add their sources to curr_target */
            //     co_edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, strat.colive_edges_, {V0});
            //     if (curr_winning_region.first.size() == curr_target.size()) {
            //         break; /* break the loop if curr_winning_region = curr_target */
            //     }
            // }
            // /* mark player1 edges from outiside to winning region as assump_colive; and add their sources to curr_target */
            // co_edges_between(curr_winning_region.second, curr_winning_region.first, curr_target, assump.colive_edges_, {V1});
            // if (curr_winning_region.first.size() == curr_target.size()) {
            //     break; /* break the loop if curr_winning_region = curr_target */
            // }
        }
        return curr_winning_region; /* return current winning region */
    }



    ///////////////////////////////////////////////////////////////
    /// assump in Parity games
    ///////////////////////////////////////////////////////////////

    /* compute the assump and strat template for parity game */
    std::pair<std::set<size_t>, std::set<size_t>> find_assump_parity(Template& assump, Template& strat) const {
        /* initialize the template */
        assump.clear();
        strat.clear();

        std::set<size_t> colive_vertices; /* initialize this set to use in recursive function */
        /* recursively solve and compute assumptions */
        auto coop_winning_region = find_colive_cond_live_parity(assump, strat , colive_vertices);
        /* assump_unsafe edges are the player 1's edges from cooperative winning region to losing region */
        assump.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V1});
        /* strat_unsafe edges are the player 0's edges from cooperative winning region to losing region */
        strat.unsafe_edges_ = edges_between(coop_winning_region.first, coop_winning_region.second, {V0});
        return coop_winning_region;
    }

    /* find colive edges and conditional live groups recursively for parity */
    std::pair<std::set<size_t>, std::set<size_t>> find_colive_cond_live_parity(Template& assump, Template& strat,               
                                std::set<size_t>& colive_vertices) const {
        if (n_vert_ == 0) /* if current region is empty, nothing to do, return empty */
            return std::make_pair(vertices_, std::set<size_t> {});
        /* vertices with maximum color */
        // std::cout << "n_vert:"<< n_vert_<<"\n";
        std::set<size_t> max_col_vertices = vertex_with_color(max_color_);
        // print_set(max_col_vertices, "max_col_vertices");
        if (max_color_ % 2 == 1){ /* when max_color_ is odd */
            // std::cout << "odd max_col:"<<max_color_ <<"\n";
            auto forced_reach_max_col = solve_reachability_game(max_col_vertices, {}); /* vertices from which max_col is visited anyway */
            Game game_new = *this; /* new game to solve parity without visiting max_col (odd) vertices */
            game_new.remove_vertices(forced_reach_max_col.first);
            /* vertices from which players can cooperatively win without visiting max_col_vertices */
            auto cobuchi_winning_region = game_new.find_colive_cond_live_parity(assump, strat, colive_vertices);
            /* a winning play eventually stays in the cobuchi_winning_region
            as it does not visit (odd) max_color vertices infinitely often*/
            // print_set(cobuchi_winning_region.first,"cobuchi_win");
            auto coop_winning_region = find_colive_edges_cobuchi(cobuchi_winning_region.first, assump, strat); /* cooperative winning region */
            // print_set(coop_winning_region.first);
            set_merge(colive_vertices, set_difference(coop_winning_region.first, cobuchi_winning_region.first)); /* the vertices outside of cobuchi_winning (but in coop_winning) should not be visited inf often */
            return coop_winning_region; /* return the winning region */
        }
        else{ /* max_color_ is even */
            // std::cout << "even max_col:"<<max_color_ <<"\n";
            /* vertices in maybe_winning_region from which one can cooperatively
            visit (even) max_color_vertices infinitely often */
            auto buchi_winning_region = solve_coop_buchi_game(max_col_vertices);
            if (!buchi_winning_region.first.empty()){
                // std::cout << "buchi winning region NOT empty!\n";
                // print_set(buchi_winning_region.first, "buchi winning region");
                Game game_new = *this;
                /* if buchi_winning_region is not empty, first we deal with that region */
                game_new.remove_vertices(buchi_winning_region.second);
                /* iterate over every odd color from 1 to max_col-1 */
                for (size_t col = 1; col < max_color_; ){
                    /* vertices with color = col */
                    std:: set<size_t> col_vertices = vertex_with_color(col, buchi_winning_region.first);
                    // std::cout << "vertices of color "<< col;
                    // print_set(col_vertices);
                    /* if col_vertices is non-empty, first compute vertices with even color
                    that is greater than col and solve buchi with them as target */
                    if (!col_vertices.empty()){
                        /* vertices with even color that is greater than col */
                        std::set<size_t> col_target;
                        /* iterate over every even color from col to max_col */
                        for (size_t win_col = col+1; win_col <= max_color_; ) {
                            std:: set<size_t> win_col_vertices = vertex_with_color(win_col, buchi_winning_region.first);
                            set_merge(col_target, win_col_vertices);
                            win_col += 2;
                        }
                        // print_set(col_vertices, "col_vertices");
                        // print_set(col_target,"col_target");
                        Template col_assump;
                        Template col_strat;
                        /* if a play visit col_vertices infinitely often,
                        it should visit col_target (vertices with even and > col color) infinitely often*/
                        auto col_winning_region = game_new.find_live_groups_reach(col_target, col_assump, col_strat);
                        if (!col_assump.live_groups_.empty()){
                            assump.cond_sets_.push_back(col_vertices);
                            assump.cond_live_groups_.push_back(col_assump.live_groups_);
                        }
                        if (!col_strat.live_groups_.empty()){
                            strat.cond_sets_.push_back(col_vertices);
                            strat.cond_live_groups_.push_back(col_strat.live_groups_);
                        }
                    }
                    col += 2; /* only consider odd colors col in the loop */
                }
            }
            if (!buchi_winning_region.second.empty()){
                // std::cout << "buchi losing region NOT empty!\n";
                // print_set(buchi_winning_region.second, "buchi losing region");
                /* if buchi_losing_region is non-empty,
                 finally we only need to deal with them */
                Game game_copy(subgame(buchi_winning_region.second)); /* the subgame restricted buchi_losing_region */
                max_col_vertices = vertex_with_color(max_color_, buchi_winning_region.second);
                /* a play can not visit (even) max_col infinitely often,
                so we change their color to 0 and solve odd degree color parity game */
                if (!max_col_vertices.empty()){
                    for (const auto& v : max_col_vertices){
                        game_copy.colors_.at(v) = 0;
                    }
                    game_copy.max_color_ = max_col(game_copy.colors_);
                }
                /* there is no unsafe edges,
                so we start with find_live_depressed_assumption directly */
                auto second_winning_region = game_copy.find_colive_cond_live_parity(assump, strat, colive_vertices);
                /* cooperative winning region is the union of buchi_winning and second winning */
                return (std::make_pair(set_union(buchi_winning_region.first, second_winning_region.first), second_winning_region.second));
            }
            /* if buchi_losing is empty then everything is cooperative winning */
            return buchi_winning_region;
        }
    }

    ///////////////////////////////////////////////////////////////
    /// Filter all edge-states out (needed for HOA formatted games)
    ///////////////////////////////////////////////////////////////

    /* replace edge-states by its (only) successor in a set of edges */
    int filter_edges(std::map<size_t, std::set<size_t>>& edges) const {
        for (auto& pair : edges){
            std::set<size_t> new_succs;
            for (const auto& succ : pair.second){
                new_succs.insert(*edges_.at(succ).begin());
            }
            pair.second = new_succs;
        }
        return 1;
    }

    /* filter out edge-states in a template */
    int filter_templates(Template& assump, const std::set<size_t>& org_vertices) const {
        filter_edges(assump.unsafe_edges_);
        filter_edges(assump.colive_edges_);
        for (auto& live_group : assump.live_groups_){
            filter_edges(live_group);
        }
        for (auto& live_groups : assump.cond_live_groups_){
            for (auto& live_group : live_groups){
                filter_edges(live_group);
            }
        }
        for (auto& cond_sets : assump.cond_sets_){
            cond_sets = set_intersection(cond_sets, org_vertices);
        }
        
        return 1;
    }

    /* filter out edge-states from a winning region, assumption and strategy template */
    int filter_out_edge_states(std::pair<std::set<size_t>, std::set<size_t>>& winning_region, Template& assump, Template& strat) const {
        if (labels_.empty()){
            return 0;
        }
        std::set<size_t> org_vertices;
        for (const auto& v : vertices_){
            if (vert_id_.at(v) != 2){
                org_vertices.insert(v);
            }
        }
        winning_region.first = set_intersection(winning_region.first, org_vertices);
        winning_region.second = set_intersection(winning_region.second, org_vertices);
        
        filter_templates(assump, org_vertices);
        filter_templates(strat, org_vertices);
        return 1;
    }

    ///////////////////////////////////////////////////////////////
    ///Basic functions
    ///////////////////////////////////////////////////////////////
    
    /* function: set_difference
     *
     * compute set difference of two sets*/

    std::set<size_t> set_difference(const std::set<size_t>& set2, const std::set<size_t>& set1) const {
        std::set<size_t> set3; /* set2 - set1 */
        for (const auto& u : set2){
            if (set1.find(u) == set1.end()){
                set3.insert(u);
            }
        }
        return set3;
    }

    /* function: set_complement
     *
     * compute complement of a set*/
    
    std::set<size_t> set_complement(const std::set<size_t>& set1) const {
        return set_difference(vertices_, set1);
    }

    /* function: check_set_inclusion
     *
     * check if set1 is included in set2 */
    template<class T>
    bool check_set_inclusion(const std::set<T>& set1, const std::set<T>& set2) const {
        if (set2.empty()){
            return false;
        }
        for (auto a = set1.begin(); a != set1.end(); ++a) {
            if (set2.find(*a) == set2.end()) {
                return false;
            }
        }
        return true;
    }

    /* function: check_set_intersection
     *
     * check if there is nonempty intersection between the two sets */
    template<class T>
    bool check_set_intersection(const std::set<T>& set1, const std::set<T>& set2) const {
        for (auto a = set1.begin(); a != set1.end(); ++a) {
            if (set2.find(*a) != set2.end()) {
                return true;
            }
        }
        return false;
    }


    /* function: set_merge
     *
     * merge one set into another */
    
    void set_merge(std::set<size_t>& set1, const std::set<size_t>& set2) const {
        set1.insert(set2.begin(), set2.end());
    }

    /* function: vertex_with_color
     *
     * returns the set of vertices (in a set) with color c */
    std::set<size_t> vertex_with_color(const size_t& c) const {
        std::set<size_t> result;
        for (const auto& v : vertices_){
            if (colors_.at(v) == c){
                result.insert(v);
            }
        }
        return result;
    }
    std::set<size_t> vertex_with_color(const size_t& c, const std::set<size_t>& set) const {
        std::set<size_t> result;
        for (const auto& v : set){
            if (colors_.at(v) == c){
                result.insert(v);
            }
        }
        return result;
    }

    /* function: map_remove_values
     *
     * remove set of values from one map (set of edges) */

    void map_remove_values(std::map<size_t, std::set<size_t>>& map, const std::set<size_t>& values) const {
        for (auto v = map.begin(); v != map.end(); v++){
            v->second = set_difference(v->second, values);
        }
    }
    void map_remove_values(std::map<size_t, size_t>& map, const std::set<size_t>& values) const {
        for (auto e = map.begin(); e != map.end();) {
            if (values.find(e->second) != values.end())
                map.erase(e++);
            else
                ++e;
        }
    }

    /* function: map_remove_keys
     *
     * remove set of keys from one map  */
    template<class T>
    void map_remove_keys(std::map<size_t, T>& map, const std::set<size_t>& keys) const {
        for (auto const & key : keys){
            map.erase(key);
        }
    }

    /* function: map_remove_vertices
     *
     * remove set of keys and values from one map */
    void map_remove_vertices(std::map<size_t, std::set<size_t>>& map, const std::set<size_t>& keys) const {
        map_remove_keys(map, keys);
        map_remove_values(map,keys);
    }

    /* function: max_col
     *
     * compute the max_color the game */
    size_t max_col(const std::map<size_t,size_t>& colors) const {
        size_t max_color = 0;
        for (const auto& pair : colors){
            if (pair.second > max_color)
                max_color = pair.second;
        }
        return max_color;
    }

    /* function: remove_vertices
     *
     * remove vertices from the game */
    
    void remove_vertices(const std::set<size_t>& set) {
        vertices_= set_complement(set);
        n_vert_ = vertices_.size();
        map_remove_keys(vert_id_, set);
        
        map_remove_keys(colors_, set);
        max_color_ = max_col(colors_);
        
        map_remove_keys(edges_, set);
        n_edge_ = 0;
        for (const auto& v : vertices_){
            edges_.at(v) = set_difference(edges_.at(v),set);
            n_edge_ += edges_.at(v).size();
        }
        /* sanity check */
        // valid_game();
    }

    /* function: subgame
     *
     * returns a subgame restricted to this set */
    
    Game subgame(const std::set<size_t>& set) const{
        Game game(*this);
        game.n_vert_ = set.size();
        game.vertices_= set;

        std::set<size_t> complement = set_complement(set);
        map_remove_keys(game.vert_id_, complement);
        map_remove_keys(game.colors_, complement);
        game.max_color_ = max_col(game.colors_);
        

        
        map_remove_keys(game.edges_, complement);
        game.n_edge_ = 0;
        for (const auto& v : game.vertices_){
            game.edges_.at(v) = set_difference(game.edges_.at(v),complement);
            game.n_edge_ += game.edges_.at(v).size();
        }
        /* sanity check */
        // game.valid_game();
        
        return game;
    }

    /* check valid_game */
    void valid_game() const {
        if (n_vert_ != vertices_.size() || n_vert_ != edges_.size() || n_vert_ != vert_id_.size() || n_vert_ != colors_.size()){
            std::cerr << "ERROR: Invalid Game ###############.\n";
            std::cout << "vertex:"<<n_vert_ <<"," << vertices_.size()<<"   edges"<< n_edge_ <<","<<edges_.size()<< "  vert_id:" << vert_id_.size() << " color:"<<colors_.size()<<"\n";
        }
    }
    
    /* function: set_intersetion
     *
     * compute intersection of sets */
    std::set<size_t> set_intersection(const std::set<size_t>& set1, const std::set<size_t>& set2, const std::set<size_t>& set3 = std::set<size_t>{}) const {
        std::set<size_t> set4;
        if (set3.empty()){
            for (const auto& a : set1) {
                if (set2.find(a) != set2.end()) {
                    set4.insert(a);
                }
            }
        }
        else{
            for (const auto& a : set1) {
                if (set2.find(a) != set2.end() && set3.find(a) != set3.end()) {
                    set4.insert(a);
                }
            }
        }
        return set4;
    }

    /* function: set_union
     *
     * compute union of two or more sets */
    std::set<size_t> set_union(const std::set<size_t>& set1, const std::set<size_t>& set2) const {
        std::set<size_t> result = set1;
        result.insert(set2.begin(), set2.end());
        return result;
    }
    std::set<size_t> set_union(const std::vector<std::set<size_t>>& sets) const {
        std::set<size_t> result;
        for (const auto& set: sets){
            set_merge(result,set);
        }
        return result;
    }

    /* function: edges_between
     *
     * return all plyaer i's edges from source to target */
    std::map<size_t, std::set<size_t>>  edges_between(const std::set<size_t>& source,
                        const std::set<size_t>& target,
                        const std::set<size_t>& players) const {
        std::map<size_t, std::set<size_t>> result_edges;
        /* include every player i edge from source to target */
        for (const auto& v : source){
            if (players.find(vert_id_.at(v)) != players.end()){
                for (const auto& u : edges_.at(v)){
                    if (target.find(u) != target.end()){
                        result_edges[v].insert(u);
                    }
                }
            }
        }
        return result_edges;
    }
    void edges_between(const std::set<size_t>& source,
                        const std::set<size_t>& target,
                        std::map<size_t, std::set<size_t>>& result_edges,
                        const std::set<size_t>& players) const {
        /* include every player i edge from source to target */
        for (const auto& v : source){
            if (players.find(vert_id_.at(v)) != players.end()){
                for (const auto& u : edges_.at(v)){
                    if (target.find(u) != target.end()){
                        result_edges[v].insert(u);
                    }
                }
            }
        }
    }
    std::map<size_t, std::set<size_t>> edges_between(const std::set<size_t>& source,
                        const std::set<size_t>& target,
                        std::set<size_t>& new_sources,
                        const std::set<size_t>& players) const {
        std::map<size_t, std::set<size_t>> result_edges;
        /* include every player i edge from source to target */
        for (const auto& v : source){
            if (players.find(vert_id_.at(v)) != players.end()){
                for (const auto& u : edges_.at(v)){
                    if (target.find(u) != target.end()){
                        result_edges[v].insert(u);
                        new_sources.insert(v);
                    }
                }
            }
        }
        return result_edges;
    }

    /* function: co_edges_between
     *
     * add all plyaer i's edges that are from source but not to target (and there is an edge from that source to target) */
    std::map<size_t, std::set<size_t>> co_edges_between(const std::set<size_t>& source,
                        const std::set<size_t>& target,
                        std::set<size_t>& new_sources,
                        std::map<size_t, std::set<size_t>>& result_edges,
                        const std::set<size_t>& players) const {
        /* include every player i edge from source but not to target when there is an edge from source to target */
        for (const auto& v : source){
            bool colive_source = false; /* if there is an edge from this source to target */
            std::set<size_t> colive_neighbours; /* all neighbours of this source that is not a target */
            if (players.find(vert_id_.at(v)) != players.end()){
                for (const auto& u : edges_.at(v)){
                    if (target.find(u) == target.end()){
                        colive_neighbours.insert(u);
                    }
                    else{
                        colive_source = true;
                    }
                }
                if(colive_source){ /* if there is an edge from source to target then add the edges to colive_neighbours */
                    set_merge(result_edges[v], colive_neighbours);
                    new_sources.insert(v);
                }
            }
        }
        return result_edges;
    }



    ///////////////////////////////////////////////////////////////
    ///Print functions
    ///////////////////////////////////////////////////////////////
    
    /* print game informations */
    int print_game(){
        if (labels_.empty()){
            std::cout << "Game constructed! #vertices:"<<n_vert_<<"  #edges:"<<n_edge_<<"  #colors:"<<max_color_+1<<"\n";
            return 0;
        }
        std::cout << "Game constructed! #vertices:"<<n_vert_-n_edge_/2<<"  #edges:"<<n_edge_/2<<"  #colors:"<<max_color_+1<<"\n";
        return 1;
    }

    /* function: print_set
     *
     * print out all elements of the set*/
    void print_set (const std::set<size_t>& set, const std::string& note = "set") const {
        std::cout << "\n" << note << ": ";
        for (auto u=set.begin(); u != set.end();){
            std::cout << *u;
            ++u;
            if (u != set.end())
                std::cout << ", ";
        }
        std::cout << "\n";
    }

}; /* close class definition */
} /* close namespace */

#endif
