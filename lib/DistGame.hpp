/*
 * Class: DistGame
 *
 * Class formalizing distributed parity games (both players have their own parity objective)
 */

#ifndef DISTGAME_HPP_
#define DISTGAME_HPP_

#include "Game.hpp"
#include <random>

// #include "FileHandler.hpp"


namespace mpa {
class DistGame: public Game {
public:
    /* all vertex ids */
    std::vector<std::map<size_t, size_t>> all_vert_id_;
    /* all maximum of colors */
    std::vector<size_t> all_max_color_;
    /* all_colors: the i-th vector represents i-th color set (for i=1,2)*/
    std::vector<std::map<size_t, size_t>> all_colors_;
    /* rounds of negotiation */
    size_t counter_;
    /* number of objectives (0th one for player 0, rest for player 1)*/
    size_t n_games_;
    
public:
    /* default constructor */
    DistGame(): Game() {
        n_games_ = 2;
        all_vert_id_ = std::vector<std::map<size_t, size_t>>(2,std::map<size_t, size_t>());
        all_max_color_ = std::vector<size_t>(2,0);
        all_colors_ = std::vector<std::map<size_t, size_t>>(2,std::map<size_t, size_t>());
        counter_ = 0;
    }

    /* copy a normal game */
    DistGame(const Game& other): Game(other){
        counter_ = 0;
        n_games_ = 1;
        all_max_color_.push_back(max_color_);
        all_colors_.push_back(colors_);
        all_vert_id_.push_back(vert_id_);
        all_vert_id_.push_back(vert_id_);
        for (const auto& v : vertices_){
            all_vert_id_[1][v] = 1 - all_vert_id_[0][v];
        }
    }

    ///////////////////////////////////////////////////////////////
    ///DistGame operators
    ///////////////////////////////////////////////////////////////

    /* find the parity game for n-th player */
    Game nthGame(const size_t& n) const{
        mpa::Game game(*this);
        if (n == 0){/* 0-th objective is player-0's; rest are player-1's */
            game.vert_id_ = all_vert_id_[0];
        }
        else{
            game.vert_id_ = all_vert_id_[1];
        }
        game.max_color_ = all_max_color_[n];
        game.colors_ = all_colors_[n];
        return game;
    }


    ///////////////////////////////////////////////////////////////
    ///Product of two games
    ///////////////////////////////////////////////////////////////

    /* compute product of two games */
    int product_games(const Game& game1, const Game& game2, const bool hoa = true) {
        *this = DistGame(); /* clear the game */
        n_games_ = 2; /* number of objective is two */
        counter_ = 0;

        std::map<std::pair<size_t,size_t>,size_t> product_verts; /* map of product of vertices to new vertieces */

        /* if one game is empty then return empty */
        if (game1.n_vert_ == 0 || game2.n_vert_ == 0){
            return 1;
        }
        /* update the aps : ids of aps in all 3 games */
        std::vector<std::vector<size_t>> common_aps; /* ids of common aps */
        std::vector<std::vector<size_t>> first_aps; /* only first ap-ids */
        std::vector<std::vector<size_t>> second_aps; /* only second ap-ids */
        std::map<size_t,size_t> second_aps_map; /* map for second aps to update controllable_aps */
        std::map<size_t,size_t> first_aps_map; /* map for first aps to update controllable_aps */
        std::set<size_t> common_in_second; /* id of commons aps in second game -- need to compute second_aps */
        for (size_t i = 0; i < game1.ap_id_.size(); i++){
            bool common = false; /* if i-th vertex is common in both game */
            for (size_t j = 0; j < game2.ap_id_.size(); j++){
                if (game1.ap_id_.at(i) == game2.ap_id_.at(j)){ /* i-th vertex is common */
                    common_aps.push_back(std::vector<size_t>{ap_id_.size(),i,j}); /* add it in common_aps vector with ids in all games */
                    second_aps_map[j] = ap_id_.size();
                    common_in_second.insert(j); 
                    common = true; /* set common to true */
                    break;
                }
            }
            if (!common){/* if not common then add it to first_aps */
                first_aps.push_back(std::vector<size_t>{ap_id_.size(),i});
            }
            /* add this ap and its label to new game */
            first_aps_map[i] = ap_id_.size();
            ap_id_[ap_id_.size()] = game1.ap_id_.at(i);
        }
        for (size_t j = 0; j < game2.ap_id_.size(); j++){ /* iterate over all aps of 2nd game */
            if (common_in_second.find(j) == common_in_second.end()){ /* if it is not common then add it to second_aps and new game */
                second_aps.push_back(std::vector<size_t>{ap_id_.size(),j});
                ap_id_[ap_id_.size()] = game2.ap_id_.at(j);
                second_aps_map[j] = ap_id_.size();
            }
        }

        /* update the controllable_aps */
        for (const size_t& a : game1.controllable_ap_){
            controllable_ap_.insert(first_aps_map.at(a));
        }
        for (const size_t& a : game2.controllable_ap_){
            controllable_ap_.insert(second_aps_map.at(a));
        }
        

        /* initial vertex should be of same player in both games */
        if (game1.vert_id_.at(game1.init_vert_) != game2.vert_id_.at(game2.init_vert_)){
            std::cerr << "Error: vertex ids of both games are not same!\n";
        }
        /* insert initial vertex and update all variables */
        n_vert_ = 1;
        size_t org_vert = 1; /* counter for normal vertices */
        size_t edge_vert = (game1.n_vert_-game1.n_edge_/2)*(game2.n_vert_-game2.n_edge_/2); /* counter for edge-vertices */
        vertices_.insert(0); 
        init_vert_ = 0;
        all_vert_id_[0][0] = game1.vert_id_.at(game1.init_vert_);
        all_vert_id_[1][0] = 1-all_vert_id_[0][0];
        all_colors_[0][0] = game1.colors_.at(game1.init_vert_);
        all_colors_[1][0] = game2.colors_.at(game2.init_vert_);
        product_verts.insert({std::make_pair(game1.init_vert_,game2.init_vert_),0});

        /* maintain a stack to explore new (product) vertices */
        std::stack<std::vector<size_t>> stack_list;
        /* initialize the stack with initial vertex */
        stack_list.push(std::vector<size_t>{0,game1.init_vert_,game2.init_vert_});

        /* explore until the stack list is empty */
        while (!stack_list.empty()){
            std::vector<size_t> curr = stack_list.top(); /* current pair of vertices */
            stack_list.pop(); /* pop the top element from stack list */
            for (const auto& u : game1.edges_.at(curr[1])){ /* for each edge-neighbour of 1st vertex */
                for (const auto& v : game2.edges_.at(curr[2])){ /* for each edge-neihbour of 2nd vertex */
                    bool valid = true; /* if product of these two edges is possible */
                    std::vector<size_t> temp_common(ap_id_.size(),2); /* needed if possible, temporarily store the ids of this new edge */
                    for (const auto& ap : common_aps){ /* first go through common aps */
                        if (game1.labels_.at(u)[ap[1]]+game2.labels_.at(v)[ap[2]] == 1){
                            valid = false; /* one of the id of this ap is 1 and other is 0, so not valid product */
                            break;
                        }
                        /* if this ap is consistent in both edges then add its label to temp_common */
                        else if (game1.labels_.at(u)[ap[1]] == 0 || game2.labels_.at(v)[ap[2]] == 0){
                            temp_common[ap[0]] = 0;
                        }
                        else if (game1.labels_.at(u)[ap[1]] == 1 || game2.labels_.at(v)[ap[2]] == 1){
                            temp_common[ap[0]] = 1;
                        }
                    }
                    if (valid){ /* if the product-edge is valid */
                        /* update its temp_common to all ids */
                        for (const auto& ap : first_aps){
                            temp_common[ap[0]] = game1.labels_.at(u)[ap[1]];
                        }
                        for (const auto& ap : second_aps){
                            temp_common[ap[0]] = game2.labels_.at(v)[ap[1]];
                        }
                        size_t newId = org_vert; 
                        if (hoa){
                            newId = edge_vert;/* new state for product-edges : use edge_vert counter */
                            edge_vert += 1;
                        }
                        else{
                            org_vert += 1;
                        }
                        /* update all variables for new game */
                        n_vert_ += 1;
                        n_edge_ += 2;
                        vertices_.insert(newId);
                        all_vert_id_[0][newId] = 2;
                        all_vert_id_[1][newId] = 2;
                        edges_[curr[0]].insert(newId);
                        labels_[newId] = temp_common;
                        all_colors_[0][newId] = game1.colors_.at(u);
                        all_colors_[1][newId] = game2.colors_.at(v);

                        /* product of end-states of these edges */
                        std::pair<size_t,size_t> new_succ = std::make_pair(*game1.edges_.at(u).begin(),*game2.edges_.at(v).begin());
                        if (product_verts.find(new_succ) != product_verts.end()){ /* if the product of end-states of edges are already present add the edge */
                            edges_[newId].insert(product_verts.at(new_succ));
                            pre_edges_[product_verts.at(new_succ)].insert(newId);
                        }
                        else{ /* if the product of end-states of edges are not present */
                            size_t succId = org_vert; /* new state for the product-successor-states : use org_vert counter */
                            /* update all variables in the new game */
                            org_vert += 1;
                            n_vert_ += 1;
                            vertices_.insert(succId);
                            all_vert_id_[0][succId] = game1.vert_id_.at(new_succ.first);
                            all_vert_id_[1][succId] = 1-all_vert_id_[0][succId];
                            edges_[newId].insert(succId);
                            pre_edges_[succId].insert(newId);
                            all_colors_[0][succId] = game1.colors_.at(new_succ.first);
                            all_colors_[1][succId] = game2.colors_.at(new_succ.second);
                            product_verts[new_succ] = succId;
                            stack_list.push(std::vector<size_t>{succId,new_succ.first,new_succ.second});
                        }
                    }
                }
            }
        }
        /* update max_color */
        all_max_color_[0] = max_col(all_colors_[0]); 
        all_max_color_[1] = max_col(all_colors_[1]); 
        return 1;
    }


    ///////////////////////////////////////////////////////////////
    ///Game to DistGame with random sets of colors
    ///////////////////////////////////////////////////////////////

    /* generate and add n_games sets of colors to get a DistGame */
    void randDistGame(const size_t& n_games, const size_t& max_col, const bool clear = true){
        if (clear){
            /* clear the colors */
            all_max_color_.clear();
            all_colors_.clear();
            n_games_ = 0;
        }
        /* update number of games */
        n_games_ += n_games;
        /* update max_color */
        max_color_ = max_col;
        /* add two set of colors and update max_color */
        for (size_t i = 0 ; i < n_games; i++){
            all_max_color_.push_back(max_col);
            all_colors_.push_back(random_colors(max_col));
        }    
    }

    /* generate a set of random colors <= max_col */
    std::map<size_t, size_t> random_colors(const size_t& max_col){
        std::map<size_t, size_t> colors;
        
        /* a vector of all vertices */
        std::vector<size_t> vertices(vertices_.begin(), vertices_.end());
        
        size_t min_num = n_vert_/(2*max_col); /* each color has atleast this many vertices */
        std::vector<size_t> remaining_vert = vertices; /* remaining vertices to be colored (initially all) */
        size_t remaining_num = n_vert_-1; /* number of remaining vertices - 1 (initially n_vert_-1) */
            
        /* randomly choose min_num vertices for each color */
        for (size_t col = 0; col <= max_col; col++){
            for (size_t i = 0; i < min_num; i++){
                size_t rand_index = random_num(remaining_num); /* generate random index */
                size_t vertex = remaining_vert[rand_index]; /* vertex at rand_index in remaining vertices */
                
                colors.insert(std::make_pair(vertex, col)); /* set color of that vertex as col */
                remaining_vert.erase(remaining_vert.begin() + rand_index); /* remove that vertex from remaining vertices */
                remaining_num -= 1; /* update remaining number of vertices */
            }
        }
        /* randomly choose color for remaining vertices */
        for (const auto& vertex : remaining_vert){
            size_t rand_col = random_num(max_col); /* generate random color */
            colors.insert(std::make_pair(vertex, rand_col)); /* set color of vertex to rand_col */
        }
        return colors;
    }

    /* random number from [0,n] */
    size_t random_num(const size_t& max, const size_t& min = 0, const size_t& seed = 42){
        static std::mt19937 gen(seed); // fixed seed for reproducibility
        std::uniform_int_distribution<> distr(min, max); // define the range

        return distr(gen); // generate numbers
    }


    ///////////////////////////////////////////////////////////////
    ///Negotiate between two players to compute AG contract
    ///////////////////////////////////////////////////////////////
    
    /* compute the AG contract by negotiation */
    std::pair<std::set<size_t>, std::set<size_t>> find_AG_contract(std::vector<Template>& assumps, std::vector<Template>& strats) {
        /* clear the template */
        assumps = std::vector<Template> (n_games_, Template {});
        strats = std::vector<Template> (n_games_, Template {});
        
        /* solve the games without changing anything in the original game */
        DistGame DistGame_copy = *this; /* copy of the multi-game */
        auto winning_region = DistGame_copy.recursive_negotiate(assumps, strats);
        /* update losing region and counter */
        winning_region.second = set_complement(winning_region.first);
        counter_ = DistGame_copy.counter_;
        /* unsafe edges are the player i's edges from winning region to losing region */
        strats[0].unsafe_edges_ = nthGame(0).edges_between(winning_region.first, winning_region.second, {V0});
        strats[1].unsafe_edges_ = nthGame(0).edges_between(winning_region.first, winning_region.second, {V1});
        /* return winning region */
        return winning_region;
    }

    /* recursively negotiate between two players */
    std::pair<std::set<size_t>, std::set<size_t>> recursive_negotiate(std::vector<Template>& assumps, std::vector<Template>& strats, std::set<size_t> colive_vertices = {}) {
        std::pair<std::set<size_t>, std::set<size_t>> winning_region; /* overall winning region of the games */
        std::vector<std::set<size_t>> losing_regions(n_games_, std::set<size_t> {}); /* losing region of each game */
        std::vector<std::set<size_t>> i_colive_vertices(n_games_,std::set<size_t>{}); /* vertices not to visit inf often in each game */
        // std::set<size_t> colive_vertices; /* union of both i_colive vertices */
        while (true){/* iterate until there is no need to solve any game again */
            /* compute template for every game */
            #pragma omp parallel
            #pragma omp for 
            for (size_t i = 0; i < n_games_; i++){
                size_t colive_color = max_odd(all_max_color_[i]); /* minimum odd color >= max color */
                for (const size_t& v : colive_vertices){/* set color of all colive vertices colive_color */
                    all_colors_[i].at(v) = colive_color;
                    all_max_color_[i] = colive_color;
                }
                Game game = nthGame(i); /* i-th game */
                
                /* compute template and winning region is intersection with the current one */
                losing_regions[i] = game.find_colive_cond_live_parity(assumps[i], strats[i], i_colive_vertices[i]).second;
            }
            
            /* compute the overall winning region */
            winning_region.second = nthGame(0).solve_reachability_game(set_union(losing_regions),{}).first; 
            winning_region.first = set_complement(winning_region.second);

            /* remove colive edges from losing regions */
            for (size_t i = 0; i < n_games_; i++){
                map_remove_keys(assumps[i].colive_edges_, winning_region.second);
                map_remove_keys(strats[i].colive_edges_, winning_region.second);
            }

            /* merge assumptions and strategies for each player */
            Template all_temps;
            std::vector<Template> temp_set = assumps;
            temp_set.insert(temp_set.end(), strats.begin(), strats.end());
            all_temps.merge_live_colive(temp_set);
                

            /* check if there is no conflicts by live groups and colive edges */
            if (!isImplemntable(all_temps, winning_region)){/* if there is no conflicts  */
                for (size_t i =2; i < n_games_; i++){
                    assumps[1].merge(assumps[i]);
                    strats[1].merge(strats[i]);
                }
                for (size_t i = 0; i < 2; i++){
                    /* remove the unsafe edges and edges from losing region from colive edge set */
                    map_remove_values(assumps[i].colive_edges_, winning_region.second);
                    map_remove_values(strats[i].colive_edges_, winning_region.second);
                    /* remove the restrictions on losing regions and colive edges from live group (as there is other choice from that source) */
                    for (auto& live_groups : assumps[i].cond_live_groups_){
                        for (auto& live_group : live_groups){
                            map_remove_vertices(live_group, winning_region.second);
                            map_remove(live_group, strats[1-i].colive_edges_);
                            map_remove(live_group, assumps[i].colive_edges_);
                        }
                    }
                    for (auto& live_groups : strats[i].cond_live_groups_){
                        for (auto& live_group : live_groups){
                            map_remove_vertices(live_group, winning_region.second);
                            map_remove(live_group, assumps[1-i].colive_edges_);
                            map_remove(live_group, strats[i].colive_edges_);
                        }
                    }
                }
                remove_vertices(winning_region.second);
                break; /* break the loop as there is no more conflict */
            }
            
            /* print to analyze the results */
            counter_ += 1;
            // if (counter_ <= 2){
            //     std::cout << "couter:"<<counter_<< "  colive:"<<colive_vertices.size()<<"  losing:"<<winning_region.second.size()<<"\n";
            // }
            
            /* remove losing region from everywhere and update the colive_vertices */
            remove_vertices(winning_region.second);
            colive_vertices = set_difference(set_union(i_colive_vertices), winning_region.second);
            
            /* clear the template and i_colive_vertices */
            for (size_t i = 0 ; i < n_games_; i++){
                assumps[i].clear();
                strats[i].clear();
                i_colive_vertices[i].clear();
            }
        }
        /* print to analyze the results */
        // std::cout << "couter:"<<counter_<< "  colive:"<<colive_vertices.size()<<"  winning:"<<winning_region.first.size()<<"\n"; 
        return winning_region; /* return winning region */
    }

    /* check if there is any conflict */
    bool isImplemntable(Template temp, const std::pair<std::set<size_t>, std::set<size_t>>& winning_region) const {
        if (conflict_colive(temp.colive_edges_) || conflict_live_colive(temp.live_groups_, temp.colive_edges_) || conflict_unsafe(temp.live_groups_, temp.colive_edges_, winning_region)){
                return true;
        }
        return false; /* return false if there is no conflict */
    }
    
    /* Check if the unsafe edges create some conflict */
    bool conflict_unsafe(std::vector<std::map<size_t, std::set<size_t>>>& live_group_set,
                            std::map<size_t, std::set<size_t>>& colive_edge_set,
                            const std::pair<std::set<size_t>, std::set<size_t>>& winning_region) const {

        for (const auto& v : winning_region.first){
            if (check_set_inclusion(edges_.at(v), set_union(winning_region.second, colive_edge_set[v]))){
                return true; /* return true if there is a conflict */
            }
        }                
        for (auto& live_group : live_group_set){ /* iterate over all live groups to compute live_unsafe_region */
            for (const auto& v : winning_region.first){
                if (!live_group[v].empty() &&  check_set_inclusion(live_group[v], set_union(winning_region.second, colive_edge_set[v]))){
                    return true; /* return true if there is a conflict */
                }
            }
        }
        return false; /* return false if there is no conflict */
    }

    /* check conflict when the union of all colive edge set contains all edges of some vertex */
    bool conflict_colive(std::map<size_t, std::set<size_t>>& colive_edge_set) const {
        for (auto it = colive_edge_set.begin(); it != colive_edge_set.end(); it++){
            auto v = it ->first;
            if (check_set_inclusion(edges_.at(v), it->second)){
                /* if all edges of a source is colive then there is a conflict */
                return true;
            }
        }
        return false; /* return false if there is no conflict */
    }

    /* solve the conflict when the intersection of colive edges and live groups is non-empty */
    bool conflict_live_colive(std::vector<std::map<size_t, std::set<size_t>>>& live_group_set,
                                const std::map<size_t, std::set<size_t>>& colive_edge_set) const {
        for (auto it = colive_edge_set.begin(); it != colive_edge_set.end(); it++){
            auto v = it ->first;
            for (auto& live_group : live_group_set){ /* iterate over all live groups */
                if (!live_group[v].empty() && check_set_inclusion(live_group[v], it->second)){
                    /* for any conflict source, if there is no other choice in the live group
                    then there is a conflict */
                    return true;
                }
            }
        }
        return false; /* return false if there is no conflict */
    }



    ///////////////////////////////////////////////////////////////
    ///Solve distributed games using AG contract
    ///////////////////////////////////////////////////////////////

    /* solve distributed parity game  */
    std::pair<std::set<size_t>, std::set<size_t>> solve_dist_game() {
        /* initialize the strategy template */
        std::vector<mpa::Template> assumps;
        std::vector<mpa::Template> strats;
        return find_AG_contract(assumps,strats);
    }

    /* solve cooperative parity game by converting it to distGame */std::pair<std::set<size_t>, std::set<size_t>> solve_coop_game_with_negotiate() {
        DistGame copy(*this);
        copy.parityToDistGame();
        counter_ = copy.counter_;
        return copy.solve_dist_game();
    }

    /* convert the parity game to multiple small games */
    void parityToDistGame(){
        /* sanity check: number of games should be 1 */
        if (colors_.size() == 0){
            std::cerr << "[ERROR] parityToMultigame: one normal game is needed.\n";
        }
        /* clear the multigame parts except gamegraph (all_colors, n_games) */
        all_colors_.clear();
        all_max_color_.clear();
        n_games_ = 0;

        /* for every odd color create a new color_vector for new game */
        for (size_t odd_col = 1; odd_col <= max_color_; odd_col+= 2){
            std::map<size_t, size_t> colors; /* new color vector */
            size_t max_color=2;
            for (const auto& vertex : vertices_){/* loop through each vertex */
                if (colors_[vertex] < odd_col){ 
                    colors.insert(std::make_pair(vertex,0)); /* new color of vertices with color < odd_col is 0 */
                }
                else if (colors_[vertex] == odd_col){
                    colors.insert(std::make_pair(vertex,1)); /* new color of vertices with color = odd_col is 1 */
                }
                else{
                    if (colors_[vertex]%2 == 0){
                        colors.insert(std::make_pair(vertex,2)); /* new color of vertices with even color > odd_col is 2 */
                    }
                    else{
                        colors.insert(std::make_pair(vertex,1)); /* new color of vertices with odd color > odd_col is 1 */
                    }
                } 
            }
            all_colors_.push_back(colors);
            all_max_color_.push_back(max_color);
            n_games_ += 1;
        }
    }




    ///////////////////////////////////////////////////////////////
    ///Incremental distributed synthesis
    ///////////////////////////////////////////////////////////////

    /* just solve the game incrementally */
    int solve_dist_game_incr() {
        std::vector<mpa::Template> assumps;
        std::vector<mpa::Template> strats;
        auto winning_region = find_AG_contact_incr(assumps, strats);
        return 0;
    }

    /* solve distGame incrementally adding specs one by one */
    std::pair<std::set<size_t>, std::set<size_t>> find_AG_contact_incr(std::vector<mpa::Template>& assumps, std::vector<mpa::Template>& strats) {

        /* copy the game to modify it */
        DistGame copy(*this);
        
        /* solve it once for 1st two objectives */
        copy.n_games_ = 2;
        /* initialize the templates */
        assumps = std::vector<Template> (n_games_, Template {});
        strats = std::vector<Template> (n_games_, Template {});
        auto winning_region = copy.recursive_negotiate(assumps,strats);

        /* solve for other objectives incrementally */
        for (size_t i = 3; i <= n_games_; i++){
            copy.n_games_ = i;
            assumps.push_back(Template {});
            strats.push_back(Template {});
            winning_region = copy.recursive_negotiate_incr(assumps,strats);
        }

        counter_ = copy.counter_;
        return winning_region;
    }

    /* recursively negotiate between two players with a newly added objective for player 1 */
    std::pair<std::set<size_t>, std::set<size_t>> recursive_negotiate_incr(std::vector<Template>& assumps, std::vector<Template>& strats) {
        /* merge all colive edges */
        std::map<size_t, std::set<size_t>> all_colive_edges;
        assumps[0].edge_merge(all_colive_edges,assumps[0].colive_edges_);
        assumps[0].edge_merge(all_colive_edges,assumps[1].colive_edges_);
        assumps[0].edge_merge(all_colive_edges,strats[0].colive_edges_);
        assumps[0].edge_merge(all_colive_edges,strats[1].colive_edges_);
    

        /* collect all colive vertices */
        std::set<size_t> colive_vertices;
        for (const auto& pair : all_colive_edges){
            set_merge(colive_vertices,pair.second);
        }
        
        /* compute template for last game */
        size_t colive_color = max_odd(all_max_color_[n_games_-1]); /* minimum odd color >= max color */
        colive_vertices = set_intersection(vertices_,colive_vertices);
        for (const size_t& v : colive_vertices){/* set color of all colive vertices colive_color */
            all_colors_[n_games_-1].at(v) = colive_color;
            all_max_color_[n_games_-1] = colive_color;
        }
        
        Game game = nthGame(n_games_-1); /* i-th game */
        
        /* compute template and winning region is intersection with the current one */
        auto winning_region = game.find_colive_cond_live_parity(assumps[2], strats[2], colive_vertices);

        /* remove colive edges from losing regions */
        for (size_t i = 0; i < 3; i++){
            map_remove_keys(assumps[i].colive_edges_, winning_region.second);
            map_remove_keys(strats[i].colive_edges_, winning_region.second);
        }

        /* merge assumptions and strategies for each player */
        Template all_temps;
        all_temps.merge_live_colive(std::vector<Template> {assumps[0],assumps[1],assumps[2],strats[0],strats[1],strats[2]});
        

        /* check if there is no conflicts by live groups and colive edges */
        if (!isImplemntable(all_temps, winning_region)){/* if there is no conflicts  */
            for (size_t i =2; i < 3; i++){
                assumps[1].merge(assumps[i]);
                strats[1].merge(strats[i]);
            }
            for (size_t i = 0; i < 2; i++){
                /* remove the unsafe edges and edges from losing region from colive edge set */
                map_remove_values(assumps[i].colive_edges_, winning_region.second);
                map_remove_values(strats[i].colive_edges_, winning_region.second);
                /* remove the restrictions on losing regions and colive edges from live group (as there is other choice from that source) */
                for (auto& live_groups : assumps[i].cond_live_groups_){
                    for (auto& live_group : live_groups){
                        map_remove_vertices(live_group, winning_region.second);
                        map_remove(live_group, strats[1-i].colive_edges_);
                        map_remove(live_group, assumps[i].colive_edges_);
                    }
                }
                for (auto& live_groups : strats[i].cond_live_groups_){
                    for (auto& live_group : live_groups){
                        map_remove_vertices(live_group, winning_region.second);
                        map_remove(live_group, assumps[1-i].colive_edges_);
                        map_remove(live_group, strats[i].colive_edges_);
                    }
                }
            }
            remove_vertices(winning_region.second);
            return winning_region;
        }

        /* print to analyze the results */
        counter_ += 1;
        // if (counter_ <= 2){
        //     std::cout << "couter:"<<counter_<< "  colive:"<<colive_vertices.size()<<"  losing:"<<winning_region.second.size()<<"\n";
        // }

        /* remove losing region from everywhere and update the colive_vertices */
        remove_vertices(winning_region.second);
        colive_vertices = set_difference(colive_vertices, winning_region.second);

        /* clear the template and i_colive_vertices */
        for (size_t i = 0 ; i < n_games_; i++){
            assumps[i].clear();
            strats[i].clear();
        }
        
        /* print to analyze the results */
        // std::cout << "couter:"<<counter_<< "  colive:"<<colive_vertices.size()<<"  winning:"<<winning_region.first.size()<<"\n"; 
        return recursive_negotiate(assumps,strats,colive_vertices); /* return winning region */
    }





    ///////////////////////////////////////////////////////////////
    /// Filter all edge-states out (needed for HOA formatted games)
    ///////////////////////////////////////////////////////////////

    /* filter out edge-states from a winning region, assumption and strategy template */
    int filter_out_edge_states(std::pair<std::set<size_t>, std::set<size_t>>& winning_region, std::vector<Template>& assumps, std::vector<Template>& strats) const {
        if (labels_.empty()){
            return 0;
        }
        std::set<size_t> org_vertices;
        for (const auto& v : vertices_){
            if (all_vert_id_[0].at(v) != 2){
                org_vertices.insert(v);
            }
        }
        winning_region.first = set_intersection(winning_region.first, org_vertices);
        winning_region.second = set_intersection(winning_region.second, org_vertices);
        
        for (size_t i = 0; i < assumps.size(); i++){
            filter_templates(assumps[i], org_vertices);
            filter_templates(strats[i], org_vertices);
        }
        return 1;
    }
    
    ///////////////////////////////////////////////////////////////
    ///Basic functions
    ///////////////////////////////////////////////////////////////

    /* function: max_odd
     *
     * return minimum odd color that is greater than or equal to max_color */
    size_t max_odd(const std::map<size_t, size_t>& colors) const {
        size_t odd_col = 1;
        for (auto const & col : colors){
            if (col.second%2 == 1 && col.second > odd_col)
                odd_col = col.second;
            else if (col.second%2 == 0 && col.second > odd_col)
                odd_col = col.second+1;
        }
        return odd_col;
    }
    size_t max_odd(const size_t& max_color) const {
        if (max_color%2 == 1){
            return max_color;
        }
        else{
            return max_color+1;
        }
    }

    /* function: vec_erase_duplicate
     *
     * erase duplicates in a vector */
    template<class T>
    void vec_erase_duplicate(std::vector<T>& vec) const {
        std::vector<T> new_vec;
        for (const auto& element : vec){
            if (std::find(new_vec.begin(), new_vec.end(), element) == new_vec.end()){
                new_vec.push_back(element);
            }
        }
        vec = new_vec;
    }


    /* function: remove_vertices
     *
     * remove vertices from the game */
    
    void remove_vertices(const std::set<size_t>& set) {
        vertices_= set_complement(set);
        n_vert_ = vertices_.size();
        
        map_remove_keys(vert_id_, set);
        
        for (auto& colors: all_colors_){
            map_remove_keys(colors, set);
            max_color_ = std::max(max_color_, max_col(colors));
        }
        colors_ = all_colors_[0];

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
    
    DistGame subgame(const std::set<size_t>& set) const{
        DistGame game = *this;
        game.n_vert_ = set.size();
        game.vertices_= set;

        std::set<size_t> complement = set_complement(set);
        map_remove_keys(game.vert_id_, complement);
        
        for (auto& colors: game.all_colors_){
            map_remove_keys(colors, set);
            game.max_color_ = std::max(game.max_color_, max_col(colors));
        }
        game.colors_ = game.all_colors_[0];

        map_remove_keys(game.edges_, set);
        game.n_edge_ = 0;
        for (const auto& v : game.vertices_){
            game.edges_.at(v) = set_difference(game.edges_.at(v),complement);
            game.n_edge_ += game.edges_.at(v).size();
        }
        /* sanity check */
        // game.valid_game();

        return game;
    }

    /* function: map_remove
     *
     * remove set of values from one map (set of edges) */
    
    void map_remove(std::map<size_t, std::set<size_t>>& map2, const std::map<size_t, std::set<size_t>>& map1) const {
        for (auto v = map1.begin(); v != map1.end(); v++){
            map2[v->first] = set_difference(map2[v->first], v->second);
        }
    }

    /* print game informations */
    int print_game(){
        if (labels_.empty()){
            std::cout << "Game constructed! #vertices:"<<n_vert_<<"  #edges:"<<n_edge_;
            if (n_games_ == 2){
                std::cout << "  #colors0:"<<all_max_color_[0]+1<<"  #colors1:"<<all_max_color_[1]+1<<"\n";
                return 0;
            }
            size_t max_color = 0;
            for (size_t i = 1; i<n_games_; i++){
                max_color = std::max(max_color,all_max_color_[i]);
            }
            std::cout << "  #games:"<<n_games_<< "  #colors0:"<<all_max_color_[0]+1<<"  #colors1:"<<max_color+1<<"\n";
            return 1;
        }
        std::cout << "Game constructed! #vertices:"<<n_vert_-n_edge_/2<<"  #edges:"<<n_edge_/2<<"  #colors0:"<<all_max_color_[0]+1<<"  #colors1:"<<all_max_color_[1]+1<<"\n";
        return 1;
    }
    
}; /* close class definition */
} /* close namespace */

#endif