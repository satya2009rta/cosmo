/*
 * Class: Maze
 *
 * Class formalizing a maze with two robots
 */

#ifndef MAZE_HPP_
#define MAZE_HPP_

#include "DistGame.hpp"
#include <random>
#include "FileHandler.hpp"

// // /* agnes files */
// #include "Component.hpp"
// #include "SafetyAutomaton.hpp"
// #include "Monitor.hpp"
// //#include "LivenessGame.hpp"
// #include "SafetyGame.hpp"
// #include "Spoilers.hpp"
// #include "Negotiate.hpp"
// #include "TicToc.hpp"



namespace mpa {
class Location{
public:
    /* x-coordinate */
    size_t x_;
    /* y-coordinate */
    size_t y_;

public:
    /* default constructor */
    Location(){
        x_ = 1;
        y_ = 1;
    }

    /* main constructor */
    Location(size_t a, size_t b){
        x_ = a;
        y_ = b;
    }

    /* if inside the boundary */
    bool inBoundary(const Location boundary) const{
        if (x_ <= boundary.x_ && y_ <= boundary.y_){
            return true;
        }
        return false;
    }

    /* all directions */
    Location up() const{
        Location new_loc = *this;
        new_loc.y_ += 1;
        return new_loc;
    }
    Location down() const{
        Location new_loc = *this;
        new_loc.y_ -= 1;
        return new_loc;
    }
    Location right() const{
        Location new_loc = *this;
        new_loc.x_ += 1;
        return new_loc;
    }
    Location left() const{
        Location new_loc = *this;
        new_loc.x_ -= 1;
        return new_loc;
    }

    /* binary operators */
    bool operator==(const Location& k) const{
        return x_ == k.x_ && y_ == k.y_;
    }
    bool operator!=(const Location& k) const{
        return (x_ != k.x_) || (y_ != k.y_);
    }
    // bool operator<(const Location& k) const{
    //     return (x_ < k.x_) || (x_ == k.x_ && y_ < k.y_);
    // }
    // bool operator>(const Location& k) const{
    //     return (x_ > k.x_) || (x_ == k.x_ && y_ > k.y_);
    // }
    // bool operator<=(const Location& k) const{
    //     return (x_ < k.x_) || (x_ == k.x_ && y_ <= k.y_);
    // }
    // bool operator>=(const Location& k) const{
    //     return (x_ > k.x_) || (x_ == k.x_ && y_ >= k.y_);
    // }

    /* print the location */
    void print(const std::string line = "") const{
        std::cout << "("<<x_<<","<<y_<<")"<<line;
    }
    std::string print2string(const std::string line = "") const{
        std::string output = "(" + std::to_string(x_) +","+std::to_string(y_)+")"+line;
        return output;
    }

    /* convert location from/to index using boundary (as default location) */
    size_t to_index(const Location loc) const{
        return (x_*(loc.x_-1)+loc.y_-1);
    }
    Location from_index(const size_t i) const{
        /* sanity check: x_ should be at least y_n */
        if (x_ < y_){
            std::cerr << "Error: from index!\n";
        }
        Location loc;
        loc.y_ = i%x_ +1;
        loc.x_ = (i-loc.y_+1)/x_ +1;
        return loc;
    }

}; /* close Location class definition */


class Maze: public DistGame {
public:
    /* boundary */
    Location boundary_;
    /* explicit directed walls for Robot i */
    std::vector<std::vector<std::vector<Location>>> walls_;
    /* initial location of each robot*/
    std::vector<Location> init_state_; 
    /* location of each robot represented by a vertex */
    std::map<size_t,std::vector<Location>> states_; 
    /* number of walls and corridors */
    std::size_t n_walls_, n_corridors_;
    
public:
    /* default constructor */
    Maze(): DistGame() {
        boundary_ = Location(1,1);
        n_walls_ = 0;
        n_corridors_ = 0;
    }


    ///////////////////////////////////////////////////////////////
    /// Basic constructor from explicit walls and boundary
    ///////////////////////////////////////////////////////////////

    /* construct from boundary and explicit (directed) walls */
    void mazeSimple(const Location boundary, 
        const std::vector<Location> init_state, 
        const std::vector<Location> targets, 
        const std::vector<std::vector<std::vector<Location>>> walls,
        const bool standard_distance = true){
        boundary_ = boundary;
        walls_ = walls;
        init_state_ = init_state;
        ap_id_ = {{0, "s2"},
              {1, "u2"},
              {2, "d2"},
              {3, "r2"},
              {4, "l2"},
              {5, "ll2"},
              {6, "lr2"},
              {7, "ul2"},
              {8, "ur2"},
              {9, "s1"},
              {10, "u1"},
              {11, "d1"},
              {12, "r1"},
              {13, "l1"},
              {14, "ll1"},
              {15, "lr1"},
              {16, "ul1"},
              {17, "ur1"}};
        controllable_ap_ = {0,1,2,3,4,5,6,7};
        /* initialize initial vertex with init_state_ */
        n_vert_ = 1;
        n_edge_ = 0;
        vertices_.insert(0);
        init_vert_ = 0;
        vert_id_[0]= 1;
        all_vert_id_[0][0] = vert_id_[0];
        all_vert_id_[1][0] = 1-vert_id_[0];
        states_[0] = init_state_;
        size_t org_vert = 1; /* counter for normal vertices */
        size_t edge_vert = 2*(boundary_.x_*boundary_.y_)*(boundary_.x_*boundary_.y_); /* counter for edge-vertices */

        /* stack to explore new vertices */
        std::stack<size_t> stack_list;
        stack_list.push(0);

        /* explore until stack_list is empty */
        while (!stack_list.empty()){
            size_t curr = stack_list.top(); /* current vertex to explore */
            stack_list.pop(); /* pop the top element */

            size_t id = vert_id_[curr]; /* vertex id of current vertex */
            auto state = states_[curr][id]; /* current state of id-th robot */

            /* if id-th robot stayed there */
            auto next_states = std::vector<Location>{2,Location()};
            next_states[id] = state;
            next_states[1-id]=states_[curr][1-id];
            new_direction(curr,next_states,0,stack_list,org_vert,edge_vert);

            /* if id-th robot can go up */
            if (safe_distance(state.up(),next_states[1-id],standard_distance) && state.up().inBoundary(boundary_) && !isWall(state,state.up(),walls[id])){
                /* next states when id-th robot go up */
                next_states[id] = state.up();
                /* add the new direction */
                new_direction(curr,next_states,1,stack_list,org_vert,edge_vert);
            }

            /* if id-th robot can go down */
            if (state.y_ > 1 && safe_distance(state.down(),next_states[1-id],standard_distance) && state.down().inBoundary(boundary_) && !isWall(state,state.down(),walls[id])){
                /* next states when id-th robot go up */
                next_states[id] = state.down();
                /* add the new direction */
                new_direction(curr,next_states,2,stack_list,org_vert,edge_vert);
            }

            /* if id-th robot can go right */
            if (safe_distance(state.right(),next_states[1-id],standard_distance) && state.right().inBoundary(boundary_) && !isWall(state,state.right(),walls[id])){
                /* next states when id-th robot go up */
                next_states[id] = state.right();
                /* add the new direction */
                new_direction(curr,next_states,3,stack_list,org_vert,edge_vert);
            }

            /* if id-th robot can go left */
            if (state.x_ > 1 && safe_distance(state.left(),next_states[1-id],standard_distance) && state.left().inBoundary(boundary_) && !isWall(state,state.left(),walls[id])){
                /* next states when id-th robot go up */
                next_states[id] = state.left();
                /* add the new direction */
                new_direction(curr,next_states,4,stack_list,org_vert,edge_vert);
            }
        }

        /* update the colors */
        all_max_color_[0] = 2;
        all_max_color_[1] = 2;
        for (auto v : vertices_){
            if (vert_id_[v]!=2){ /* colors are only on edge-vertices */
                all_colors_[0][v]=1;
                all_colors_[1][v]=1;
            }
            else{/* if succ-state of this edge-vertex is in i-th target, then its ith color is 2 else it is 1 */
                auto succ = *edges_[v].begin();
                if (states_[succ][0] == targets[0]){
                    all_colors_[0][v] = 2;
                }
                else{
                    all_colors_[0][v] = 1;
                }
                if (states_[succ][1] == targets[1]){
                    all_colors_[1][v] = 2;
                }
                else{
                    all_colors_[1][v] = 1;
                }
            }
        }

        /* update the state_names_ to locations of the robots */
        for (auto pair : states_){
            state_names_[pair.first] = "R1:"+pair.second[1].print2string()+" & R2:"+pair.second[0].print2string();
        }

    }

    /* check if both robot's locations are in safe distance */
    bool safe_distance(const Location A, const Location B, const bool standard = true){
        if (standard){ /* if we follow the standard rule, both location should not be equal */
            if (A == B){
                return false;
            }
            return true;
        }
        else{ /* else we need to check if the locations are blocking a direction */
            if(A == B || A.up() == B || A.down() == B || A.left() == B || A.right() == B){
                return false;
            }
            return true;
        }
    }

    /* add new direction */
    void new_direction(const size_t curr, const std::vector<Location> next_states, const size_t direction, std::stack<size_t>& stack_list, size_t& org_vert, size_t& edge_vert){
        size_t id = vert_id_.at(curr);
        size_t next; /* next vertex */
        bool next_exists = false; /* if next_vert already present */
        for (auto v : vertices_){
            if (vert_id_.at(v) == 1-id && states_.at(v) == next_states){
                next = v;
                next_exists = true;
                break;
            }
        }
        if (!next_exists){
            next = org_vert;
            org_vert += 1;
            /* insert new vertex for the new_state */ 
            n_vert_ += 1;
            vertices_.insert(next);
            vert_id_[next] = 1-id;
            all_vert_id_[0][next] = vert_id_[next];
            all_vert_id_[1][next] = 1- vert_id_[next];
            states_[next] = next_states;

            /* insert it to stack to explore later */
            stack_list.push(next);
        }
        auto mid = edge_vert;
        edge_vert +=1;
        n_vert_ += 1;
        vertices_.insert(mid);
        vert_id_[mid] = 2;
        all_vert_id_[0][mid] = 2;
        all_vert_id_[1][mid] = 2;
        edges_[curr].insert(mid);
        edges_[mid].insert(next);
        n_edge_ += 2;
        pre_edges_[next].insert(mid);

        labels_[mid] = std::vector<size_t> (ap_id_.size(),2);
        labels_[mid][direction+id*9] = 1;
        if (next_states[id] == Location(1,1)){
            labels_[mid][5+id*9] = 1;
        }
        else if (next_states[id] == Location(boundary_.x_,1)){
            labels_[mid][6+id*9] = 1;
        }
        else if (next_states[id] == Location(1,boundary_.y_)){
            labels_[mid][7+id*9] = 1;
        }
        else if (next_states[id] == Location(boundary_.x_,boundary_.y_)){
            labels_[mid][8+id*9] = 1;
        }
    }


    /* Is there a wall between two location */
    bool isWall(const Location l1, const Location l2, const std::vector<std::vector<Location>> walls) const{
        for (auto wall : walls){
            if (wall[0]==l1 && wall[1]==l2){
                return true;
            }
        }
        return false;
    }



    ///////////////////////////////////////////////////////////////
    /// Generate walls and corridors for a dummy maze
    ///////////////////////////////////////////////////////////////

    /* construct dummy maze with fixed targets, init_states, and generated walls, corridors */
    void dummyMaze(const size_t max_x, const size_t max_y, 
                const std::vector<std::vector<std::vector<Location>>> allWalls,
                const bool standard_safe_distance = true){
        /* initialize boundary, init_state and walls */
        Location boundary(max_x,max_y); 
        std::vector<Location> init_state {mpa::Location(1,1),mpa::Location(max_x,1)}; 
        std::vector<Location> targets{mpa::Location(max_x,max_y),mpa::Location(1,max_y)};

        /* construct the maze */
        mazeSimple(boundary,init_state,targets,allWalls, standard_safe_distance);
    }
    /* same but generate the walls first */
    void dummyMaze(const size_t max_x, const size_t max_y, 
            const size_t n_walls = 0, 
            const size_t  max_corridors = 0){
        /* initialize boundary, init_state and walls */
        Location boundary(max_x,max_y); 
        std::vector<Location> init_state {mpa::Location(1,1),mpa::Location(max_x,1)}; 
        std::vector<Location> targets{mpa::Location(max_x,max_y),mpa::Location(1,max_y)};

        /* generate walls and corridors */
        auto allWalls = genWalls(max_x,max_y,n_walls,max_corridors);

        /* construct the maze */
        mazeSimple(boundary,init_state,targets,allWalls);
    }
    

    /* generate walls and corridors from #walls, #max_corridors */
    std::vector<std::vector<std::vector<Location>>> genWalls(const size_t max_x, const size_t max_y, 
                const size_t n_walls = 0, 
                const size_t  max_corridors = 0){
        
        /* initialize walls and corridors */
        std::vector<std::vector<Location>> walls;
        std::vector<std::vector<Location>> corridors;
                
        /* counters */
        size_t c_rows = max_y-1;
        size_t c_walls = n_walls;
        /* maximum/minimum #walls in a row */
        size_t k_max, k_min;

        while (c_rows > 0){
            k_max = std::min(c_walls,max_x-1);
            k_min = 0;
            if (c_walls>k_max*(c_rows-1)){
                k_min = c_walls-(k_max*(c_rows-1));
                /* sanity check: k_min should not be more than k_max */
                if (k_min>k_max){
                    std::cerr << "Error:dummyMaze:k_min>k_max! Resetting #walls to maximum one.\n";
                    c_walls = k_max*c_rows;
                    k_min = c_walls-k_max*(c_rows-1);
                }
            }
            
            /* generate k = #walls in c_rows-th row */
            size_t k = random_num(k_max,k_min);

            /* generate k radom walls */
            auto walls_indices = random_knum(k,max_x);
            /* add the walls */
            for (size_t w : walls_indices){
                std::vector<mpa::Location> w_wall{mpa::Location(w,c_rows),mpa::Location(w,c_rows+1)};
                walls.push_back(w_wall);
                // std::cout << "w:"<<w <<","<< c_rows << "-->"<<w << ","<<c_rows+1 << "\n";
            }

            /* generate corridors when possible with randomization */
            if (max_corridors > corridors.size() && k < max_x-1){
                /* to randomly choose the corridor: 0:up, 1:none, 2:down */
                size_t range_min = 0;
                size_t range_max = 2;
                for (size_t i = 1; i <= max_x; i++){
                    if (walls_indices.find(i) == walls_indices.end()){
                        size_t decide = random_num(range_max,range_min);
                        if (decide != 1 && decide == 0){ /* add up corridor */
                            std::vector<mpa::Location> i_wall{mpa::Location(i,c_rows),mpa::Location(i,c_rows+1)};
                            corridors.push_back(i_wall);
                            // range_min = 1;
                            // std::cout << "cu:"<<i <<","<< c_rows << "-->"<<i << ","<<c_rows+1 << "\n";
                        }
                        else if (decide != 1 && decide == 2){ /* add down corridor */
                            std::vector<mpa::Location> i_wall{mpa::Location(i,c_rows+1),mpa::Location(i,c_rows)};
                            corridors.push_back(i_wall);
                            // range_max = 1;
                            // std::cout << "cd:"<<i <<","<< c_rows+1 << "-->"<<i << ","<<c_rows << "\n";
                        }
                    }
                }
            }

            /* update counters */
            if (c_walls > k){
                c_walls -= k;
            }
            else{
                c_walls = 0;
            }
            c_rows -= 1;
        }

        /* merge and return all explicit walls and corridors */
        std::vector<std::vector<Location>> explicit_walls = corridors;
        for (auto wall :  walls){
            /* add new both directed wall for each wall */
            auto new_wall = std::vector<Location> {wall[1],wall[0]};
            explicit_walls.push_back(new_wall);
            auto new_wall2 = std::vector<Location> {wall[0],wall[1]};
            explicit_walls.push_back(new_wall2);
        }

        auto allWalls = std::vector<std::vector<std::vector<Location>>> {explicit_walls,explicit_walls};
        n_walls_ = walls.size();
        n_corridors_ = corridors.size();
        return allWalls;
    }

    /* random k number from [1,n] */
    std::set<size_t> random_knum(const size_t k, const size_t max, const size_t min = 1){
        std::set<size_t> result;
        std::vector<size_t> v(max);
        /* fill the vector with the values min, min+1, min+2, ... */
        for (size_t i = 0; i < max; i++){
            v[i] = i + min;
        }
        /* generate k random numbers */
        for (size_t i = 0; i < k; i++){
            size_t rand = random_num(v.size()-1);
            result.insert(v[rand]);
            /* remove rand-th element from v */
            std::swap(v[rand],v[v.size()-1]);
            v.pop_back();
        }
        return result;
    }


    ///////////////////////////////////////////////////////////////
    /// Product of a maze and a dist-game (hoa)
    ///////////////////////////////////////////////////////////////
    
    /* compute product of two games */
    int product_maze_game(const DistGame game1) {
        auto game2 = *this;
        *this = Maze(); /* clear the maze */
        boundary_ = game2.boundary_;
        walls_ = game2.walls_;
        init_state_ = game2.init_state_;
        n_walls_ = game2.n_walls_;
        n_corridors_ = game2.n_corridors_;

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
        for (size_t a : game1.controllable_ap_){
            controllable_ap_.insert(first_aps_map.at(a));
        }
        for (size_t a : game2.controllable_ap_){
            controllable_ap_.insert(second_aps_map.at(a));
        }

        /* initial vertex should be of same player in both games */
        if (game1.all_vert_id_[0].at(game1.init_vert_) != game2.all_vert_id_[0].at(game2.init_vert_)){
            std::cerr << "Error: vertex ids of both games are not same!\n";
        }
        /* insert initial vertex and update all variables */
        n_vert_ = 1;
        size_t org_vert = 1; /* counter for normal vertices */
        size_t edge_vert = (game1.n_vert_-game1.n_edge_/2)*(game2.n_vert_-game2.n_edge_/2); /* counter for edge-vertices */
        vertices_.insert(0); 
        init_vert_ = 0;
        all_vert_id_[0][0] = game1.all_vert_id_[0].at(game1.init_vert_);
        all_vert_id_[1][0] = 1-all_vert_id_[0][0];
        all_colors_[0][0] = game1.all_colors_[0].at(game1.init_vert_);
        all_colors_[1][0] = game1.all_colors_[1].at(game1.init_vert_);
        product_verts.insert({std::make_pair(game1.init_vert_,game2.init_vert_),0});
        state_names_[0] = game2.state_names_.at(game2.init_vert_);

        /* maintain a stack to explore new (product) vertices */
        std::stack<std::vector<size_t>> stack_list;
        /* initialize the stack with initial vertex */
        stack_list.push(std::vector<size_t>{0,game1.init_vert_,game2.init_vert_});

        /* explore until the stack list is empty */
        while (!stack_list.empty()){
            std::vector<size_t> curr = stack_list.top(); /* current pair of vertices */
            stack_list.pop(); /* pop the top element from stack list */
            for (auto u : game1.edges_.at(curr[1])){ /* for each edge-neighbour of 1st vertex */
                for (auto v : game2.edges_.at(curr[2])){ /* for each edge-neihbour of 2nd vertex */
                    bool valid = true; /* if product of these two edges is possible */
                    std::vector<size_t> temp_common(ap_id_.size(),2); /* needed if possible, temporarily store the ids of this new edge */
                    for (auto ap : common_aps){ /* first go through common aps */
                        if (game1.labels_.at(u)[ap[1]] == 0 && game2.labels_.at(v)[ap[2]] == 1){
                            valid = false; /* one of the id of this ap is 1 and other is 0, so not valid product */
                            break;
                        }
                        else if (game1.labels_.at(u)[ap[1]] == 1 && game2.labels_.at(v)[ap[2]] == 2){
                            valid = false; /* one of the id of this ap is 1 and other is 0, so not valid product */
                            break;
                        }
                        /* if this ap is consistent in both edges then add its label to temp_common */
                        else if (game1.labels_.at(u)[ap[1]] == 1 || game2.labels_.at(v)[ap[2]] == 1){
                            temp_common[ap[0]] = 1;
                        }
                        else{
                            temp_common[ap[0]] = 2;
                        }
                    }
                    if (valid){ /* if the product-edge is valid */
                        /* update its temp_common to all ids */
                        for (auto ap : first_aps){
                            temp_common[ap[0]] = game1.labels_.at(u)[ap[1]];
                        }
                        for (auto ap : second_aps){
                            temp_common[ap[0]] = game2.labels_.at(v)[ap[1]];
                        }

                        auto newId = edge_vert;/* new state for product-edges : use edge_vert counter */
                        edge_vert += 1;
                        
                        /* update all variables for new game */
                        n_vert_ += 1;
                        n_edge_ += 2;
                        vertices_.insert(newId);
                        all_vert_id_[0][newId] = 2;
                        all_vert_id_[1][newId] = 2;
                        edges_[curr[0]].insert(newId);
                        labels_[newId] = temp_common;
                        all_colors_[0][newId] = game1.all_colors_[0].at(u);
                        all_colors_[1][newId] = game1.all_colors_[1].at(u);

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
                            all_vert_id_[0][succId] = game1.all_vert_id_[0].at(new_succ.first);
                            all_vert_id_[1][succId] = 1-all_vert_id_[0][succId];
                            edges_[newId].insert(succId);
                            pre_edges_[succId].insert(newId);
                            all_colors_[0][succId] = game1.all_colors_[0].at(new_succ.first);
                            all_colors_[1][succId] = game1.all_colors_[1].at(new_succ.first);
                            product_verts[new_succ] = succId;
                            stack_list.push(std::vector<size_t>{succId,new_succ.first,new_succ.second});
                            state_names_[succId] = game2.state_names_.at(new_succ.second);
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
    /// Buchi objective maze
    ///////////////////////////////////////////////////////////////
    /* construct a color map for a buchi objective */
    std::map<size_t, size_t> buchi_for_maze(const size_t p, const Location loc) const{
        std::map<size_t, size_t> colors;
        for (auto v : vertices_){
            if (vert_id_.at(v)!=2 && states_.at(v)[p] == loc){
                colors[v] = 2;
            }
            else{
                colors[v] = 1;
            }
        }
        return colors;
    }

    /* generate random buchi objective */
    std::map<size_t, size_t> generate_buchi_maze() {
        size_t p = random_num(1);
        return generate_buchi_maze(p);
    }
    std::map<size_t, size_t> generate_buchi_maze(const size_t p) {
        size_t l_x = random_num(boundary_.x_,1);
        size_t l_y = random_num(boundary_.y_,1);
        return buchi_for_maze(p, Location(l_x,l_y));
    }


//     ///////////////////////////////////////////////////////////////
//     /// Basic **Agnes-type** constructor from explicit walls and boundary
//     ///////////////////////////////////////////////////////////////

//     /* construct dummy agnes-format maze with fixed targets, init_states, and generated walls, corridors */
//     negotiation::Negotiate dummyAgnesMaze(const size_t max_x, const size_t max_y, 
//                 const std::vector<std::vector<std::vector<Location>>> allWalls,
//                 const bool standard_safe_distance = true){
//         /* initialize boundary, init_state and walls */
//         Location boundary(max_x,max_y); 
//         std::vector<Location> init_state {mpa::Location(1,1),mpa::Location(max_x,1)}; 
//         std::vector<Location> targets{mpa::Location(max_x,max_y),mpa::Location(1,max_y)};

//         /* construct the maze */
//         return mazeAgnesSimple(boundary,init_state,targets,allWalls, standard_safe_distance);
//     }

//     /* construct from boundary and explicit (directed) walls */
//     negotiation::Negotiate mazeAgnesSimple(const Location boundary, 
//         const std::vector<Location> init_state, 
//         const std::vector<Location> targets, 
//         const std::vector<std::vector<std::vector<Location>>> walls,
//         const bool standard_distance = false){
    
//     /* initialize to use later */
//     boundary_ = boundary;
//     walls_ = walls;
//     init_state_ = init_state;

//     /* initialize component for 0-th robot */
//     negotiation::Component* C0 = new negotiation::Component;
//     std::map<size_t,negotiation::abs_type> temp_posts; /* temporarily store all posts of C0 */
    
//     /* initialize fixed parts for a maze */
//     C0->no_control_inputs = 5;
//     C0->no_dist_inputs = boundary_.x_*boundary_.y_;
//     C0->no_outputs = C0->no_dist_inputs;
//     C0->output_to_state = std::vector<negotiation::abs_type>(C0->no_outputs);

//     /* add initial vertex with init_state_ */
//     C0->no_states = 1;
//     C0->init_.insert(0);
//     states_[0] = init_state_;
//     C0->state_to_output.push_back(boundary_.to_index(init_state_[0]));
//     C0->output_to_state[boundary_.to_index(init_state_[0])] = 0;
    

//     /* stack to explore new vertices */
//     std::stack<size_t> stack_list;
//     stack_list.push(0);

//     /* explore until stack_list is empty */
//     while (!stack_list.empty()){
//         size_t curr = stack_list.top(); /* current vertex to explore */
//         stack_list.pop(); /* pop the top element */

//         /* direction for 0-th robot */
//         for (size_t d0 = 0; d0 < 5; d0++){
//             /* direction for 1-st robot */
//             for (size_t d1 = 0; d1 < 5; d1++){
//                 new_direction_agnes(C0,temp_posts, curr, std::vector<size_t>{d0,d1},stack_list, standard_distance);
//             }
//         }
//     }

//     /* update post using temp_posts */
//     negotiation::abs_type no_post_elems = C0->no_states*C0->no_control_inputs*C0->no_dist_inputs;
//     C0->post = new std::vector<negotiation::abs_type>*[no_post_elems];
//     for (size_t i=0; i<no_post_elems; i++) {
//         std::vector<negotiation::abs_type> *v = new std::vector<negotiation::abs_type>;
//         C0->post[i]=v;
//     }
//     for (auto pair : temp_posts){
//         C0->post[pair.first]->push_back(pair.second);
//     }

//     /* initialize negotiate N and push two identical components C0 */
//     negotiation::Negotiate N;
//     N.components_.push_back(C0);
//     N.components_.push_back(C0);

//     /* initialize safe states to have all vertices */
//     std::unordered_set<negotiation::abs_type>* safes = new std::unordered_set<negotiation::abs_type>;
//     for (size_t v =0; v < C0->no_states; v++){
//         safes->insert(v);
//     }
//     N.safe_states_.push_back(safes);
//     N.safe_states_.push_back(safes);

//     /* initialize target states to have all vertices with states_ in target */
//     std::unordered_set<negotiation::abs_type>* target = new std::unordered_set<negotiation::abs_type>;
//     N.target_states_.push_back(target);
//     N.target_states_.push_back(target);
//     for (size_t v =0; v < C0->no_states; v++){
//         for (size_t i =0; i<2; i++){
//             if (states_[v][i] == targets[i]){
//                 N.target_states_[i]->insert(v);
//             }
//         }
//     }

//     /* initialize the sets of guarantees as all accepting safety automata */
//     for (int c=0; c<2; c++) {
//         negotiation::SafetyAutomaton* s=new negotiation::SafetyAutomaton(N.components_[c]->no_outputs);
//         N.guarantee_.push_back(s);
//     }

//     return N;    
//     }

//     /* add new direction for agnes-type component */
//     int new_direction_agnes(negotiation::Component*& C0, std::map<size_t,negotiation::abs_type>& temp_posts, const size_t curr, const std::vector<size_t> d, std::stack<size_t>& stack_list, const bool standard_distance = false){
//         std::vector<Location> next_states(2); /* new location after the move */
//         for (size_t i = 0; i < 2; i++){ /* for i-th robot */
//             auto state = states_[curr][i];
//             if(d[i]==0){ /* stay */
//                 next_states[i] = states_[curr][0];
//             }
//             else if (d[i]==1 && state.up().inBoundary(boundary_) && !isWall(state,state.up(),walls_[i])){ /* go up */
//                 next_states[i] = state.up();
//             }
//             else if (d[i]==2 && state.y_ > 1 && state.down().inBoundary(boundary_) && !isWall(state,state.down(),walls_[i])){ /* go down */
//                 next_states[i] = state.down();
//             }
//             else if (d[i]==3 && state.right().inBoundary(boundary_) && !isWall(state,state.right(),walls_[i])){ /* go right */
//                 next_states[i] = state.right();
//             }
//             else if (d[i]==4 && state.x_ > 1 && state.left().inBoundary(boundary_) && !isWall(state,state.left(),walls_[i])){ /* go left */
//                 next_states[i] = state.left();
//             }
//             else{
//                 return 0;
//             }
//         }
//         if (!safe_distance(next_states[0],states_[curr][1],standard_distance) || !safe_distance(next_states[1],states_[curr][0],standard_distance)){ /* not in safe distance, so this edge is not possible */
//             return 0;
//         }
        
//         size_t next; /* next vertex */
//         bool next_exists = false; /* if next_vert already present */
//         for (size_t v = 0; v < C0->no_states; v++){
//             if (states_.at(v) == next_states){
//                 next = v;
//                 next_exists = true;
//                 break;
//             }
//         }
//         if (!next_exists){
//             next = C0->no_states;
//             /* insert new vertex for the new_state */ 
//             C0->no_states += 1;
//             states_[next] = next_states;
//             C0->state_to_output.push_back(boundary_.to_index(next_states[0]));
//             C0->output_to_state[boundary_.to_index(next_states[0])] = next;

//             /* insert it to stack to explore later */
//             stack_list.push(next);
//         }
//         temp_posts[C0->addr(curr,d[0],boundary_.to_index(next_states[1]))]=next;
//         return 1;
//     }




}; /* close Maze class definition */
}/* close namespace */

/*! output a maze distgame to gpg (gpg) format 
 * \param[in] DistGame  */
int maze2gpg(const mpa::DistGame G, std::ostream& ostr = std::cout){
    size_t n_vertices = G.n_vert_ - G.n_edge_/2;
    /* print first line */
    ostr<< "parity "<< G.n_vert_-1 <<";\n"; 
        
    for (size_t v = 0; v < n_vertices; v++){ /* print the following for each vertex */
        ostr << v << " "; /* vertex name (number) */
        ostr<<G.all_colors_[0].at(v); /* print color in 1st game separately to avoid comma */
        for (size_t i = 1; i < G.n_games_; i++){/* for each other game print color of v with comma */
            ostr<<","<<G.all_colors_[i].at(v);
        }

        ostr << " " << G.all_vert_id_[0].at(v) << " "; /* print vertex id (which player it belongs to) */
        
        if (!G.edges_.at(v).empty()){ /* if v has neighbours then print them */
            size_t counter = 0;
            for (auto u : G.edges_.at(v)){ /* print all neighbours */
                if (counter == 0){
                    ostr <<*G.edges_.at(u).begin();
                    counter = 1;
                }
                else{
                    ostr <<","<<*G.edges_.at(u).begin();
                }
            }
        }
        ostr<<"\n";
    }
    return 0;
}
#endif