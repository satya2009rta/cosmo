/* FileHandler.hpp */

#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <locale>
#include <vector>
#include <set>
#include <dirent.h>
#include <sys/stat.h>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <functional>

// #include "Maze.hpp"
#include "DistGame.hpp"
#include "hoa_consumer_build_parity_game.hh"


///////////////////////////////////////////////////////////////
/// File/Folder operations
///////////////////////////////////////////////////////////////

/* Check if input is a .hoa file, then we use the hoa game constructor */
inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/*! Ensures that a specified subdirectory exists.
*  \param[in]  dirName     Name of the desired subdirectory.
*/
template <class dir_type>
void checkMakeDir(dir_type dirName) {
    DIR* dir = opendir(dirName);
    if (dir) {
        closedir(dir);
    }
    else if (ENOENT == errno) {
        int result = mkdir(dirName, 0777);
        (void) result;
    }
 }

 /*! Create a file OR erase previous data written to a file
 * \param[in] filename  The name of the file */
void create(const std::string& filename) {
    std::ofstream file;
    file.open(filename, std::ofstream::out | std::ofstream::trunc);
    if (file.is_open()) {
        file.close();
    }
}


///////////////////////////////////////////////////////////////
/// Read games from different formats
///////////////////////////////////////////////////////////////

/*! read a game in extended hoa format from a file/input and convert it to normal game
 * \param[in] filename  Name of the file OR
 * \param[in] issr input stream */
mpa::Game hoa2game(std::istream& issr){
    mpa::Game G;
    cpphoafparser::HOAConsumer::ptr consumer;
    cpphoafparser::parity_game_data data;
    consumer.reset(new cpphoafparser::HOAConsumerBuildParityGame(&data));
    cpphoafparser::HOAParser::parse(issr, consumer);

    G.n_vert_ = data.n_vert;
    G.n_edge_ = data.n_edge;
    G.vertices_ = data.vertices;
    G.edges_ = data.edges;
    G.vert_id_ = data.vert_id;
    G.max_color_ = data.max_color;
    G.colors_ = data.colors;
    G.pre_edges_ = data.pre_edges;
    G.init_vert_ = data.init_vert;
    G.ap_id_ = data.ap_id;
    G.labels_ = data.labels;
    G.controllable_ap_ = data.controllable_ap;

    return G;
}
mpa::Game hoa2game(const std::string filename){
    std::ifstream file(filename);
    return hoa2game(file);
    
}


/*! read a game in pgsolver format from a file/input and convert it to normal game
 * \param[in] filename  Name of the file OR
 * \param[in] issr input stream */
mpa::Game pg2game(std::istream& issr){
    mpa::Game G;
    std::string line;
    /* go through all the lines until a match with arr_name is found */
    while(std::getline(issr,line)) {
        std::stringstream line_stream(line);
        std::string name;
        line_stream >> name;
        /* check for the header line starting with "parity" */
        if (name == "parity"){
            break; /* break the loop when header line found */
        }
    }
    /* initialize everything to be 0 */
    G.n_vert_ = 0;
    G.n_edge_ = 0;
    G.init_vert_ = 0;
    G.max_color_ = 0;
    while(std::getline(issr,line)) {
        std::stringstream line_stream(line);
        G.n_vert_ += 1;
        /* first number is (index of) the vertex */
        size_t vertex;
        line_stream >> vertex;
        G.vertices_.insert(vertex);

        /* second number is colors of that vertex */
        size_t color; 
        line_stream >> color;
        G.colors_[vertex] = color;
        /* update the max_color */
        if (G.max_color_ < color){
            G.max_color_ = color;
        }
        
        /* third number is owner of that vertex */
        size_t vert_id;
        line_stream >> vert_id;
        /* set the owner of that vertex */
        G.vert_id_[vertex] = vert_id;
        
        std::string successors; /* string of successors of that vertex */
        line_stream >> successors;
        std::stringstream ss(successors);
        /* update transitions of that vertex */
        for (size_t i; ss >> i; ++ i) {
            G.edges_[vertex].insert(i);
            G.n_edge_ += 1;
            if (ss.peek() == ',')
                ss.ignore();
        }
    }
    return G;
}
mpa::Game pg2game(const std::string& filename){
    std::ifstream file(filename);
    return pg2game(file);
}

/*! read a game in pgsolver/ehoa format from a file and convert it to normal game
 * \param[in] filename  Name of the file */
mpa::Game file2game(const std::string& filename){
    if(ends_with(filename, ".gm") || ends_with(filename, ".pg")){/* if the file is a .pg (pgsolver) file  */
        return pg2game(filename);
    }
    else if (ends_with(filename, ".hoa") || ends_with(filename, ".ehoa")){ /* if it is a HOA file */
        return hoa2game(filename);
    }
    else{
        std::cerr << "Error: extension of the file is not clear!\n";
        return hoa2game(filename);
    }
}

/* read a game from std::cin */
mpa::Game std2game(){
    size_t pg = 2; /* determine if the format is psolver (default: 2 = undecided) */
    /* construct the game from stdin */
    std::string str, line;
    while (std::getline(std::cin, line)){
        if (pg == 2){
            std::stringstream line_stream(line);
            std::string name;
            line_stream >> name;
            /* check for the header line starting with "parity" then it is pgsolver format */
            if (name == "parity"){
                pg = 1;
                str += "\n" + line;
            }
            /* check for the header line starting with "HOA:" then it is hoa format */
            else if (name == "HOA:"){
                pg = 0;
                str += "\n" + line;
            }
        }
        else{
            str += "\n" + line;
            if (pg == 0 && line.find("END") != std::string::npos){
                break;
            }
            else if (pg == 1 && line == ""){
                break;
            }
        }
    }
    std::istringstream issr(str);
    if (pg == 1){
        return pg2game(issr);
    }
    return hoa2game(issr);
}


///////////////////////////////////////////////////////////////
/// Read dist-games from different formats
///////////////////////////////////////////////////////////////


/*! read a distributed parity game in gpgsolver format from a file/input and convert it to normal dist-game
 * \param[in] filename  Name of the file OR
 * \param[in] issr input stream */
mpa::DistGame gpg2distgame(std::istream& issr){
    mpa::DistGame G;
    G.all_colors_.clear();
    std::string line;
    /* go through all the lines until a match with arr_name is found */
    while(std::getline(issr,line)) {
        std::stringstream line_stream(line);
        std::string name;
        line_stream >> name;
        /* check for the header line starting with "parity" */
        if (name == "parity"){
            break; /* break the loop when header line found */
        }
    }
    /* initialize everything */
    G.n_vert_ = 0;
    G.n_edge_ = 0;
    G.init_vert_ = 0;
    while(std::getline(issr,line)) {
        std::stringstream line_stream(line);
        G.n_vert_ += 1;
        /* first number is (index of) the vertex */
        size_t vertex;
        line_stream >> vertex;
        G.vertices_.insert(vertex);
        /* second string is colors of that vertex */
        std::string colors; /* string of colors of that vertex */
        line_stream >> colors;
        std::stringstream ss_col(colors);
        /* add the colors of that vertex */
        size_t counter_games = 0;
        for (size_t col; ss_col >> col; ++ col){
            if (counter_games == G.all_colors_.size()){
                G.all_colors_.push_back(std::map<size_t, size_t>());
                G.all_max_color_.push_back(0);
            }
            G.all_colors_[counter_games][vertex] = col;
            /* update the max_colors */
            if (G.all_max_color_[counter_games] < col){
                G.all_max_color_[counter_games] = col;
            }
            if (ss_col.peek() == ',')
                ss_col.ignore();
            counter_games += 1;
        }
        
        /* third number is owner of that vertex */
        size_t vert_id;
        line_stream >> vert_id;
        /* set the owner of that vertex */
        G.all_vert_id_[0][vertex] = std::min(vert_id, 2-vert_id);
        G.all_vert_id_[1][vertex] = 1-G.all_vert_id_[0][vertex];
        
        std::string successors; /* string of successors of that vertex */
        line_stream >> successors;
        std::stringstream ss(successors);
        /* update transitions of that vertex */
        for (size_t i; ss >> i; ++ i) {
            G.edges_[vertex].insert(i);
            G.n_edge_ += 1;
            if (ss.peek() == ',')
                ss.ignore();
        }
    }
    G.n_games_ = G.all_colors_.size();
    return G;
}
mpa::DistGame gpg2distgame(const std::string& filename){
    std::ifstream file(filename);
    return gpg2distgame(file);
}

/*! read a game in extended hoa format from a file/input and convert it to dist game
 * \param[in] filename  Name of the file OR
 * \param[in] issr input stream */
mpa::DistGame hoa2distgame(std::istream& issr){
    mpa::DistGame G;
    cpphoafparser::HOAConsumer::ptr consumer;
    cpphoafparser::parity_game_data data;
    consumer.reset(new cpphoafparser::HOAConsumerBuildParityGame(&data));
    cpphoafparser::HOAParser::parse(issr, consumer);

    G.n_vert_ = data.n_vert;
    G.n_edge_ = data.n_edge;
    G.vertices_ = data.vertices;
    G.edges_ = data.edges;
    G.pre_edges_ = data.pre_edges;
    G.init_vert_ = data.init_vert;
    G.ap_id_ = data.ap_id;
    G.labels_ = data.labels;
    G.controllable_ap_ = data.controllable_ap;
    G.n_games_ = 2;

    if (data.all_colors.size() > 1){
        G.n_games_ = data.all_colors.size();
        for (auto v : G.vertices_){
            for (size_t i = 0; i < G.n_games_; i++){
                if (data.all_colors[i].find(v) == data.all_colors[i].end()){
                    data.all_colors[i][v] = data.minCol;
                }
            }
        }
        G.all_colors_ = data.all_colors;
        for (size_t i = 0; i < G.n_games_; i++){
            G.all_max_color_.push_back(G.max_col(G.all_colors_[i]));
        }       
        G.all_vert_id_[0] = data.vert_id;
        for (auto pair : data.vert_id){
            G.all_vert_id_[1][pair.first] = 1 - pair.second;
        }
        return G;
    }
    G.vert_id_ = data.vert_id;
    G.max_color_ = data.max_color;
    G.colors_ = data.colors;
    return G;
}
mpa::Game hoa2distgame(const std::string filename){
    std::ifstream file(filename);
    return hoa2distgame(file);
}
mpa::Game hoa2distgame(const std::string filename, const std::string filename2){
    mpa::DistGame G;
    G.product_games(hoa2game(filename),hoa2game(filename2));
    return std::move(G);
}

/* read a dist-game from std::cin */
mpa::DistGame std2distgame(){
    size_t pg = 2; /* determine if the format is psolver (default: 2 = undecided) */
    /* construct the game from stdin */
    std::string str, line;
    while (std::getline(std::cin, line)){
        if (pg == 2){
            std::stringstream line_stream(line);
            std::string name;
            line_stream >> name;
            /* check for the header line starting with "parity" then it is pgsolver format */
            if (name == "parity"){
                pg = 1;
                str += "\n" + line;
            }
            /* check for the header line starting with "HOA:" then it is hoa format */
            else if (name == "HOA:"){
                pg = 0;
                str += "\n" + line;
            }
        }
        else{
            str += "\n" + line;
            if (pg == 0 && line.find("END") != std::string::npos){
                break;
            }
            else if (pg == 1 && line == ""){
                break;
            }
        }
    }
    std::istringstream issr(str);
    if (pg == 1){
        return gpg2distgame(issr);
    }

    mpa::DistGame G = hoa2distgame(issr);
    if (G.vert_id_.empty()){ /* if it is a dist-game */
        return G;
    }

    /* other wise construct the second game from stdin and take product */
    std::string str2, line2;
    bool start = false;
    while (std::getline(std::cin, line2)){
        if (!start){
            if (line2.find("HOA:") != std::string::npos){
                start = true;
                str2 += "\n" + line2;
            }
        }
        else{
            str2 += "\n" + line2;
            if (line2.find("END") != std::string::npos){
                break;
            }
        }
    }
    if (!start){ /* if it only one game */
        G.n_games_ = 1;
        G.all_colors_[0] = G.colors_;
        G.all_max_color_[0] = G.max_col(G.colors_);
        G.all_vert_id_[0] = G.vert_id_;
        return G;
    }
    std::istringstream issr2(str2);
    mpa::Game G2 = hoa2game(issr2);
    mpa::Game G1 = G;
    G.product_games(G1,G2);
    return G;
}

/*! read a game in pgsolver/ehoa format from a file and convert it to dist game
 * \param[in] filename  Name of the file */
mpa::DistGame file2distgame(const std::string& filename){
    if(ends_with(filename, ".gm") || ends_with(filename, ".pg")){/* if the file is a .pg (pgsolver) file  */
        return gpg2distgame(filename);
    }
    else if (ends_with(filename, ".hoa") || ends_with(filename, ".ehoa")){ /* if it is a HOA file */
        return hoa2distgame(filename);
    }
    else{
        std::cerr << "Error: extension of the file is not clear!\n";
        return hoa2distgame(filename);
    }
}
mpa::DistGame file2distgame(const std::string& filename, const std::string& filename2){
    return hoa2distgame(filename, filename2);
}


///////////////////////////////////////////////////////////////
/// Write games to different formats
///////////////////////////////////////////////////////////////

/* function: print_acceptance */
std::string print_acceptance(const size_t k) {
    std::string output;
    if (k == 0){
        output = "Inf(0)";
    }
    else if (k == 1){
        output = "Fin(1) & Inf(0)";
    }
    else if (k%2 == 1){
        output = "Fin("+std::to_string(k)+") & (" + print_acceptance(k-1) + ")";
    }
    else{
        output = "Inf("+std::to_string(k)+") | (" + print_acceptance(k-1) + ")";
    }
    return output;
}

/* function: print_label
*
* print out the label of state */
std::string print_label(const mpa::Game G, size_t u, const bool num = true) {
    std::string output;
    std::string output_num;
    size_t start = 1;
    for (size_t i = 0; i < G.ap_id_.size(); i++){
        if (G.labels_.at(u)[i] == 0){
            if (start == 1){
                start = 0;
                output += "!"+G.ap_id_.at(i);
                output_num += "!" + std::to_string(i);
            }
            else{
                output += " & !"+G.ap_id_.at(i);
                output_num += "&!" + std::to_string(i);
            }
        }
        else if (G.labels_.at(u)[i] == 1){
            if (start == 1){
                start = 0;
                output += G.ap_id_.at(i);
                output_num += std::to_string(i);
            }
            else{
                output += " & "+ G.ap_id_.at(i);
                output_num += "&"+std::to_string(i);
            }
        }
    }
    if (output.empty()){
        output = "*";
        output_num = "t";
    }
    if (!num){
        return output;
    }
    return output_num;
}

/*! print out (to a file if given) the game in hoa format 
 * \param[in] Game  */  
int game2hoa(const mpa::Game G, std::ostream& ostr = std::cout) {
    size_t n_vertices = G.n_vert_ - G.n_edge_/2;
    ostr << "HOA: v1\n";
    ostr << "States: "<< n_vertices << "\n";
    ostr << "Start: "<< G.init_vert_ << "\n";
    ostr << "AP: "<< G.ap_id_.size();
    for (auto pair : G.ap_id_){
        ostr << " \"" << pair.second << "\"";
    }
    ostr << "\nacc-name: parity max even "<< G.max_color_+1<< "\n";
    ostr << "Acceptance: " << G.max_color_+1 << " " << print_acceptance(G.max_color_)<<"\n";
    ostr << "properties: trans-labels explicit-labels trans-acc deterministic\n";
    ostr << "spot-state-player:";
    for (size_t i = 0; i < n_vertices; i++){
        ostr << " " << 1-G.vert_id_.at(i);
    }
    ostr << "\ncontrollable-AP:";
    for (auto a : G.controllable_ap_){
        ostr << " " << a;
    }
    /* print the body */
    ostr << "\n--BODY--\n";
    for (size_t u = 0; u < n_vertices; u++){
        if (G.vert_id_.at(u) != 2){
            ostr << "State: "<< u;
            if (!G.state_names_.empty()){
                ostr << " \""<<G.state_names_.at(u)<<"\"";
            }
            ostr << "\n";
            for (auto v : G.edges_.at(u)){
                ostr << "[" << print_label(G,v,true) << "] ";
                ostr << *G.edges_.at(v).begin();
                ostr << " {" << G.colors_.at(v)<<"}";
                ostr << "\n";
            }
        }
    }
    ostr << "--END--\n";
    return 1;
}


///////////////////////////////////////////////////////////////
/// Write dist-games to different formats
///////////////////////////////////////////////////////////////


/*! print out (to a file if given) the dist-game in hoa format 
 * \param[in] DistGame  */
int distgame2hoa(mpa::DistGame& G, std::ostream& ostr = std::cout) {
    size_t n_vertices = G.n_vert_ - G.n_edge_/2;
    ostr << "HOA: v1\n";
    ostr << "States: "<< n_vertices << "\n";
    ostr << "Start: "<< G.init_vert_ << "\n";
    ostr << "AP: "<< G.ap_id_.size();
    for (auto pair : G.ap_id_){
        ostr << " \"" << pair.second << "\"";
    }
    size_t max_color = 0;
    for (auto c : G.all_max_color_){
        max_color = std::max(max_color, c);
    }
    
    ostr << "\nacc-name: parity max even "<< max_color +1<< "\n";
    ostr << "Acceptance: " << max_color +1 << " " << print_acceptance(max_color)<<"\n";
    ostr << "properties: trans-labels explicit-labels trans-acc deterministic\n";
    ostr << "spot-state-player:";
    for (size_t i = 0; i < n_vertices; i++){
        ostr << " " << 1-G.all_vert_id_[0].at(i);
    }
    ostr << "\ncontrollable-AP:";
    for (auto a : G.controllable_ap_){
        ostr << " " << a;
    }
    /* print the body */
    ostr << "\n--BODY--\n";
    for (size_t u = 0; u < n_vertices; u++){
        if (G.all_vert_id_[0].at(u) != 2){
            ostr << "State: "<< u;
            if (!G.state_names_.empty()){
                ostr << " \""<<G.state_names_.at(u)<<"\"";
            }
            ostr << "\n";
            for (auto v : G.edges_.at(u)){
                ostr << "[" << print_label(G,v,true) << "] ";
                ostr << *G.edges_.at(v).begin();
                ostr << " {";
                // G.all_colors_[0].at(v) << " "<<G.all_colors_[1].at(v);
                for (size_t i = 0; i < G.n_games_; i++){
                    ostr << G.all_colors_[i].at(v);
                    if (i < G.n_games_-1){
                        ostr << " ";
                    }
                }
                ostr << "}\n";
            }
        }
    }
    ostr << "--END--\n";
    return 1;
}

/*! output a distgame to gpg (gpg) format 
 * \param[in] DistGame  */
int distgame2gpg(mpa::DistGame& G, std::ostream& ostr = std::cout){
    /* print first line */
    ostr<< "parity "<< G.n_vert_-1 <<";\n"; 
        
    for (size_t v = 0; v < G.n_vert_; v++){ /* print the following for each vertex */
        ostr << v << " "; /* vertex name (number) */
        ostr<<G.all_colors_[0].at(v); /* print color in 1st game separately to avoid comma */
        if (G.n_games_ > 1){
            for (size_t i = 1; i < G.n_games_; i++){/* for each other game print color of v with comma */
                ostr<<","<<G.all_colors_[i].at(v);
            }
        }

        ostr << " " << G.all_vert_id_[0].at(v) << " "; /* print vertex id (which player it belongs to) */
        
        if (!G.edges_.at(v).empty()){ /* if v has neighbours then print them */
            size_t counter = 0;
            for (auto u : G.edges_.at(v)){ /* print all neighbours */
                if (counter == 0){
                    ostr <<u;
                    counter = 1;
                }
                else{
                    ostr <<","<<u;
                }
            }
        }
        ostr<<"\n";
    }
    return 0;
}

/*! output a distgame to file/output
 * \param[in] DistGame  */
int distgame2std(mpa::DistGame& G, std::ostream& ostr = std::cout, const std::string format = "hoa"){
    if (G.labels_.empty() || format != "hoa"){
        return distgame2gpg(G, ostr);
    }
    return distgame2hoa(G, ostr);
}