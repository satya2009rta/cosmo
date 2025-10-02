/*
 * Class: Template
 *
 *  Class formalizing a template for assumptions and strategies
 */

#ifndef TEMPLATE_HPP_
#define TEMPLATE_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <algorithm>
#include <stack> 

namespace mpa {
class Template {
public:
    /* set of unsafe edges */
    std::map<size_t, std::set<size_t>> unsafe_edges_;
    /* set of colive edges */
    std::map<size_t, std::set<size_t>> colive_edges_;
    /* vector of condition sets */
    std::vector<std::set<size_t>> cond_sets_;
    /* a vector of live groups */
    std::vector<std::map<size_t, std::set<size_t>>> live_groups_;
    /* set of conditional live groups (mapping from the conditions sets) */
    std::vector<std::vector<std::map<size_t, std::set<size_t>>>> cond_live_groups_;

    /* set of unsafe strings */
    std::map<size_t, std::set<std::string>> unsafe_strings_;
    /* set of colive strings */
    std::map<size_t, std::set<std::string>> colive_strings_;
    /* a vector of live groups for strings */
    std::vector<std::map<size_t, std::set<std::string>>> live_groups_strings_;
    /* set of conditional live groups for strings (mapping from the conditions sets) */
    std::vector<std::vector<std::map<size_t, std::set<std::string>>>> cond_live_groups_strings_;
public:
    /* default constructor */
    Template() {}


    ///////////////////////////////////////////////////////////////
    /// Operators
    ///////////////////////////////////////////////////////////////

    /* merge a new template */
    void merge(const Template& new_temp){
        edge_merge(unsafe_edges_, new_temp.unsafe_edges_);
        edge_merge(colive_edges_, new_temp.colive_edges_);
        live_groups_.insert(live_groups_.end(), new_temp.live_groups_.begin(), new_temp.live_groups_.end());
        cond_sets_.insert(cond_sets_.end(),new_temp.cond_sets_.begin(), new_temp.cond_sets_.end());
        cond_live_groups_.insert(cond_live_groups_.end(),new_temp.cond_live_groups_.begin(),new_temp.cond_live_groups_.end());
    }
    /* merge new set of edges to this set of edges */
    void edge_merge(std::map<size_t, std::set<size_t>>& edges, const std::map<size_t, std::set<size_t>>& new_edges) const {
        for (auto v = new_edges.begin(); v != new_edges.end(); v++){
            edges[v->first].insert(v->second.begin(), v->second.end());
        }
    }

    /* merge all colive edges and live groups (in cond_live_groups_) with another template */
    void merge_live_colive(const Template& new_temp){
        edge_merge(colive_edges_, new_temp.colive_edges_);
        live_groups_.insert(live_groups_.end(), new_temp.live_groups_.begin(), new_temp.live_groups_.end());
        for (const auto& live_groups : cond_live_groups_){
            live_groups_.insert(live_groups_.end(), live_groups.begin(), live_groups.end());
        }
        for (const auto& live_groups : new_temp.cond_live_groups_){
            live_groups_.insert(live_groups_.end(), live_groups.begin(), live_groups.end());
        }
    }

    /* merge all colive edges and live groups (in cond_live_groups_) in a set of templates */
    void merge_live_colive(const std::vector<Template>& new_temps){
        for (const auto& new_temp : new_temps){
            merge_live_colive(new_temp);
        }
    }

    ///////////////////////////////////////////////////////////////
    /// Size of the template
    ///////////////////////////////////////////////////////////////
    size_t size_unsafe() const {
        return unsafe_edges_.size();
    }
    size_t size_colive() const {
        return colive_edges_.size();
    }
    size_t size_live() const {
        return live_groups_.size();
    }
    size_t size_cond_live() const {
        size_t result = 0;
        for (const auto& live_groups : cond_live_groups_){
            result += live_groups.size();
        }
        return result;
    }

    ///////////////////////////////////////////////////////////////
    /// Clean/clear elements
    ///////////////////////////////////////////////////////////////

    /* clear the template */
    void clear(){
        unsafe_edges_.clear();
        colive_edges_.clear();
        live_groups_.clear();
        cond_sets_.clear();
        cond_live_groups_.clear();
    }

    /* clean the keys with empty value in a map of edges */
    void clean_edges(std::map<size_t, std::set<size_t>>& edges) const {
        for (auto v = edges.begin(); v != edges.end();) {
            if (v->second.empty()){
                edges.erase(v++);
            }
            else{
                ++v;
            }
        }
    }

    /* clean keys with empty value in unsafe edges */
    void clean_unsafe_edges(){
        clean_edges(unsafe_edges_);
    }

    /* clean keys with empty value in colive edges */
    void clean_colive_edges(){
        clean_edges(colive_edges_);
    }

    /* clean empty live group in a vector of live groups */
    void clean_live_groups(std::vector<std::map<size_t, std::set<size_t>>>& live_groups) const {
        std::vector<std::map<size_t, std::set<size_t>>> new_live_groups;
        for (auto& live_group : live_groups){
            clean_edges(live_group);
            if (!live_group.empty()){
                new_live_groups.push_back(live_group);
            }
        }
        live_groups = new_live_groups;
        /* remove duplicate live groups */
        std::sort(live_groups.begin(), live_groups.end());
        live_groups.erase(std::unique(live_groups.begin(), live_groups.end()), live_groups.end());
    }
    /* clean empty live group in live_groups_ */
    void clean_live_groups() {
        clean_live_groups(live_groups_);
    }

    /* clean empty live groups from cond_live_groups_ */
    void clean_cond_live_groups(){
        std::vector<std::set<size_t>> new_cond_sets;
        std::vector<std::vector<std::map<size_t, std::set<size_t>>>> new_cond_live_groups;
        for (size_t i = 0; i < cond_sets_.size(); i++){
            if (!cond_sets_[i].empty()){
                clean_live_groups(cond_live_groups_[i]);
                if (!cond_live_groups_[i].empty()){
                    new_cond_sets.push_back(cond_sets_[i]);
                    new_cond_live_groups.push_back(cond_live_groups_[i]);
                }
            }
        }
        cond_sets_ = new_cond_sets;
        cond_live_groups_ = new_cond_live_groups;
        /* merge cond_sets if the corresponding live groups are the same (and remove the duplicate) */
        for (size_t i = 0; i < cond_sets_.size(); i++){
            for (size_t j = i+1; j < cond_sets_.size(); j++){
                if (cond_live_groups_[i] == cond_live_groups_[j]){
                    cond_sets_[i].insert(cond_sets_[j].begin(), cond_sets_[j].end());
                    cond_sets_.erase(cond_sets_.begin()+j);
                    cond_live_groups_.erase(cond_live_groups_.begin()+j);
                    j -= 1;
                }
            }
        }
    }

    /* clean every empty things from the template */
    void clean(){
        clean_unsafe_edges();
        clean_colive_edges();
        clean_live_groups();
        clean_cond_live_groups();
    }

    ///////////////////////////////////////////////////////////////
    /// Compute strings from template if needed
    ///////////////////////////////////////////////////////////////

    /* convert set of size_t to set of strings */
    std::set<std::string> set_to_string(const std::set<size_t>& set) const {
        std::set<std::string> result;
        for (auto element : set){
            result.insert(std::to_string(element));
        }
        return result;
    }

    /* compute strings from the template */
    void compute_strings() {
        for (auto& pair : unsafe_edges_){
            unsafe_strings_[pair.first] = set_to_string(pair.second);
        }
        for (auto& pair : colive_edges_){
            colive_strings_[pair.first] = set_to_string(pair.second);
        }
        for (auto& live_group : live_groups_){
            std::map<size_t, std::set<std::string>> live_group_strings;
            for (auto& pair : live_group){
                live_group_strings[pair.first] = set_to_string(pair.second);
            }
            live_groups_strings_.push_back(live_group_strings);
        }
        for (auto& live_groups : cond_live_groups_){
            std::vector<std::map<size_t, std::set<std::string>>> cond_live_strings;
            for (auto& live_group : live_groups){
                std::map<size_t, std::set<std::string>> live_group_strings;
                for (auto& pair : live_group){
                    live_group_strings[pair.first] = set_to_string(pair.second);
                }
                cond_live_strings.push_back(live_group_strings);
            }
            cond_live_groups_strings_.push_back(cond_live_strings);
        }
    }


    
    ///////////////////////////////////////////////////////////////
    /// Print functions
    ///////////////////////////////////////////////////////////////

    /* function: print_size
     *
     * print out size of all template */
    void print_size(){
        clean();
        std::cout << "\n";
        std::cout << "#unsafe_edges:     "<< size_unsafe()<<"\n";
        std::cout << "#colive_edges:     "<< size_colive()<<"\n";
        std::cout << "#cond_live_groups: "<< size_cond_live() <<"\n";
    }

    
    /* function: print_edges
     *
     * print out all edges of the map */
    void print_edges (const std::map<size_t, std::set<std::string>>& edges_string, const std::string note = "edges", const int print_empty = 1) const {
        if (print_empty == 1 || edges_string.size()!=0){
            std::cout << "\n" << note << ": \n";
            for (auto v = edges_string.begin(); v != edges_string.end(); ++v){
                for (const auto& u : v->second){
                    std::cout << v->first << " -> " << u << "\n";
                }
            }
        }
    }

    /* function: print_unsafe_edges
     *
     * print out all unsafe edges */
    void print_unsafe_edges (const std::string note = "Unsafe edges", 
                            const int print_empty = 0) const {
    print_edges(unsafe_strings_, note, print_empty);
    }

    /* function: print_colive_edges
     *
     * print out all colive edges */
    void print_colive_edges (const std::string note = "Colive edges", 
                            const int print_empty = 0) const {
    print_edges(colive_strings_, note, print_empty);
    }


    /* function: print_live_groups
     *
     * print out all live_groups of the live_groups */
    void print_live_groups (const std::vector<std::map<size_t, std::set<std::string>>>& live_groups_strings, const std::string note = "live groups", const int print_empty = 1) const {
        if (print_empty == 1 || live_groups_strings.size()!=0){
            std::cout << "\n" << note << ": \n";
            for (const auto& live_group : live_groups_strings){
                std::cout << "{";
                size_t counter = 0;
                for (auto v = live_group.begin(); v != live_group.end(); v++){
                    for (const auto& u : v->second){
                        if (counter == 0){
                            std::cout << "("<<v->first << " -> " << u<<")";
                            counter += 1;
                        }
                        else{
                            std::cout << ", " <<"("<< v->first << " -> " << u<<")";
                        }
                    }
                }
                std::cout << "}\n";
            }
        }
    }
    void print_live_groups(const std::string note = "Live groups", const int print_empty = 0) const {
        print_live_groups(live_groups_strings_, note, print_empty);
    }


    /* function: print_cond_live_groups
     *
     * print out all conditional live_groups of the cond_live_groups_ */
    void print_cond_live_groups (const std::string note = "Conditional live groups", 
                                const int print_empty = 0) const {
        if (print_empty == 1 || cond_live_groups_strings_.size()!=0){
            std::cout << "\n" << note << ": \n";
            for (size_t i = 0; i < cond_sets_.size(); i++){
                const auto& set = cond_sets_[i];
                const auto& live_groups_strings = cond_live_groups_strings_[i];
                for (const auto& live_group_string : live_groups_strings){
                    std::cout << "{";
                    for (auto u=set.begin(); u != set.end();){
                        std::cout << *u;
                        ++u;
                        if (u != set.end())
                            std::cout << ", ";
                    }
                    std::cout << "}: ";
                
                    std::cout << "{";
                    size_t counter = 0;
                    for (const auto& v : live_group_string){
                        for (const auto& u : v.second){
                            if (counter == 0){
                                std::cout << "("<<v.first << " -> " << u<<")";
                                counter += 1;
                            }
                            else{
                                std::cout << ", " <<"("<< v.first << " -> " << u<<")";
                            }
                        }
                    }
                    std::cout << "}\n";
                }
            }
        }
    }

    /* function: print_template
     *
     * print out the whole template */
    void print_template(){
        if (unsafe_strings_.empty() && colive_strings_.empty() && live_groups_strings_.empty() && cond_live_groups_strings_.empty()){
            compute_strings();
        }
        print_unsafe_edges();
        print_colive_edges();
        print_live_groups();
        print_cond_live_groups();
        if (unsafe_edges_.empty() && colive_edges_.empty() && live_groups_.empty() && cond_live_groups_.empty()){
            std::cout << "\nTRUE\n";
        }
    }

}; /* close class definition */
} /* close namespace */

#endif
