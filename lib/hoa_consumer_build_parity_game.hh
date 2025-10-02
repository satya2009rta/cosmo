//==============================================================================
//
// Author: Satya Prakash Nayak
//
// An interface to the cpphoafparser library for parsing a ParityGame automaton in
// HOA format
//==============================================================================

#ifndef CPPHOAFPARSER_HOACONSUMER_BUILDPARITYGAME_H
#define CPPHOAFPARSER_HOACONSUMER_BUILDPARITYGAME_H

#include <iostream>
#include <limits>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <cstring>
#include <string>

#include "cpphoafparser/consumer/hoa_consumer.hh"
#include "cpphoafparser/parser/hoa_parser.hh"

#define V0 0 /* vertices belonging to player 0 */
#define V1 1 /* vertices belonging to player 1 */

#define UNUSED(x) (void)(x)

namespace cpphoafparser {


    struct parity_game_data {
        /* number of vertices */
        size_t n_vert;
        /* number of edges */
        size_t n_edge;
        /* vertices */
        std::set<size_t> vertices;
        /* vertex id: V0, when the vertex belongs to player 0 and V1 when it belongs to player 1 */
        std::map<size_t, size_t> vert_id;
        /* edges as a map from vertices to set of its neighbours */
        std::map<size_t, std::set<size_t>> edges;
        /* number of colors */
        size_t max_color;
        /* colors of vertices */
        std::map<size_t, size_t> colors;
        /* vector of colors of vertices for multi-objective games */
        std::vector<std::map<size_t, size_t>> all_colors;

        /* pre of edges for original vertices (incoming edges) */
        std::map<size_t, std::set<size_t>> pre_edges;
        
        /* initial vertex */
        size_t init_vert;
        /* ids of atomic proposition */
        std::map<size_t, std::string> ap_id;
        /* labels of edges (mid-states) */
        std::map<size_t, std::vector<size_t>> labels;
        /* controllable APs */
        std::set<size_t> controllable_ap;
        /* name of states */
        std::map<size_t, std::string> state_names;
        
        std::string acc_name;
        size_t minCol;
    };

    /**
 * A HOAConsumer implementation that works as an interface to convert a parity_game automaton from HOA format
 * to the format used in Mascot-SDS.
 *
 */

    class HOAConsumerBuildParityGame : public HOAConsumer {
    public:
        parity_game_data *data_;

        /** Constructor, providing a reference to the output stream */
        HOAConsumerBuildParityGame(parity_game_data *data) : data_(data), out(std::cout) {}

        virtual bool parserResolvesAliases() override {
            return false;
        }

        virtual void notifyHeaderStart(const std::string &version) override {
            UNUSED(version);
            UNUSED(out);
        }

        virtual void setNumberOfStates(unsigned int numberOfStates) override {
            data_->n_vert = numberOfStates;
            data_->max_color = 1;
            data_->n_edge = 0;
        }

        virtual void addStartStates(const int_list &stateConjunction) override {
            if (stateConjunction.size() > 0) {
                data_->init_vert = stateConjunction[0];
            }
        }

        virtual void addAlias(const std::string &name, label_expr::ptr labelExpr) override {
            UNUSED(name);
            UNUSED(labelExpr);
        }

        virtual void setAPs(const std::vector<std::string> &aps) override {
            data_->ap_id.clear();
            for (size_t i = 0; i < aps.size(); i++) {
                data_->ap_id.insert({i, aps[i]});
            }
        }

        virtual void setAcceptanceCondition(unsigned int numberOfSets, acceptance_expr::ptr accExpr) override {
            // for (size_t i = 0; i < numberOfSets; i++) {
            //     std::vector<size_t> acc_label;
            //     data_->acc_signature.push_back(acc_label);
            // }

            data_->minCol = 0;
            std::stack<acceptance_expr::ptr> nodes;
            nodes.push(accExpr);
            while (nodes.size() != 0) {
                acceptance_expr::ptr curr_node = nodes.top();
                nodes.pop();
                if (curr_node->isOR()) {
                    nodes.push(curr_node->getLeft());
                    nodes.push(curr_node->getRight());
                } 
                else if (curr_node->isAND()) {
                    nodes.push(curr_node->getLeft());
                    nodes.push(curr_node->getRight());
                }
                else if (curr_node->isAtom()){
                    if (curr_node->getAtom().getType() == AtomAcceptance::TEMPORAL_FIN){
                        if (curr_node->getAtom().getAcceptanceSet()%2 == 1){
                            data_->minCol = 1;
                        }
                        else{
                            data_->minCol = 0;
                        }
                    }
                    else if (curr_node->getAtom().getType() == AtomAcceptance::TEMPORAL_INF){
                        if (curr_node->getAtom().getAcceptanceSet()%2 == 0){
                            data_->minCol = 1;
                        }
                        else{
                            data_->minCol = 0;
                        }
                    }
                    break;
                }
            }
            data_->max_color = data_->minCol;
            UNUSED(numberOfSets);
            // UNUSED(accExpr);
        }

        virtual void provideAcceptanceName(const std::string &name, const std::vector<IntOrString> &extraInfo) override {
            data_->acc_name = name;
            UNUSED(extraInfo);
        }

        virtual void setName(const std::string &name) override {
            UNUSED(name);
        }

        virtual void setTool(const std::string &name, std::shared_ptr<std::string> version) override {
            UNUSED(name);
            UNUSED(version);
        }

        virtual void addProperties(const std::vector<std::string> &properties) override {
            UNUSED(properties);
        }

        // virtual void addMiscHeader(const std::string &name, const std::vector<IntOrString> &content) override {
        //     if (strcmp(name.c_str(), "controllable-AP")) {
        //         throw std::runtime_error("The set controllable-AP is not provided");
        //     } else {
        //         data_->controllable_ap = content;
        //     }
        // }

        virtual void addMiscHeader(const std::string &name, const std::vector<IntOrString> &content) override {
            if (strcmp(name.c_str(), "spot-state-player") == 0) {
                for (size_t i = 0; i < data_->n_vert; i++){
                    data_->vert_id.insert({i, 1-content[i].getInteger()});
                }
            }
            
            if (strcmp(name.c_str(), "controllable-AP") == 0) {
                for (auto u : content){
                    data_->controllable_ap.insert(u.getInteger());
                }
            }
        }

        virtual void notifyBodyStart() override {}

        virtual void addState(unsigned int id,
                              std::shared_ptr<std::string> info,
                              label_expr::ptr labelExpr,
                              std::shared_ptr<int_list> accSignature) override {
            data_->vertices.insert(id);
            size_t color;
            std::vector<size_t> all_color;
            if (accSignature) {
                std::vector<unsigned int> accSignatureValue = *accSignature.get();
                color = accSignatureValue[0]+data_->minCol+1;
                if (color > data_->max_color){
                    data_->max_color = color;
                }
                /* continue collecting colors as long as possible */
                for (size_t i = 0; i < accSignatureValue.size(); i++){
                    all_color.push_back(accSignatureValue[i]+data_->minCol+1);
                }
                if (all_color.size() > 0 && all_color.size() != data_->all_colors.size()){
                    data_->all_colors = std::vector<std::map<size_t, size_t>>(all_color.size());
                }
                
                
            }
            else{
                color = data_->minCol;
            }
            data_->colors.insert({id,color});
            for (size_t i = 0; i < all_color.size(); i++){
                data_->all_colors[i].insert({id,all_color[i]});
            }
            
            // if (data_->vert_id[id] == V0){
            //     UNUSED(labelExpr);
            // } else{
            //     std::stack<label_expr::ptr> nodes;
            //     nodes.push(labelExpr);
            //     while (nodes.size() != 0) {
            //         label_expr::ptr curr_node = nodes.top();
            //         nodes.pop();
            //         if (curr_node->isAND()) {
            //             nodes.push(curr_node->getLeft());
            //             nodes.push(curr_node->getRight());
            //         } else if (curr_node->isNOT()){
            //             nodes.push(curr_node->getLeft());
            //         } else if (curr_node->isAtom()){
            //             unsigned int atom_id = curr_node->getAtom().getAPIndex();
            //             if (atom_id == data_->controllable_ap){
            //                 data_->vert_id[id] = V0;
            //                 break;
            //             }
            //         }
            //     }
            // }
            UNUSED(labelExpr);
            // UNUSED(info);
            if (info){
                data_->state_names.insert({id,*info});
            }
            else{
                data_->state_names.insert({id,""});
            }
        }

        virtual void addEdgeImplicit(unsigned int stateId,
                                     const int_list &conjSuccessors,
                                     std::shared_ptr<int_list> accSignature) override {
            UNUSED(stateId);
            UNUSED(conjSuccessors);
            UNUSED(accSignature);
        }

        virtual void addEdgeWithLabel(unsigned int stateId,
                                      label_expr::ptr labelExpr,
                                      const int_list &conjSuccessors,
                                      std::shared_ptr<int_list> accSignature) override {
            size_t nbrId = conjSuccessors[0];
            size_t color;
            std::vector<size_t> all_color;
            if (accSignature) {
                std::vector<unsigned int> accSignatureValue = *accSignature.get();
                color = accSignatureValue[0]+data_->minCol+1;
                if (color > data_->max_color){
                    data_->max_color = color;
                }
                /* continue collecting colors as long as possible */
                for (size_t i = 0; i < accSignatureValue.size(); i++){
                    all_color.push_back(accSignatureValue[i]+data_->minCol+1);
                }
                if (all_color.size() > 0 && all_color.size() != data_->all_colors.size()){
                    data_->all_colors = std::vector<std::map<size_t, size_t>>(all_color.size());
                }
                
                
            }
            else{
                color = data_->minCol;
            }
            std::stack<label_expr::ptr> nodes;
            nodes.push(labelExpr);
            if (nodes.size() == 0){
                size_t newId = data_->n_vert;
                data_->n_vert += 1;
                data_->n_edge += 2;
                data_->vertices.insert(newId);
                data_->vert_id[newId] = 2;
                data_->edges[stateId].insert(newId);
                data_->edges[newId].insert(nbrId);
                data_->pre_edges[nbrId].insert(newId);
                data_->colors.insert({newId,color});
                for (size_t i = 0; i < all_color.size(); i++){
                    data_->all_colors[i].insert({newId,all_color[i]});
                }
                
                std::vector<size_t> edgeLabel(data_->ap_id.size(),2);
                data_->labels.insert({newId, edgeLabel});
            }else{
                std::stack<label_expr::ptr> conjunct_only_nodes;
                while (nodes.size() != 0) {
                    label_expr::ptr curr_node = nodes.top();
                    nodes.pop();
                    if (curr_node->isOR()) {
                        nodes.push(curr_node->getLeft());
                        nodes.push(curr_node->getRight());
                    } else {
                        conjunct_only_nodes.push(curr_node);
                    }
                }
                
                while (conjunct_only_nodes.size() != 0) {
                    std::stack<label_expr::ptr> single_node;
                    single_node.push(conjunct_only_nodes.top());
                    conjunct_only_nodes.pop();

                    size_t newId = data_->n_vert;
                    data_->n_vert += 1;
                    data_->n_edge += 2;
                    data_->vertices.insert(newId);
                    data_->vert_id[newId] = 2;
                    data_->edges[stateId].insert(newId);
                    data_->edges[newId].insert(nbrId);
                    data_->pre_edges[nbrId].insert(newId);
                    data_->colors.insert({newId,color});
                    for (size_t i = 0; i < all_color.size(); i++){
                        data_->all_colors[i].insert({newId,all_color[i]});
                    }

                    std::vector<size_t> edgeLabel(data_->ap_id.size(),2);
                    while (single_node.size() != 0) {
                        label_expr::ptr curr_node = single_node.top();
                        single_node.pop();
                        if (curr_node->isAND()) {
                            single_node.push(curr_node->getLeft());
                            single_node.push(curr_node->getRight());
                        } else if (curr_node->isNOT()){
                            unsigned int atom_id = curr_node->getLeft()->getAtom().getAPIndex();
                            edgeLabel[atom_id] = 0;
                        } else if (curr_node->isAtom()){
                            unsigned int atom_id = curr_node->getAtom().getAPIndex();
                            edgeLabel[atom_id] = 1;
                        }
                    }
                    data_->labels[newId] = edgeLabel;
                }
            }
        }

        virtual void notifyEndOfState(unsigned int stateId) override {
            UNUSED(stateId);
        }

        virtual void notifyEnd() override {}

        virtual void notifyAbort() override {}

        virtual void notifyWarning(const std::string &warning) override {
            std::cerr << "Warning: " << warning << std::endl;
        }

    private:
        /** Reference to the output stream */
        std::ostream &out;
    };

}// namespace cpphoafparser

#endif