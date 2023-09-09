#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct _suffix_tree_node;
struct _suffix_tree_edge {
    size_t start_pos, end_pos;
    _suffix_tree_node* target_node;
    _suffix_tree_edge(size_t i, size_t j, _suffix_tree_node* target) : 
        start_pos(i), end_pos(j), target_node(target) {}
    size_t length() { return end_pos - start_pos; }
};

using _edges_dict = std::unordered_map<char, _suffix_tree_edge*>;

enum lcsLabel
{
    UNKNOWN = 0,
    FIRST,
    SECOND,
    BOTH
};

struct lcsHelper
{
    size_t firstSize, maxHeight = 0, substrStart = 0;
    std::unordered_map<_suffix_tree_node*, lcsLabel> map;
};

struct _suffix_tree_node {
    _edges_dict edges;
    _suffix_tree_node* parent = nullptr;
    _suffix_tree_node* suff_link = nullptr;
    _suffix_tree_edge* default_map_value = nullptr;

    _suffix_tree_node() {}
    _suffix_tree_node(_suffix_tree_node* suff_link_to_set) { this->suff_link = suff_link_to_set; }
    _suffix_tree_edge* get_edge(char c);
    void set_default(_suffix_tree_edge* edge) { this->default_map_value = edge; }
    void set_parent(_suffix_tree_node* node) { this->parent = node; }
    bool isTerminal() { return this->edges.size() == 0; }
};

class suffix_tree {
public:
    suffix_tree(const std::string&);
    ~suffix_tree();
    void get_edges(std::ofstream& ost, _suffix_tree_node* node = nullptr);
    void lrep(std::ofstream& ost);
    void lcs(std::ofstream& ost, size_t first_size);
private:
    _suffix_tree_node* build_tree(_suffix_tree_node* node, size_t n, size_t size, size_t skip = 0);
    void lrep_traverse(_suffix_tree_node*, size_t, size_t&, size_t&);
    lcsLabel lcs_traverse(_suffix_tree_node*, size_t, lcsHelper&);
    void destroy_node(_suffix_tree_node* node);
private:
    size_t infty;
    std::string _text;
    _suffix_tree_edge* joker_edge;
    _suffix_tree_node* joker;
    _suffix_tree_node* root;
};

