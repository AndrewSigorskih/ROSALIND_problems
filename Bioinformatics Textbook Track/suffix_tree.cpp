#include "suffix_tree.hpp"

_suffix_tree_edge* _suffix_tree_node::get_edge(char c)
{
    if (this->edges.find(c) == this->edges.end())
        return this->default_map_value;
    return this->edges[c];
}

suffix_tree::suffix_tree(const std::string& text)
{
    _text = text;
    infty = _text.size();
    joker_edge = new _suffix_tree_edge(0, 1, nullptr);
    joker = new _suffix_tree_node();
    root = new _suffix_tree_node(joker);
    joker->set_default(joker_edge);
    joker_edge->target_node = root;
    std::ignore = this->build_tree(root, 0, this->infty);
}

suffix_tree::~suffix_tree()
{
    destroy_node(root);
    destroy_node(joker);
    delete joker_edge;
}

_suffix_tree_node* suffix_tree::build_tree(_suffix_tree_node* node, size_t n, size_t size, size_t skip)
{
    while (n < size)
    {
        char c = this->_text[n];
        _suffix_tree_edge* edge = node->get_edge(c);
        if (edge)
        {
            size_t first = edge->start_pos, last = edge->end_pos;
            size_t i = first, n0 = n;
            if (skip > 0)
            {
                size_t can_skip = std::min(skip, last-first);
                i += can_skip;
                n += can_skip;
                skip -= can_skip;
            }

            while ((i < last) && \
                   (n < size) && \
                   (
                       (this->_text[i] == this->_text[n]) ||
                       (edge == this->joker_edge)
                   )) 
            {
                i++;
                n++;
            }

            if (i == last)
            { // go to the next node
                node = edge->target_node;
            } else { // splitting edge
                _suffix_tree_node* middle_node = new _suffix_tree_node();
                middle_node->set_parent(node);
                middle_node->edges[this->_text[i]] = edge;
                node->edges[c] = new _suffix_tree_edge(first, i, middle_node);
                edge->start_pos = i;
                edge->target_node->set_parent(middle_node);
                middle_node->suff_link = this->build_tree(node->suff_link, n0, n, i-first);
                node = middle_node;
            }

        } else {
            // nowhere to go; creating a new leaf
            _suffix_tree_node* new_leaf = new _suffix_tree_node;
            new_leaf->set_parent(node);
            node->edges[c] = new _suffix_tree_edge(n, this->infty, new_leaf);
            node = node->suff_link;
        }
    }
    return node;
}

void suffix_tree::destroy_node(_suffix_tree_node *node)
{
    if (!node) return;
    for (auto [key, edge]: node->edges)
    {
        if (edge->target_node)
            destroy_node(edge->target_node);
        delete edge;
    }
    delete node;
    node = nullptr;
}

void suffix_tree::get_edges(std::ofstream& ost, _suffix_tree_node* node)
{
    if (!node)
        node = this->root;
    for (const auto [key, edge]: node->edges)
    {
        size_t substr_len = edge->length();
        ost.write(this->_text.c_str()+edge->start_pos, substr_len);
        ost << '\n';
        this->get_edges(ost, edge->target_node);
    }
}

void suffix_tree::lrep(std::ofstream& ost)
{
    size_t maxHeight = 0, substrStart = 0;
    this->lrep_traverse(root, 0, maxHeight, substrStart);
    ost.write(this->_text.c_str()+substrStart, maxHeight);
    ost << '\n';
}

void suffix_tree::lrep_traverse(_suffix_tree_node* node,
                                size_t currHeight,
                                size_t& maxHeight,
                                size_t& suffixIndex)
{
    if (!node) return;

    for (const auto [key, edge]: node->edges)
    {
        size_t substr_len = edge->length();

        if (edge->target_node->isTerminal())
        {
            if (maxHeight < currHeight)
            {
                maxHeight = currHeight;
                suffixIndex = this->_text.size() - currHeight - substr_len;
            }
        } else {
            lrep_traverse(edge->target_node, currHeight + substr_len, maxHeight, suffixIndex);
        }
    }
}
