// suffix array class
// construction time complexity: O(n*logn*logn) implementation
// source: https://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm/
#pragma once
#include <algorithm> // sort
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct suffix
{
    size_t index;
    size_t rank[2];
    bool operator < (const suffix& other)
    {
        return (this->rank[0] == other.rank[0]) ? (this->rank[1] < other.rank[1] ? 1: 0) :
               (this->rank[0] < other.rank[0] ? 1: 0);
    }
};

class suffix_array
{
public:
    suffix_array(const std::string&);
    void print_array(std::ofstream&);
    void all_matches(const std::string&, std::vector<size_t>&);
    const std::vector<size_t>& get_array() const { return this->array; }
private:
    void build_array();
private:
    std::vector<size_t> array;
    std::string text;
};
