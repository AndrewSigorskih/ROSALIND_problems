#include <iostream>
#include <iomanip> // precision
#include <fstream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm> // min element

using std::vector;
using std::unordered_map;

struct Edge
{
    int end;
    double weight;
};

// hash function for storing pair as key in unordered map
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, should use sth. like boost.hash_combine
        return h1 ^ h2;  
    }
};
// function to find key for minimal value in hasmap
template <class Key, class Value, class HashFn>
std::pair<Key, Value> findMinValuePair(
    std::unordered_map<Key, Value, HashFn> &x)
{
    return *std::min_element(x.begin(), x.end(),
                             [](const std::pair<Key, Value> &p1,
                                const std::pair<Key, Value> &p2)
                             {
                                 return p1.second < p2.second;
                             });
}

using key = std::pair<int, int>;
using distMatrix = unordered_map<key, double, pair_hash>;

// to store only one triangle from symmetrical matrix
inline key orderPair(int const a, int const b)
{
    return (a < b) ? key(a, b) : key(b, a);
}