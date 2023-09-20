#pragma once
#include <map>
#include <set>
#include "suffix_array.hpp"

using position = std::pair<char, size_t>;

std::string BWTfromArray(const std::string& text, const std::vector<size_t> array);

void index_seq(const std::string& text, std::vector<position>& vec);
std::string sort_string(const std::string& text);
std::string inverse_bwt(const std::string& transform);

size_t last2first(const std::string& transform, size_t ind);

class BWMatcher
{
public:
    BWMatcher(const std::string& transform);
    size_t match(const std::string& pattern);
private:
    std::string lastColumn;
    std::vector<position> first, last;
};

class CountSymbols
{
public:
    CountSymbols() {}
    CountSymbols(const std::string& transform);
    size_t get_count(size_t ind, char symbol) { return counts[ind][symbol]; }
private:
    std::vector<std::map<char, size_t>> counts;
};

class BetterBWMatcher
{
public:
    BetterBWMatcher(const std::string& transform);
    size_t match(const std::string& pattern);
private:
    std::string lastColumn;
    CountSymbols symbolCounts;
    std::map<char, size_t> firstOccurence;
};

class PatternMatcher
{
public:
    PatternMatcher(const std::string& text);
    void match(const std::string& pattern, std::set<size_t>& result);
    void match_approx(const std::string& pattern, std::vector<size_t>& result, size_t d);
private:
    std::string text;
    std::string bwt;
    std::vector<size_t> sarray;
    CountSymbols symbolCounts;
    std::map<char, size_t> firstOccurence;
};
