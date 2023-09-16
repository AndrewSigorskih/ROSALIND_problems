#pragma once
#include <unordered_map>
#include "suffix_array.hpp"

using position = std::pair<char, size_t>;

std::string BWTfromArray(const std::string& text, const std::vector<size_t> array);

void index_seq(const std::string& text, std::vector<position>& vec);
std::string inverse_bwt(const std::string& transform);

size_t last2first(const std::string& transform, size_t ind);
