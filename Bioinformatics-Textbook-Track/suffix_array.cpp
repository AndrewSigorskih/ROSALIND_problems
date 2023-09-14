#include "suffix_array.hpp"
#include <cstring>

suffix_array::suffix_array(const std::string& text)
{
    this->text = text;
    this->build_array();
}

void suffix_array::build_array()
{
    size_t n = this->text.size();
    this->array.clear();
    suffix suffixes[n];

    for (size_t i = 0; i < n; ++i)
    {
        suffixes[i].index = i;
        suffixes[i].rank[0] = text[i] - 'a';
        suffixes[i].rank[1] = ((i+1) < n) ? (text[i+1] - 'a') : -1;
    }

    std::sort(suffixes, suffixes+n);

    size_t ind[n];
    for (size_t k = 4; k < 2*n; k = k*2)
    {
        size_t rank = 0;
        size_t prev_rank = suffixes[0].rank[0];
        suffixes[0].rank[0] = rank;
        ind[suffixes[0].index] = 0;

        for (size_t i = 1; i < n; ++i)
        {
            if (suffixes[i].rank[0] == prev_rank &&
                    suffixes[i].rank[1] == suffixes[i-1].rank[1])
            {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = rank;
            } else {
                prev_rank = suffixes[i].rank[0];
                suffixes[i].rank[0] = ++rank;
            }
            ind[suffixes[i].index] = i;
        }

        for (size_t i = 0; i < n; i++)
        {
            size_t nextindex = suffixes[i].index + k/2;
            suffixes[i].rank[1] = (nextindex < n) ?
                                  suffixes[ind[nextindex]].rank[0] : -1;
        }
        std::sort(suffixes, suffixes+n);
    }

    this->array.reserve(n);
    for (size_t i = 0; i < n; ++i)
        this->array.push_back(suffixes[i].index);
}

void suffix_array::print_array(std::ofstream& ost)
{
    size_t n = this->array.size();
    for (size_t i = 0; i < n - 1; ++i)
        ost << this->array[i] << ", ";
    ost << this->array[n-1] << "\n";
}

void suffix_array::all_matches(const std::string& pat, std::vector<size_t>& result)
{
    size_t minIdx = 0, maxIdx = this->text.size();

    while (minIdx < maxIdx)
    {
        size_t midIdx = (minIdx + maxIdx) >> 1;
        if (
            strncmp(
                pat.c_str(), 
                this->text.c_str() + this->array[midIdx],
                pat.size()
            ) > 0
        ) {
            minIdx = midIdx + 1;
        } else {
            maxIdx = midIdx;
        }
    }
    size_t first = minIdx;
    maxIdx = this->text.size();
    while (minIdx < maxIdx)
    {
        size_t midIdx = (minIdx + maxIdx) >> 1;
        if (
            strncmp(
                pat.c_str(), 
                this->text.c_str() + this->array[midIdx],
                pat.size()
            ) < 0
        ) {
            maxIdx = midIdx;
        } else {
            minIdx = midIdx + 1;
        }
    }
    size_t last = maxIdx;
    if (first > last)
    {
        return;
    } else {
        for (size_t i = first; i < last; ++i)
            result.push_back(this->array[i]);
    }
}
