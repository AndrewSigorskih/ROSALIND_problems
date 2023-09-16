#include "bwt.hpp"

std::string BWTfromArray(const std::string& text, const std::vector<size_t> array)
{
    size_t n = text.size();
    std::vector<char> bwt(n);

    for(size_t i = 0; i < n; ++i)
    {
        size_t ind = (array[i] + n - 1) % n;
        bwt[i] = text[ind];
    }

    return std::string(bwt.begin(), bwt.end());
}

void index_seq(const std::string& text, std::vector<position>& vec)
{
    std::unordered_map<char, size_t> counts;
    for (const char c: text)
    {
        vec.push_back({c, counts[c]});
        counts[c]++;
    }
}

std::string inverse_bwt(const std::string& transform)
{
    size_t n = transform.size();
    std::string sort_seq = transform;
    std::sort(sort_seq.begin(), sort_seq.end());

    std::vector<position> first, last;
    first.reserve(n);
    last.reserve(n);
    index_seq(sort_seq, first);
    index_seq(transform, last);
    
    std::vector<char> result;
    result.reserve(n);
    position current = {'$', 0};

    for (size_t i = 0; i < n; ++i)
    {
        auto ind = (std::find(last.begin(), last.end(), current) - last.begin());
        current = first[ind];
        result.push_back(current.first);
    }
    return std::string(result.begin(), result.end());
}

size_t last2first(const std::string& transform, size_t ind)
{
    size_t n = transform.size();
    std::string sort_seq = transform;
    std::sort(sort_seq.begin(), sort_seq.end());

    std::vector<position> first, last;
    first.reserve(n);
    last.reserve(n);
    index_seq(sort_seq, first);
    index_seq(transform, last);

    return static_cast<size_t>(
        std::find(first.begin(), first.end(), last[ind]) - first.begin()
    );
}
