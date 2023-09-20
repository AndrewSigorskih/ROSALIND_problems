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
    std::map<char, size_t> counts;
    for (const char c: text)
    {
        vec.push_back({c, counts[c]});
        counts[c]++;
    }
}

std::string sort_string(const std::string& text)
{
    std::string result = text;
    std::sort(result.begin(), result.end());
    return result;
}

std::string inverse_bwt(const std::string& transform)
{
    size_t n = transform.size();
    std::string sort_seq = sort_string(transform);
    
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
    std::string sort_seq = sort_string(transform);

    std::vector<position> first, last;
    first.reserve(n);
    last.reserve(n);
    index_seq(sort_seq, first);
    index_seq(transform, last);

    return static_cast<size_t>(
        std::find(first.begin(), first.end(), last[ind]) - first.begin()
    );
}

void _fillFirstColumn(const std::string& transform, std::map<char, size_t>& firstOccurence)
{
    if (firstOccurence.size() > 0)
        firstOccurence.clear();
    
    std::string firstColumn = sort_string(transform);
    for (const char c: firstColumn)
    {
        if (firstOccurence.count(c) == 0)
        {
            firstOccurence[c] =\
                std::find(firstColumn.begin(), firstColumn.end(), c) - firstColumn.begin();
        }
    }
}

BWMatcher::BWMatcher(const std::string& transform)
{
    this->lastColumn = transform;
    size_t n = transform.size();
    std::string sorted_transform = sort_string(transform);

    this->first.reserve(n);
    this->last.reserve(n);
    index_seq(sorted_transform, this->first);
    index_seq(transform, this->last);
}

#pragma GCC diagnostic ignored "-Wreturn-type"
size_t BWMatcher::match(const std::string& pattern)
{
    size_t top = 0, bottom = this->lastColumn.size() - 1;
    std::vector<char> pat(pattern.begin(), pattern.end());

    while (top <= bottom)
    {
        if (pat.size() > 0)
        {
            char symbol = pat.back();
            pat.pop_back();

            std::vector<char> positions(this->lastColumn.begin() + top,
                                        this->lastColumn.begin() + bottom + 1);
            
            auto found = std::find(positions.begin(), positions.end(), symbol);
            if (found != positions.end())
            {
                size_t topIndex = found - positions.begin() + top;
                size_t bottomIndex = std::distance(std::find(positions.rbegin(),
                                                             positions.rend(),
                                                             symbol),
                                                   positions.rend()) - 1 + top;
                top = std::find(first.begin(), first.end(), last[topIndex]) - first.begin();
                bottom = std::find(first.begin(), first.end(), last[bottomIndex]) - first.begin();
            } else {
                return 0;
            }
        } else {
            return bottom - top + 1;
        }
    }
}
#pragma GCC diagnostic pop

CountSymbols::CountSymbols(const std::string& transform)
{
    this->counts.reserve(transform.size());
    this->counts.emplace_back(std::map<char, size_t>());
    for (size_t i = 0; i < transform.size(); ++i)
    {
        this->counts.push_back(this->counts[i]);
        this->counts[i+1][transform[i]] += 1;
    }
}

BetterBWMatcher::BetterBWMatcher(const std::string& transform)
{
    this->lastColumn = transform;
    // first Occurence of symbols in firstCol
    _fillFirstColumn(transform, this->firstOccurence);
    // symbol counts
    this->symbolCounts = CountSymbols(transform);
}

#pragma GCC diagnostic ignored "-Wreturn-type"
size_t BetterBWMatcher::match(const std::string& pattern)
{
    size_t top = 0, bottom = this->lastColumn.size() - 1;
    std::vector<char> pat(pattern.begin(), pattern.end());

    while (top <= bottom)
    {
        if (pat.size() > 0)
        {
            char symbol = pat.back();
            pat.pop_back();

            std::vector<char> positions(this->lastColumn.begin() + top,
                                        this->lastColumn.begin() + bottom + 1);
            
            auto found = std::find(positions.begin(), positions.end(), symbol);
            if (found != positions.end())
            {
                top = firstOccurence[symbol] + symbolCounts.get_count(top, symbol);
                bottom = firstOccurence[symbol] + symbolCounts.get_count(bottom+1, symbol) - 1;
            } else { 
                return 0;
            }
        } else {
            return bottom - top + 1;
        }
    }
}
#pragma GCC diagnostic pop

PatternMatcher::PatternMatcher(const std::string& text)
{
    this->text = text;
    this->sarray = suffix_array(text).get_array();
    this->bwt = BWTfromArray(text, this->sarray);
    this->symbolCounts = CountSymbols(this->bwt);
    _fillFirstColumn(this->bwt, this->firstOccurence);
}

void PatternMatcher::match(const std::string& pattern, std::set<size_t>& result)
{
    size_t top = 0, bottom = this->bwt.size() - 1;
    ssize_t currIndex = pattern.size() - 1;

    while (top <= bottom)
    {
        if (currIndex >= 0)
        {
            char symbol = pattern[currIndex--];
            if (symbolCounts.get_count(bottom+1, symbol) > symbolCounts.get_count(top, symbol))
            {
                top = firstOccurence[symbol] + symbolCounts.get_count(top, symbol);
                bottom = firstOccurence[symbol] + symbolCounts.get_count(bottom+1, symbol) - 1;
            } else {
                return;
            }
        } else {
            for (size_t i = top; i <= bottom; ++i)
                result.insert(sarray[i]);
            return;
        }
    }
}

bool _matchesWithMismatches(const std::string& text,
                            size_t start,
                            const std::string& pattern,
                            size_t d)
{
    size_t mismatches = 0;
    for (size_t i = 0; i < pattern.size(); ++i)
    {
        if (pattern[i] != text[start+i])
            ++mismatches;
        if (mismatches > d)
            return false;
    }
    return true;
}

void PatternMatcher::match_approx(const std::string& pattern, std::vector<size_t>& result, size_t d)
{
    size_t n = pattern.size();
    size_t k = n / (d + 1);

    std::vector<std::pair<std::string, size_t>> seeds;
    seeds.reserve(d+1);
    for (size_t i = 0; i < d; ++i)
    {
        seeds.push_back( { pattern.substr(i*k, k), i*k } );
    }
    seeds.push_back( { pattern.substr(d*k, n-d*k+1), d*k } );

    std::set<size_t> posMatches;
    for (const auto& seed: seeds)
    {
        size_t top = 0, bottom = this->text.size() - 1;
        ssize_t currIndex = seed.first.size() - 1;
        while (top <= bottom)
        {
            if (currIndex >= 0)
            {
                char symbol = seed.first[currIndex--];
                if (symbolCounts.get_count(bottom+1, symbol) > symbolCounts.get_count(top, symbol))
                {
                    top = firstOccurence[symbol] + symbolCounts.get_count(top, symbol);
                    bottom = firstOccurence[symbol] + symbolCounts.get_count(bottom+1, symbol) - 1;
                } else {
                    break;
                }
            } else {
                for (size_t i = top; i <= bottom; ++i)
                    posMatches.insert(sarray[i] - seed.second);
                break;
            }
        }
    }
    for (auto val: posMatches)
    {
        if (_matchesWithMismatches(text, val, pattern, d))
        {
            result.push_back(val);
        }
    }
}
