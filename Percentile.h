//helper class for finding percentiles
#pragma once

#include <vector>
#include <algorithm>
#include <functional>

template<class T>
class IterativePercentile {
public:
  /// Percentile has to be in range [0, 1(
  IterativePercentile(double percentile)
    : _percentile(percentile)
  { }

  // Adds a number in O(log(n))
  void add(const T& x) {
    if (_lower.empty() || x <= _lower.front()) {
      _lower.push_back(x);
      std::push_heap(_lower.begin(), _lower.end(), std::less<T>());
    } else {
      _upper.push_back(x);
      std::push_heap(_upper.begin(), _upper.end(), std::greater<T>());
    }

    unsigned size_lower = (unsigned)((_lower.size() + _upper.size()) * _percentile) + 1;
    if (_lower.size() > size_lower) {
      // lower to upper
      std::pop_heap(_lower.begin(), _lower.end(), std::less<T>());
      _upper.push_back(_lower.back());
      std::push_heap(_upper.begin(), _upper.end(), std::greater<T>());
      _lower.pop_back();
    } else if (_lower.size() < size_lower) {
      // upper to lower
      std::pop_heap(_upper.begin(), _upper.end(), std::greater<T>());
      _lower.push_back(_upper.back());
      std::push_heap(_lower.begin(), _lower.end(), std::less<T>());
      _upper.pop_back();
    }            
  }

  /// Access the percentile in O(1)
  const T& get() const {
    return _lower.front();
  }

  void clear() {
    _lower.clear();
    _upper.clear();
  }

private:
  double _percentile;
  std::vector<T> _lower;
  std::vector<T> _upper;
};