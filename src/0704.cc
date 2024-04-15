#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <vector>

namespace std {

class Solution {
public:
  int search(vector<int> &nums, int target) {
    auto it = lower_bound(nums.begin(), nums.end(), target);
    return (it != nums.end() && *it == target) ? distance(nums.begin(), it)
                                               : -1;
  }
};

TEST_CASE("704. Binary Search","[0704]") {
  Solution s;
  vector<int> nums = {-1, 0, 3, 5, 9, 12};
  int target = 9;
  int output = 4;
  REQUIRE(s.search(nums, target) == output);
  nums = {-1, 0, 3, 5, 9, 12};
  target = 2;
  output = -1;
  REQUIRE(s.search(nums, target) == output);
  nums = {-1};
  target = 2;
  output = -1;
  REQUIRE(s.search(nums, target) == output);
}
} // namespace std
