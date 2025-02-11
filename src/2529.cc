#include <algorithm>

#include <vector>
using namespace std;

class Solution {
public:
  int maximumCount(vector<int> &nums) {
    auto lo = std::lower_bound(nums.begin(), nums.end(), 0) - nums.begin();
    auto hi = nums.end() - std::lower_bound(nums.begin(), nums.end(), 1);
    return std::max(lo, hi);
  }
};

#include <catch2/catch_test_macros.hpp>

TEST_CASE("2529. Maximum Count of Positive Integer and Negative Integer",
          "[2529]") {
  Solution s;
  vector<int> in = {-2, -1, -1, 1, 2, 3};
  int ans = 3;
  REQUIRE(s.maximumCount(in) == ans);

  in = {-3, -2, -1, 0, 0, 1, 2};
  ans = 3;
  REQUIRE(s.maximumCount(in) == ans);

  in = {5, 20, 66, 1314};
  ans = 4;
  REQUIRE(s.maximumCount(in) == ans);
}