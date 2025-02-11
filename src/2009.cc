#include <algorithm>
#include <unordered_set>
#include <vector>
using namespace std;

class Solution {
public:
  int minOperations(vector<int> &nums) {
    std::unordered_set<int> unique(nums.begin(), nums.end());
    std::vector<int> unique_sorted(unique.begin(), unique.end());
    std::sort(unique_sorted.begin(), unique_sorted.end());
    int n = nums.size();
    int res = n;
    int j = 0;
    for (int i = 0; i < unique_sorted.size(); i++) {
      int left = unique_sorted[i];
      int right = left + n - 1;
      while (j < unique_sorted.size() and unique_sorted[j] <= right) {
        res = std::min(res, n - (j - i + 1));
        j++;
      }
    }
    return res;
  }
};

#include <catch2/catch_test_macros.hpp>

TEST_CASE("2009. Minimum Number of Operations to Make Array Continuous",
          "[2009]") {
  Solution s;
  vector<int> in = {4, 2, 5, 3};
  int ans = 0;
  REQUIRE(s.minOperations(in) == ans);

  in = {1, 2, 3, 5, 6};
  ans = 1;
  REQUIRE(s.minOperations(in) == ans);

  in = {1, 10, 100, 1000};
  ans = 3;
  REQUIRE(s.minOperations(in) == ans);
}