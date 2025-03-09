#include <algorithm>
#include <vector>
using namespace std;

class Solution {
public:
  int minimumSize(vector<int> &nums, int maxOperations) {
    int left = 1, right = *max_element(nums.begin(), nums.end());
    int ans = 0;
    while (left <= right) {
      int y = (left + right) / 2;
      long long ops = 0;
      for (int x : nums) {
        ops += (x - 1) / y;
      }
      if (ops <= maxOperations) {
        ans = y;
        right = y - 1;
      } else {
        left = y + 1;
      }
    }
    return ans;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1760. Minimum Limit of Balls in a Bag", "[1760]") {
  Solution s;
  vector<int> nums = {9};
  int maxOperations = 2;
  REQUIRE(s.minimumSize(nums, maxOperations) == 3);

  nums = {2, 4, 8, 2};
  maxOperations = 4;
  REQUIRE(s.minimumSize(nums, maxOperations) == 2);
}