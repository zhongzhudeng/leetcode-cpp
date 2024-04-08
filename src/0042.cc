#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <vector>
using namespace std;

class Solution {
public:
  int trap1(vector<int> &height) {
    int ans = 0;
    for (int i = 0; i < height.size(); i++) {
      int left_max = 0, right_max = 0;
      for (int left = 0; left < i; left++) {
        left_max = std::max(height[left], left_max);
      }
      for (int right = height.size() - 1; right > i; right--) {
        right_max = std::max(height[right], right_max);
      }
      ans += std::max(0, std::min(left_max, right_max) - height[i]);
    }
    return ans;
  }

  int trap(vector<int> &height) {
    int ans = 0;
    size_t left = 0, right = height.size() - 1;
    int left_max = 0, right_max = 0;
    while (left < right) {
      if (height[left] < height[right]) {
        height[left] < left_max ? ans += left_max - height[left]
                                : left_max = height[left];
        left++;
      } else {
        height[right] < right_max ? ans += right_max - height[right]
                                  : right_max = height[right];
        right--;
      }
    }
    return ans;
  }
};

TEST_CASE("42. Trapping Rain Water","[42]") {
  Solution s;
  vector<int> in = {0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 2, 1};
  int ans = 6;
  REQUIRE(s.trap1(in) == ans);
  REQUIRE(s.trap(in) == ans);

  in = {4, 2, 0, 3, 2, 5};
  ans = 9;
  REQUIRE(s.trap1(in) == ans);
  REQUIRE(s.trap(in) == ans);
}