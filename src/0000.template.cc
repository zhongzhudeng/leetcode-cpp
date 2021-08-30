#include <vector>

#include "catch2/catch.hpp"
using namespace std;

class Solution {
 public:
  int isValid(vector<int> &nums) { return 1; }
};

TEST_CASE("template") {
  Solution s;
  vector<int> in = {1, 2, 3};
  int ans = 1;
  REQUIRE(s.isValid(in) == ans);
}