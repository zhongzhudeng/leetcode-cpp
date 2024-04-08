#include <vector>

#include <catch2/catch_test_macros.hpp>
using namespace std;

class Solution {
public:
  int isValid(vector<int> &nums) { return 1; }
};

TEST_CASE("0. template", "[template]") {
  Solution s;
  vector<int> in = {1, 2, 3};
  int ans = 1;
  REQUIRE(s.isValid(in) == ans);
}