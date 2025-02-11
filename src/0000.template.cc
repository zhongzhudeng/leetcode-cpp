#include <vector>
using namespace std;
class Solution {

public:
  int isValid(vector<int> &nums) { return 1; }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("0. template", "[template]") {
  Solution s;
  vector<int> in = {1, 2, 3};
  int ans = 1;
  REQUIRE(s.isValid(in) == ans);
}