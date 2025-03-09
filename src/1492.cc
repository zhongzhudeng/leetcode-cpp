#include <vector>

using namespace std;

class Solution {
public:
  int kthFactor(int n, int k) {
    vector<int> facts;
    for (int i = 1; i <= n; i++)
      if (n % i == 0)
        facts.push_back(i);
    return k <= facts.size() ? facts[k - 1] : -1;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1492. The kth Factor of n", "[1492]") {
  Solution s;
  REQUIRE(s.kthFactor(12, 3) == 3);
  REQUIRE(s.kthFactor(7, 2) == 7);
  REQUIRE(s.kthFactor(4, 4) == -1);
}