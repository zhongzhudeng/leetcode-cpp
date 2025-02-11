
#include <vector>
using namespace std;

class Solution {
public:
  int minMalwareSpread(vector<vector<int>> &graph, vector<int> &initial) {}
};

#include <catch2/catch_test_macros.hpp>

TEST_CASE("924. Minimize Malware Spread", "[0924]") {
  Solution s;

  std::vector<std::vector<int>> graph = {{1, 1, 0}, {1, 1, 0}, {0, 0, 1}};
  std::vector<int> initial = {0, 1};
  int ans = 0;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  initial = {0, 2};
  ans = 0;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  initial = {1, 2};
  ans = 1;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);
}
