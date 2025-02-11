
#include <vector>
using namespace std;

class Solution {
public:
  int minMalwareSpread(vector<vector<int>> &graph, vector<int> &initial) {}
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("928. Minimize Malware Spread II", "[0928]") {
  Solution s;

  std::vector<std::vector<int>> graph = {{1, 1, 0}, {1, 1, 0}, {0, 0, 1}};
  std::vector<int> initial = {0, 1};
  int ans = 0;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 1, 0}, {1, 1, 1}, {0, 1, 1}};
  initial = {0, 1};
  ans = 1;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 1, 0, 0}, {1, 1, 1, 0}, {0, 1, 1, 1}, {0, 0, 1, 1}};
  initial = {0, 1};
  ans = 1;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);
}
