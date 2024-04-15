#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
using namespace std;

class Solution {
public:
  int findChampion(vector<vector<int>> &grid) {
    for (size_t i = 0; i < grid.size(); i++) {
      auto sum = std::accumulate(grid[i].begin(), grid[i].end(), 0);
      if (sum == grid.size() - 1)
        return i;
    }
    return 0;
  }
};

TEST_CASE("2923. Find Champion I", "[2923]") {
  Solution s;
  std::vector<std::vector<int>> in = {{0, 0, 1}, {1, 0, 1}, {0, 0, 0}};
  int ans = 1;
  REQUIRE(s.findChampion(in) == ans);
}
