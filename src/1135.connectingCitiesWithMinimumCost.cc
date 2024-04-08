#include <algorithm>
#include <vector>

#include <catch2/catch_all.hpp>

#include "union_find.hpp"
using namespace std;
class Solution {
 public:
  int minimumCost(int n, vector<vector<int>>& connections) {
    UnionFind uf(n);
    int ans = 0, edgeCount = 0;
    sort(connections.begin(), connections.end(),
         [](const auto& u, const auto& v) { return u[2] < v[2]; });
    for (auto& e : connections) {
      int c1 = e[0] - 1, c2 = e[1] - 1, w = e[2];
      if (uf.find(c1) != uf.find(c2)) {
        uf.unite(c1, c2);
        ans += w;
        edgeCount++;
      }
    }

    if (edgeCount != n - 1) {
      return -1;
    }

    return ans;
  }
};

TEST_CASE("1135. Connecting Cities With Minimum Cost") {
  int n = 3;
  vector<vector<int>> connections = {{1, 2, 5}, {1, 3, 6}, {2, 3, 1}};
  int ans = 6;
  Solution s;
  REQUIRE(s.minimumCost(n, connections) == ans);
  n = 4;
  connections = {{1, 2, 3}, {3, 4, 4}};
  ans = -1;
  REQUIRE(s.minimumCost(n, connections) == ans);
}