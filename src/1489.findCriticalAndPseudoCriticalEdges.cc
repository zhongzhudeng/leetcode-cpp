#include <vector>

#include "catch2/catch.hpp"
#include "union_find.hpp"
using namespace std;

class Solution {
 public:
  vector<vector<int>> findCriticalAndPseudoCriticalEdges(
      int n, vector<vector<int>>& edges) {
    int m = edges.size();
    for (int i = 0; i < m; ++i) {
      edges[i].push_back(i);
    }
    sort(edges.begin(), edges.end(),
         [](const auto& u, const auto& v) { return u[2] < v[2]; });

    // 计算 value value: 权值
    UnionFind uf_std(n);
    int value = 0;
    for (int i = 0; i < m; ++i) {
      if (uf_std.unite(edges[i][0], edges[i][1])) {
        value += edges[i][2];
      }
    }

    vector<vector<int>> ans(2);

    for (int i = 0; i < m; ++i) {
      // 判断是否是关键边
      UnionFind uf(n);
      int v = 0;
      for (int j = 0; j < m; ++j) {
        if (i != j && uf.unite(edges[j][0], edges[j][1])) {
          v += edges[j][2];
        }
      }
      if (uf.setCount != 1 || (uf.setCount == 1 && v > value)) {
        ans[0].push_back(edges[i][3]);
        continue;
      }

      // 判断是否是伪关键边
      uf = UnionFind(n);
      uf.unite(edges[i][0], edges[i][1]);
      v = edges[i][2];
      for (int j = 0; j < m; ++j) {
        if (i != j && uf.unite(edges[j][0], edges[j][1])) {
          v += edges[j][2];
        }
      }
      if (v == value) {
        ans[1].push_back(edges[i][3]);
      }
    }

    return ans;
  }
};

TEST_CASE("1489") {
  Solution s;
  int n;
  vector<vector<int>> edges;
  vector<vector<int>> ans;

  SECTION("n=5") {
    n = 5;
    edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 2}, {0, 3, 2},
             {0, 4, 3}, {3, 4, 3}, {1, 4, 6}};
    ans = {{0, 1}, {2, 3, 4, 5}};
    REQUIRE(s.findCriticalAndPseudoCriticalEdges(n, edges) == ans);
  }
  SECTION("n=4") {
    n = 4;
    edges = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {0, 3, 1}};
    ans = {{}, {0, 1, 2, 3}};
    REQUIRE(s.findCriticalAndPseudoCriticalEdges(n, edges) == ans);
  }
}