#include <algorithm>
#include <unordered_map>
#include <vector>

#include <catch2/catch_all.hpp>

#include "tarjan.hpp"
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

    UnionFind uf(n);
    vector<vector<int>> ans(2);
    vector<int> label(m);
    for (int i = 0; i < m;) {
      // 找出所有权值为 w 的边，下标范围为 [i, j)
      int w = edges[i][2];
      int j = i;
      while (j < m && edges[j][2] == edges[i][2]) {
        ++j;
      }

      // 存储每个连通分量在图 G 中的编号
      unordered_map<int, int> compToId;
      // 图 G 的节点数
      int gn = 0;

      for (int k = i; k < j; ++k) {
        int x = uf.find(edges[k][0]);
        int y = uf.find(edges[k][1]);
        if (x != y) {
          if (!compToId.count(x)) {
            compToId[x] = gn++;
          }
          if (!compToId.count(y)) {
            compToId[y] = gn++;
          }
        } else {
          // 将自环边标记为 -1
          label[edges[k][3]] = -1;
        }
      }

      // 图 G 的边
      vector<vector<int>> gm(gn), gmid(gn);

      for (int k = i; k < j; ++k) {
        int x = uf.find(edges[k][0]);
        int y = uf.find(edges[k][1]);
        if (x != y) {
          int idx = compToId[x], idy = compToId[y];
          gm[idx].push_back(idy);
          gmid[idx].push_back(edges[k][3]);
          gm[idy].push_back(idx);
          gmid[idy].push_back(edges[k][3]);
        }
      }

      vector<int> bridges = TarjanRecursive(gn, gm, gmid).getCuttingEdge();
      // 将桥边（关键边）标记为 1
      for (int id : bridges) {
        ans[0].push_back(id);
        label[id] = 1;
      }

      for (int k = i; k < j; ++k) {
        uf.unite(edges[k][0], edges[k][1]);
      }

      i = j;
    }

    // 未标记的边即为非桥边（伪关键边）
    for (int i = 0; i < m; ++i) {
      if (!label[i]) {
        ans[1].push_back(i);
      }
    }

    for (auto& v : ans) {
      sort(v.begin(), v.end());
    }
    sort(ans.begin(), ans.end());

    return ans;
  }
};

TEST_CASE(
    "1489. Find Critical and Pseudo-Critical Edges in Minimum Spanning Tree") {
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