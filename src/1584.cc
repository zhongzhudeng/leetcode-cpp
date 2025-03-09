#include "union_find.hpp"
#include "vector"
#include <algorithm>
#include <array>

using namespace std;
class Solution {
public:
  int minCostConnectPoints(vector<vector<int>> &points) {
    vector<std::array<int, 3>> edges;
    int len = static_cast<int>(points.size());
    UnionFind uf(len);
    int ans = 0, edgeCount = 0;
    edges.reserve(len * (len - 1) / 2);
    for (int i = 0; i < len - 1; i++) {
      for (int j = i + 1; j < len; j++) {
        edges.push_back({i, j,
                         std::abs(points[i][0] - points[j][0]) +
                             std::abs(points[i][1] - points[j][1])});
      }
    }
    sort(edges.begin(), edges.end(),
         [](const auto &u, const auto &v) { return u[2] < v[2]; });

    for (auto &e : edges) {
      if (edgeCount == len - 1) {
        return ans;
      }
      auto p1 = e[0], p2 = e[1], w = e[2];
      if (uf.find(p1) != uf.find(p2)) {
        uf.merge(p1, p2);
        ans += w;
        edgeCount++;
      }
    }
    return ans;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1584. Min Cost to Connect All Points", "[1584]") {
  Solution s;
  vector<vector<int>> points = {{0, 0}, {2, 2}, {3, 10}, {5, 2}, {7, 0}};
  int output = 20;
  REQUIRE(s.minCostConnectPoints(points) == output);
  points = {{0, 0}, {1, 1}, {1, 0}, {-1, 1}};
  output = 4;
  REQUIRE(s.minCostConnectPoints(points) == output);
  points = {{-1000000, -1000000}, {1000000, 1000000}};
  output = 4000000;
  REQUIRE(s.minCostConnectPoints(points) == output);
  points = {{0, 0}};
  output = 0;
  REQUIRE(s.minCostConnectPoints(points) == output);
}