#include <algorithm>
#include <tarjan.hpp>
#include <vector>

#include "catch2/catch.hpp"

using namespace std;

class Solution {
 public:
  vector<vector<int>> criticalConnections(int n,
                                          vector<vector<int>>& connections);
};

vector<vector<int>> Solution::criticalConnections(
    int n, vector<vector<int>>& connections) {
  std::vector<std::vector<Edge>> edges(n);
  int i = 0;
  for (auto& e : connections) {
    auto first = e[0], second = e[1];
    Edge edge = {second, i};
    edges[first].push_back(edge);
    edge.node = first;
    edges[second].push_back(edge);
    i++;
  }
  TarjanNonRecursive tj(n, edges);
  auto res = tj.getCuttingEdge();
  vector<vector<int>> ans;
  for (auto& e : res) {
    ans.push_back(connections[e]);
  }
  for (auto& e : ans) {
    sort(e.begin(), e.end());
  }
  return ans;
}

TEST_CASE("criticalConnections") {
  Solution s;
  int n = 4;
  vector<vector<int>> in = {{0, 1}, {1, 2}, {2, 0}, {1, 3}};
  vector<vector<int>> ans = {{1, 3}};
  REQUIRE(s.criticalConnections(n, in) == ans);

  n = 5;
  in = {{1, 0}, {2, 0}, {3, 2}, {4, 2}, {4, 3}, {3, 0}, {4, 0}};
  ans = {{0, 1}};
  REQUIRE(s.criticalConnections(n, in) == ans);
}