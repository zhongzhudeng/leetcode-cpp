#include "tarjan.hpp"

#include <algorithm>
#include <stack>
#include <vector>

#include "catch2/catch.hpp"

TarjanSCC::TarjanSCC(int n_, const std::vector<std::vector<int>>& edges_,
                     const std::vector<std::vector<int>>& edgesId_)
    : edges(edges_),
      edgesId(edgesId_),
      low(n_, -1),
      dfn(n_, -1),
      n(n_),
      ts(-1) {}

std::vector<int> TarjanSCC::getCuttingEdge() {
  for (int i = 0; i < n; ++i) {
    if (dfn[i] == -1) {
      getCuttingEdge_(i, -1);
    }
  }
  return ans;
}

void TarjanSCC::getCuttingEdge_(int u, int parentEdgeId) {
  low[u] = dfn[u] = ++ts;
  for (int i = 0; i < edges[u].size(); ++i) {
    int v = edges[u][i];
    int id = edgesId[u][i];
    if (dfn[v] == -1) {
      getCuttingEdge_(v, id);
      low[u] = std::min(low[u], low[v]);
      if (low[v] > dfn[u]) {
        ans.push_back(id);
      }
    } else if (id != parentEdgeId) {
      low[u] = std::min(low[u], dfn[v]);
    }
  }
}

std::vector<int> TarjanSCC::getCuttingVertex(int u) {
  std::vector<int> ans;
  std::stack<int> stack;
  std::vector<int> visited(n);
  return ans;
}

TEST_CASE("TarjanSCC") {
  int n = 12;
  std::vector<std::vector<int>> edges{
      {1},    {0, 2, 9},    {1, 3}, {2, 4, 8}, {3, 5, 10}, {4, 6, 7},
      {5, 7}, {6, 8, 9, 5}, {3, 7}, {1, 7},    {4, 11},    {10}};
  std::vector<std::vector<int>> edgesID{
      {0},    {0, 1, 10},    {1, 2}, {2, 3, 8}, {3, 4, 11}, {4, 5, 13},
      {5, 6}, {6, 7, 9, 13}, {8, 7}, {10, 9},   {11, 12},   {12}};

  TarjanSCC t(n, edges, edgesID);

  SECTION("getCuttingEdge") {
    std::vector<int> bridges = t.getCuttingEdge();
    std::sort(bridges.begin(), bridges.end());
    std::vector<int> ans = {0, 11, 12};
    REQUIRE(bridges == ans);
  }

  /*
  SECTION("getCuttingVertex") {
    std::vector<int> vertex = t.getCuttingVertex(0);
    std::vector<int> ans = {1, 4, 10};
    REQUIRE(vertex == ans);
  }
  */
}