#include "tarjan.hpp"

#include <algorithm>
#include <stack>
#include <tuple>
#include <vector>

#include "catch2/catch.hpp"

using namespace std;

TarjanRecursive::TarjanRecursive(int n_, const vector<vector<int>>& edges_,
                                 const vector<vector<int>>& edgesId_)
    : edges(edges_),
      edgesId(edgesId_),
      low(n_, -1),
      dfn(n_, -1),
      n(n_),
      ts(-1) {}

vector<int> TarjanRecursive::getCuttingEdge() {
  for (int i = 0; i < n; ++i) {
    if (dfn[i] == -1) {
      getCuttingEdge_(i, -1);
    }
  }
  return ans;
}

void TarjanRecursive::getCuttingEdge_(int u, int parentEdgeId) {
  low[u] = dfn[u] = ++ts;
  for (int i = 0; i < edges[u].size(); ++i) {
    int v = edges[u][i];
    int id = edgesId[u][i];
    if (dfn[v] == -1) {
      getCuttingEdge_(v, id);
      low[u] = min(low[u], low[v]);
      if (low[v] > dfn[u]) {
        ans.push_back(id);
      }
    } else if (id != parentEdgeId) {
      low[u] = min(low[u], dfn[v]);
    }
  }
}

vector<int> TarjanRecursive::getCuttingVertex(int u) {
  vector<int> ans;
  stack<int> stack;
  vector<int> visited(n);
  return ans;
}

TEST_CASE("TarjanSCC") {
  int n = 12;
  vector<vector<int>> edges{{1},        {0, 2, 9}, {1, 3},  {2, 4, 8},
                            {3, 5, 10}, {4, 6, 7}, {5, 7},  {6, 8, 9, 5},
                            {3, 7},     {1, 7},    {4, 11}, {10}};
  vector<vector<int>> edgesID{{0},        {0, 1, 10}, {1, 2},   {2, 3, 8},
                              {3, 4, 11}, {4, 5, 13}, {5, 6},   {6, 7, 9, 13},
                              {8, 7},     {10, 9},    {11, 12}, {12}};

  TarjanRecursive t(n, edges, edgesID);

  SECTION("getCuttingEdge") {
    vector<int> bridges = t.getCuttingEdge();
    sort(bridges.begin(), bridges.end());
    vector<int> ans = {0, 11, 12};
    REQUIRE(bridges == ans);
  }

  /*
  SECTION("getCuttingVertex") {
    vector<int> vertex = t.getCuttingVertex(0);
    vector<int> ans = {1, 4, 10};
    REQUIRE(vertex == ans);
  }
  */
}

TarjanNonRecursive::TarjanNonRecursive(int n_,
                                       const vector<vector<Edge>>& edges_)
    : n(n_), ts(-1), edges(edges_), nodes(n_, make_tuple(-1, -1)) {}

vector<int> TarjanNonRecursive::getCuttingEdge() {
  int low;
  for (size_t i = 0; i < n; i++) {
    tie(low, ignore) = nodes[i];
    if (low == -1) {
      getCuttingEdge_(i);
    }
  }
  sort(ans.begin(), ans.end());
  return ans;
}

void TarjanNonRecursive::getCuttingEdge_(int u) {
  vector<vector<Edge>> init = {{{-1, 0}}};
  stack.push(make_tuple(-1, -1, init[0].begin()));
  int cn = u, pe = -1, low = -1, dfn = -1, nlow = -1, ndfn = -1;
  low = dfn = ++ts;
  nodes[cn] = make_tuple(low, dfn);
  auto ne = edges[cn].begin();
  stack.push(make_tuple(cn, pe, ne));

  while (!stack.empty()) {
    if (ne == edges[cn].end()) {
      stack.pop();
      tie(nlow, ndfn) = nodes[cn];
      tie(cn, pe, ne) = stack.top();
      if (cn == -1) {
        stack.pop();
        continue;
      }
      tie(low, dfn) = nodes[cn];
      low = min(low, nlow);
      nodes[cn] = make_tuple(low, dfn);
      if (dfn < nlow) {
        ans.push_back(ne->id);
      }
      ne++;
      stack.pop();
      stack.push(make_tuple(cn, pe, ne));
      continue;
    }

    if (ne->id == pe) {
      ne++;
      stack.pop();
      stack.push(make_tuple(cn, pe, ne));
      continue;
    }
    tie(nlow, ndfn) = nodes[ne->node];
    tie(low, dfn) = nodes[cn];
    if (ndfn == -1) {
      pe = ne->id;
      cn = ne->node;
      nlow = ndfn = ++ts;
      nodes[cn] = make_tuple(nlow, ndfn);
      ne = edges[cn].begin();
      stack.push(make_tuple(cn, pe, ne));
    } else {
      low = min(low, nlow);
      nodes[cn] = make_tuple(low, dfn);
      if (dfn < nlow) {
        ans.push_back(ne->id);
      }
      ne++;
      stack.pop();
      stack.push(make_tuple(cn, pe, ne));
    }
  }
}

TEST_CASE("TJ") {
  SECTION("getCuttingEdge") {
    vector<vector<Edge>> edges{{{1, 0}, {12, 13}},                 // 0
                               {{0, 0}, {2, 1}, {9, 10}},          // 1
                               {{1, 1}, {3, 2}},                   // 2
                               {{2, 2}, {4, 3}, {8, 8}},           // 3
                               {{3, 3}, {5, 4}, {10, 11}},         // 4
                               {{4, 4}, {6, 5}, {7, 13}},          // 5
                               {{5, 5}, {7, 6}},                   // 6
                               {{6, 6}, {8, 7}, {9, 9}, {5, 13}},  // 7
                               {{3, 8}, {7, 7}},                   // 8
                               {{1, 10}, {7, 9}},                  // 9
                               {{4, 11}, {11, 12}},                // 10
                               {{10, 12}},                         // 11
                               {{0, 13}}};

    TarjanNonRecursive t(edges.size(), edges);
    vector<int> bridges = t.getCuttingEdge();
    sort(bridges.begin(), bridges.end());
    vector<int> ans = {0, 11, 12, 13};
    REQUIRE(bridges == ans);
  }
  SECTION("getCuttingEdge") {
    vector<vector<Edge>> edges = {
        {{1, 0}},                           // 0
        {{0, 0}, {2, 1}, {9, 10}},          // 1
        {{1, 1}, {3, 2}},                   // 2
        {{2, 2}, {4, 3}, {8, 8}},           // 3
        {{3, 3}, {5, 4}, {10, 11}},         // 4
        {{4, 4}, {6, 5}, {7, 13}},          // 5
        {{5, 5}, {7, 6}},                   // 6
        {{6, 6}, {8, 7}, {9, 9}, {5, 13}},  // 7
        {{3, 8}, {7, 7}},                   // 8
        {{1, 10}, {7, 9}},                  // 9
        {{4, 11}, {11, 12}},                // 10
        {{10, 12}}                          // 11
    };
    TarjanNonRecursive t(edges.size(), edges);
    vector<int> bridges = t.getCuttingEdge();
    sort(bridges.begin(), bridges.end());
    vector<int> ans = {0, 11, 12};
    REQUIRE(bridges == ans);
  }
}
