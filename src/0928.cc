#include "union_find.hpp"
#include <vector>
using namespace std;

class SolutionII {
public:
  int minMalwareSpread(vector<vector<int>> &graph, vector<int> &initial) {
    int n = graph.size();
    vector<int> initialSet(n);
    for (int v : initial)
      initialSet[v] = 1;

    UnionFind uf(n);
    for (int u = 0; u < n; u++) {
      if (initialSet[u] == 1)
        continue;
      for (int v = 0; v < n; v++) {
        if (initialSet[v] == 1)
          continue;
        if (graph[u][v] == 1)
          uf.merge(u, v);
      }
    }

    vector<vector<int>> infectedBy(n);
    for (int v : initial) {
      vector<int> infectedSet(n);
      for (int u = 0; u < n; u++) {
        if (initialSet[u] == 1 || graph[u][v] == 0)
          continue;
        infectedSet[uf.find(u)] = 1;
      }
      for (int u = 0; u < n; u++)
        if (infectedSet[u] == 1)
          infectedBy[u].push_back(v);
    }

    vector<int> count(n);
    for (int u = 0; u < n; u++) {
      if (infectedBy[u].size() != 1)
        continue;
      int v = infectedBy[u][0];
      for (int w = 0; w < n; w++)
        if (uf.find(w) == uf.find(u))
          count[v]++;
    }

    int res = initial[0];
    for (int v : initial)
      if (count[v] > count[res] || count[v] == count[res] && v < res)
        res = v;

    return res;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("928. Minimize Malware Spread II", "[0928]") {
  SolutionII s;

  vector<vector<int>> graph = {{1, 1, 0}, {1, 1, 0}, {0, 0, 1}};
  vector<int> initial = {0, 1};
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
