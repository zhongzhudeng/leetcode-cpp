#include "tarjan.hpp"

#include <vector>

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