#pragma once
#include <vector>

class TarjanSCC {
 private:
  const std::vector<std::vector<int>>& edges;
  const std::vector<std::vector<int>>& edgesId;
  std::vector<int> low;
  std::vector<int> dfn;
  std::vector<int> ans;
  int n;
  int ts;

 private:
  void getCuttingEdge_(int u, int parentEdgeId);

 public:
  TarjanSCC(int n_, const std::vector<std::vector<int>>& edges_,
            const std::vector<std::vector<int>>& edgesId_);
  std::vector<int> getCuttingEdge();
};