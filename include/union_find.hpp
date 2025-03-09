#pragma once

#include <vector>

class UnionFind {
  std::vector<int> fa;
  std::vector<int> rank;

public:
  UnionFind(int n);
  int find(int i);
  void merge(int u, int v);
};
