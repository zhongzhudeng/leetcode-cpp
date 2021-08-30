#include "union_find.hpp"

#include <numeric>

UnionFind::UnionFind(int _n) : n(_n), setCount(_n), parent(_n), rank(_n, 1) {
  std::iota(parent.begin(), parent.end(), 0);
}

int UnionFind::find(int x) {
  return parent[x] == x ? x : (parent[x] = find(parent[x]));
}

bool UnionFind::unite(int x, int y) {
  x = find(x);
  y = find(y);
  if (x == y) {
    return false;
  }
  if (rank[x] < rank[y]) {
    std::swap(x, y);
  }
  parent[y] = x;
  rank[x] += rank[y];
  --setCount;
  return true;
}

bool UnionFind::connected(int x, int y) {
  x = find(x);
  y = find(y);
  return x == y;
}
