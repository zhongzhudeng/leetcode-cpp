#include "union_find.hpp"
#include <numeric>

UnionFind::UnionFind(int n) : fa(n), rank(n, 1) {
  std::iota(fa.begin(), fa.end(), 0);
}

int UnionFind::find(int i) { return i == fa[i] ? i : (fa[i] = find(fa[i])); }

void UnionFind::merge(int u, int v) {
  int x = find(u), y = find(v);
  if (rank[x] <= rank[y])
    fa[x] = y;
  else
    fa[y] = x;
  if (rank[x] == rank[y] && x != y)
    rank[y]++;
}
