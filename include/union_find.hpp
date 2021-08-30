#include <vector>

class UnionFind {
 public:
  std::vector<int> parent;
  std::vector<int> rank;
  int n;
  // 当前连通分量数目
  int setCount;

 public:
  UnionFind(int _n);
  int find(int x);
  bool unite(int x, int y);
  bool connected(int x, int y);
};
