#pragma once

#include <stack>
#include <vector>

class Tarjan {
public:
  // Tarjan() = default;
  // virtual ~Tarjan() = default;
  virtual std::vector<int> getCuttingEdge() = 0;
  // virtual std::vector<int> getCuttingVertex();
};

class TarjanRecursive : public Tarjan {
private:
  const std::vector<std::vector<int>> &edges;
  const std::vector<std::vector<int>> &edgesId;
  std::vector<int> low;
  std::vector<int> dfn;
  std::vector<int> ans;
  int n;
  int ts;

private:
  void getCuttingEdge_(int u, int parentEdgeId);

public:
  explicit TarjanRecursive(int n_, const std::vector<std::vector<int>> &edges_,
                           const std::vector<std::vector<int>> &edgesId_);
  std::vector<int> getCuttingEdge() override;
  std::vector<int> getCuttingVertex(int u);
};

struct Edge {
  int node;
  int id;
};

class TarjanNonRecursive : public Tarjan {
private:
  const std::vector<std::vector<Edge>> &edges;
  // low, dfn
  std::vector<std::tuple<int, int>> nodes;
  std::vector<int> cutEdge;
  std::vector<int> cutNode;
  // 当前节点,父边,连接点 it
  std::stack<std::tuple<int, int, std::vector<Edge>::const_iterator>> stack;
  int n;
  int ts;

  void run();
  void run_(int u);

public:
  explicit TarjanNonRecursive(int n_,
                              const std::vector<std::vector<Edge>> &edges_);
  std::vector<int> getCuttingEdge() override;
  std::vector<int> getCuttingNode();
};