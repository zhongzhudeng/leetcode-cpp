#include <stack>
#include <unordered_map>
#include <vector>
using namespace std;

class SCC {
private:
  vector<vector<int>> &adj; // the graph
  vector<int> &initial;

  stack<int> S; // stack
  int dfs_counter = 0;
  const int UNVISITED = -1;

  vector<int> dfs_num;
  vector<int> dfs_min;
  vector<bool> on_stack;
  // node -> scc_id
  vector<int> ids;
  // scc_id -> scc size
  unordered_map<int, int> id_size;
  int id;

  void dfs(int u) {
    dfs_min[u] = dfs_num[u] = dfs_counter++;
    S.push(u);
    on_stack[u] = true;
    for (int v = 0; v < adj.size(); v++) {
      if (adj[u][v] == 1) {
        if (dfs_num[v] == UNVISITED) {
          dfs(v);
          dfs_min[u] = min(dfs_min[u], dfs_min[v]);
        } else if (on_stack[v]) // only on_stack can use back edge
          dfs_min[u] = min(dfs_min[u], dfs_num[v]);
      }
    }
    if (dfs_min[u] == dfs_num[u]) { // output result
      vector<int> scc_;
      int v = -1;
      while (v != u) { // output SCC starting in u
        v = S.top();
        S.pop();
        on_stack[v] = false;
        scc_.push_back(v);
      }
      id++;
      for (auto i : scc_)
        ids[i] = id;
      id_size[id] = scc_.size();
    }
  }

public:
  SCC(vector<vector<int>> &adj, vector<int> &initial)
      : adj(adj), initial(initial), dfs_num(adj.size(), UNVISITED),
        dfs_min(adj.size(), UNVISITED), on_stack(adj.size(), false),
        ids(adj.size(), UNVISITED) {
    id = 0;
  }
  int scc() {
    for (int i = 0; i < adj.size(); i++) {
      if (dfs_num[i] == UNVISITED)
        dfs(i);
    }
    unordered_map<int, int> id_init;
    for (auto i : initial)
      id_init[ids[i]]++;

    int ans = adj.size() + 1, ans_removed = 0;
    for (int u : initial) {
      int removed = (id_init[ids[u]] == 1 ? id_size[ids[u]] : 0);
      if (removed > ans_removed || (removed == ans_removed && u < ans)) {
        ans = u;
        ans_removed = removed;
      }
    }
    return ans;
  }
};

class SolutionI {
public:
  int minMalwareSpread(vector<vector<int>> &graph, vector<int> &initial) {
    SCC s(graph, initial);
    return s.scc();
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("924. Minimize Malware Spread", "[0924]") {
  SolutionI s;

  std::vector<std::vector<int>> graph = {{1, 1, 0}, {1, 1, 0}, {0, 0, 1}};
  std::vector<int> initial = {0, 1};
  int ans = 0;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  initial = {0, 2};
  ans = 0;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  initial = {1, 2};
  ans = 1;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);

  graph = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 1, 1}};
  initial = {3, 1};
  ans = 3;
  REQUIRE(s.minMalwareSpread(graph, initial) == ans);
}
