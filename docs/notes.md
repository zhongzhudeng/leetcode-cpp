# Competitive Programming

## Files

```cmake
cmake_minimum_required(VERSION 3.0)
project(main VERSION 0.1.0)

add_executable(main main.cpp)
target_compile_options(main
  PUBLIC -Wall
  PUBLIC -std=c++17
)

{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "build",
      "type": "shell",
      "command": "cmake --build ${workspaceFolder}/build --config Debug --target all --",
      "problemMatcher": [],
      "group": {
        "kind": "build",
        "isDefault": true
      }
    }
  ]
}

{
  "version": "0.2.0",
  "configurations": [
    {
      "type": "lldb",
      "request": "launch",
      "name": "Debug",
      "program": "${workspaceFolder}/build/main",
      "args": [],
      "cwd": "${workspaceFolder}",
      "stdio": ["samples-marsmaps/1.in", null, null],
      "preLaunchTask": "build",
    }
  ]
}
```

# Final

```cpp
// 1
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

class segtree {
 private:
  vector<int> values;
  size_t size;
  int parent(int i) { return i / 2; }
  int left(int i) { return 2 * i; }
  int right(int i) { return 2 * i + 1; }
  void update(int i) {
    values[i] = min(values[left(i)], values[right(i)]);
    if (i > 1) update(parent(i));
  }
  bool query(int i, int j, int l, int r, int current_node) {
    if (r <= i || j <= l) return 1;
    if (r <= j && i <= l) return values[current_node];
    int m = (l + r) / 2;
    return min(query(i, j, l, m, left(current_node)),
               query(i, j, m, r, right(current_node)));
  }

 public:
  segtree(size_t n) {
    size = 1 << (int)ceil(log2(n));
    values.assign(2 * size, 0);
  }
  int64_t query(int i, int j) { return query(i, j, 0, size, 1); }
  void update(int i, int val) {
    values[i + size] += val;
    update(parent(i + size));
  }
};

int main(int argc, const char** argv) {
  int n, add, val, s, cnt = 0;
  string op;
  scanf("%d\n", &n);
  segtree seg(1 << 16);
  for (int i = 0; i < n; i++) {
    cin >> op;
    if (op == "malloc") {
      cin >> add;
      seg.update(add, 1);
    } else if (op == "free") {
      cin >> add;
      seg.update(add, -1);
    } else {
      cin >> add >> val >> s;
      if (seg.query(add, add + s) == 0) {
        printf("%d\n", i + 1);
        cnt++;
      }
    }
  }
  if (cnt==0) printf("OK\n");
  return 0;
}

// 2
#include <cstdio>
#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>
using namespace std;

int main(int argc, const char** argv) {
  int n;
  int64_t s, e;
  scanf("%d\n", &n);
  vector<pair<int64_t, int64_t>> conf;
  priority_queue<int64_t, vector<int64_t>, greater<int64_t>> pq;

  for (int i = 0; i < n; i++) {
    cin >> s >> e;
    conf.emplace_back(s, e);
  }
  sort(conf.begin(), conf.end());
  pq.emplace(conf[0].second);
  int cnt = 1;
  for (int i = 1; i < n; i++) {
    s = conf[i].first;
    e = conf[i].second;
    int64_t fe = pq.top();
    if (s < fe)
      cnt++;
    else
      pq.pop();
    pq.emplace(e);
  }
  printf("%d\n", cnt);
}

// 3
#include <cstdio>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#define node_req 16
#define num_type 4

using namespace std;

class Dinic {
 private:
  static constexpr int inf = numeric_limits<int>::max() / 2;
  int n, s, t;
  vector<vector<int>> adj;             // adjacency list
  vector<vector<int>> capacity, flow;  // adjacency matrix
  vector<int> d, cur;
  vector<bool> vis;

  bool bfs() {
    vis.assign(n, false);
    queue<int> Q;
    Q.push(s);
    d[s] = 0;
    vis[s] = true;
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (int v : adj[u]) {
        if (!vis[v] && capacity[u][v] > flow[u][v]) {
          vis[v] = true;
          d[v] = d[u] + 1;
          Q.push(v);
        }
      }
    }
    return vis[t];
  }
  int dfs(int u, int a) {
    if (u == t || a == 0) return a;
    int fl = 0, f;
    for (auto v : adj[u]) {
      if (d[u] + 1 == d[v] &&
          (f = dfs(v, min(a, capacity[u][v] - flow[u][v]))) > 0) {
        flow[u][v] += f;
        flow[v][u] -= f;
        fl += f;
        a -= f;
        if (a == 0) break;
      }
    }
    return fl;
  }

 public:
  Dinic(int n) {
    this->n = n;
    adj.resize(n);
    capacity.assign(n, vector<int>(n, 0)), flow.assign(n, vector<int>(n, 0));
    d.resize(n), cur.resize(n), vis.resize(n);
  }
  void addEdge(int u, int v, int c) {
    adj[u].push_back(v);
    adj[v].push_back(u);
    capacity[u][v] += c;
    // capacity[v][u] += c;
  }
  int maxFlow(int s, int t) {
    this->s = s;
    this->t = t;
    int flow = 0;
    while (bfs()) {
      cur.assign(n, 0);
      flow += dfs(s, inf);
    }
    return flow;
  }
};

int maxFlow(vector<int> &request, vector<int> &idxtocnt) {
  Dinic fl(node_req + num_type + 1);  // 21
  for (int i = 1; i < node_req; i++) {
    fl.addEdge(i, 20, request[i]);
    for (int j = 0; j < num_type; j++)
      if (i & 1 << j) fl.addEdge(j + node_req, i, idxtocnt[j]);
  }
  for (int i = 0; i < num_type; i++) fl.addEdge(0, node_req + i, idxtocnt[i]);
  return fl.maxFlow(0, 20);
}

int main(int argc, const char **argv) {
  int h, di;
  map<char, int> chtoidx{{'l', 0}, {'o', 1}, {'e', 2}, {'g', 3}};
  scanf("%d", &h);
  vector<int> idxtocnt(4);
  for (int i = 0; i < num_type; i++) scanf("%d", &idxtocnt[i]);

  vector<int> request(node_req, 0);
  for (int i = 0; i < h; i++) {
    scanf("%d", &di);
    int idx = 0;
    for (int j = 0; j < di; j++) {
      char ch;
      cin >> ch;
      idx |= 1 << chtoidx[ch];
    }
    request[idx] += 1;
  }

  int ans1, ans2, ans3;
  ans1 = maxFlow(request, idxtocnt);

  idxtocnt[3] = numeric_limits<int>::max() / 2;
  ans2 = maxFlow(request, idxtocnt);

  int l = -1, r = numeric_limits<int>::max() / 2;
  while (l < r - 1) {
    int m = (l + r) / 2;
    idxtocnt[3] = m;
    int mf = maxFlow(request, idxtocnt);
    if (mf < ans2)
      l = m;
    else
      r = m;
  }
  ans3 = r;

  printf("%d %d %d\n", ans1, ans2, ans3);

  return 0;
}

// 4
#include <cstdio>
#include <iostream>
#include <queue>
#include <set>
#include <vector>
using namespace std;

class Tarjan {
  static constexpr int inf = numeric_limits<int>::max() / 2;
  const int UNVISITED = -1;
  int dfs_counter;
  vector<int> dfs_num, dfs_min, dfs_parent;
  vector<vector<int>> &adj_;
  set<pair<int, int>> bridges;
  int n;

  void dfs(int u) {
    dfs_min[u] = dfs_num[u] = dfs_counter++;
    for (auto v : adj_[u]) {
      if (dfs_num[v] == UNVISITED) {  // Tree Edge
        dfs_parent[v] = u;
        dfs(v);
        if (dfs_num[u] < dfs_min[v]) bridges.emplace(u, v);
        dfs_min[u] = min(dfs_min[u], dfs_min[v]);
      } else if (v != dfs_parent[u])  // Back Edge
        dfs_min[u] = min(dfs_min[u], dfs_num[v]);
    }
  }

 public:
  Tarjan(vector<vector<int>> &adj) : adj_(adj) { n = adj_.size(); }
  set<pair<int, int>> articulation_points_and_bridges() {
    dfs_num.assign(n, UNVISITED), dfs_min.assign(n, UNVISITED),
        dfs_parent.assign(n, -1), bridges.clear();
    dfs_counter = 0;

    for (int i = 0; i < n; i++)
      if (dfs_num[i] == UNVISITED) dfs(i);

    return bridges;
  }
};

int main(int argc, const char **argv) {
  int n, e, m, idx, u, v;
  scanf("%d %d %d\n", &n, &e, &m);
  vector<bool> entr(n + 1, false);
  vector<vector<int>> adj_in(n + 1, vector<int>(0));
  for (int i = 0; i < e; i++) {
    scanf("%d", &idx);
    entr[idx] = true;
    adj_in[0].push_back(idx), adj_in[idx].push_back(0);
  }

  for (int i = 0; i < m; i++) {
    scanf("%d %d\n", &u, &v);
    adj_in[u].push_back(v), adj_in[v].push_back(u);
  }

  Tarjan tj(adj_in);
  auto b = tj.articulation_points_and_bridges();

  vector<vector<int>> adj(n + 1, vector<int>(0));
  for (int u = 1; u < n + 1; u++) {
    for (auto v : adj_in[u]) {
      if (v != 0 && b.find(make_pair(u, v)) == b.end() &&
          b.find(make_pair(v, u)) == b.end()) {
        adj[u].push_back(v);
      }
    }
  }

  for (int i = 1; i < n + 1; i++)
    if (entr[i]) adj[0].push_back(i), adj[i].push_back(0);

  vector<bool> visited(n + 1, false);
  queue<int> Q;       // FIFO
  Q.push(0);          //
  visited[0] = true;  // start vertex
  while (!Q.empty()) {
    int v = Q.front();
    Q.pop();
    for (int u : adj[v])
      if (!visited[u]) {
        Q.push(u);
        visited[u] = true;
      }
  }
  int cnt = 0;
  for (int i = 1; i < n + 1; i++)
    if (visited[i]) cnt++;

  printf("%d\n", cnt);

  return 0;
}

// 5
#include <cstdint>
#include <cstdio>
#include <vector>
using namespace std;
int64_t fexp(int64_t m, int64_t n, int64_t p) {
  if (n == 0)
    return 1;
  else if (n % 2 == 1)
    return (m * fexp(m, n - 1, p)) % p;
  else {  // n is even
    int64_t r = fexp(m, n / 2, p);
    return (r * r) % p;
  }
}

int64_t comb(int64_t n, int64_t k, int64_t p) {
  int64_t ni = 1, ki = 1;
  for (int64_t i = n; i > n - k; i--) {
    ni = (ni * i) % p;
  }

  for (int64_t i = k; i > 0; i--) {
    ki = (ki * i) % p;
  }
  int kni = fexp(ki, p - 2, p);

  return (ni * kni) % p;
}
int main(int argc, const char** argv) {
  int n, k;
  scanf("%d %d\n", &n, &k);
  vector<int64_t> dp(k + 1, 1);
  const int64_t p = 1'000'000'007;

  for (int i = 2; i < k + 1; i++) {
    dp[i] = fexp(i, n, p);
    for (int j = 1; j < i; j++) dp[i] = (dp[i] - comb(i, j, p) * dp[j] % p) % p;
  }
  printf("%lld\n", dp[k] < 0 ? p + dp[k] : dp[k]);

  return 0;
}

```

## UnionFind

```cpp
struct UnionFind {
  vector<int> fa;
  vector<int> rank;

  UnionFind(int n) {
    fa.resize(n);
    iota(fa.begin(), fa.end(), 0);
    rank.assign(n, 1);
  }
  int find(int i) { return i == fa[i] ? i : fa[i] = find(fa[i]); }
  void merge(int u, int v) {
    int x = find(u), y = find(v);
    if (rank[x] <= rank[y])
      fa[x] = y;
    else
      fa[y] = x;
    if (rank[x] == rank[y] && x != y) rank[y]++;
  }
};
```

## SegmentTree

```cpp
class segtree {
 private:
  vector<int64_t> values;
  size_t size;
  int parent(int i) { return i / 2; }
  int left(int i) { return 2 * i; }
  int right(int i) { return 2 * i + 1; }
  void update(int i) {
    values[i] = values[left(i)] + values[right(i)];
    if (i > 1) update(parent(i));
  }
  int64_t query(int i, int j, int l, int r, int current_node) {
    if (r <= i || j <= l) return 0;
    if (r <= j && i <= l) return values[current_node];
    int m = (l + r) / 2;
    return query(i, j, l, m, left(current_node)) +
           query(i, j, m, r, right(current_node));
  }

 public:
  segtree(size_t n) {
    size = 1 << (int)ceil(log2(n));
    values.assign(2 * size, 0);
  }
  segtree(vector<int64_t> v) : segtree(v.size()) {
    for (size_t i = 0; i < v.size(); ++i) values[i + size] = v[i];
    for (size_t i = size - 1; i > 0; --i)
      values[i] = values[left(i)] + values[right(i)];
  }
  int64_t query(int i, int j) { return query(i, j, 0, size, 1); }
  void update(int i, int64_t val) {
    values[i + size] = val;
    update(parent(i + size));
  }
};
```

### SegmentTree range update

```cpp
class segtree {
 private:
  vector<int64_t> value;
  vector<int> count;
  vector<int>& X;
  size_t size;

  int parent(int i) { return i / 2; }
  int left(int i) { return 2 * i; }
  int right(int i) { return 2 * i + 1; }
  bool is_leaf(int i) { return i >= size; }

  void add_impl(size_t i, size_t j, int val, size_t node, size_t l, size_t r) {
    if (i <= l && r <= j) {
      count[node] += val;
      if (count[node])
        value[node] = X[r] - X[l];
      else if (is_leaf(node))
        value[node] = 0;
      else
        value[node] = value[left(node)] + value[right(node)];
    } else if (i < r && l < j) {
      auto m = (l + r) / 2;
      add_impl(i, j, val, left(node), l, m);
      add_impl(i, j, val, right(node), m, r);
      if (count[node] == 0)
        value[node] = value[left(node)] + value[right(node)];
    }
  }

 public:
  segtree(size_t n, vector<int>& x) : X(x) {
    size = 1 << (int)ceil(log2(n));
    count.assign(2 * size, 0);
    value.assign(2 * size, 0);
  }
  int64_t add(size_t i, size_t j, int val) {
    add_impl(i, j, val, 1, 0, size);
    return value[1];
  }
};

int main(int argc, const char** argv) {
  int n;
  scanf("%d\n", &n);
  vector<tuple<int, int, int, int>> events;  // x, y1, y2, value
  set<int> Xvals;
  for (int i = 0; i < n; i++) {
    int x1, y1, x2, y2;
    scanf("%d %d %d %d\n", &x1, &y1, &x2, &y2);
    events.emplace_back(y1, x1, x2, 1);
    events.emplace_back(y2, x1, x2, -1);
    Xvals.emplace(x1), Xvals.emplace(x2);
  }
  sort(events.begin(), events.end());
  vector<int> X(Xvals.begin(), Xvals.end());  // i->x
  map<int, int> Xi;                           // x->i
  for (int i = 0; i < X.size(); i++) Xi[X[i]] = i;

  segtree s(Xvals.size() - 1, X);

  int64_t ans = 0;
  int64_t cur_x_sum = 0;
  int cur_y = get<0>(events[0]);
  for (auto [y, x1, x2, v] : events) {
    ans += cur_x_sum * (y - cur_y);
    cur_x_sum = s.add(Xi[x1], Xi[x2], v);
    cur_y = y;
  }
  cout << ans << "\n";
}
```

## Binary Search

```cpp
int l = -1, r = numeric_limits<int>::max() / 2;
while (l < r - 1) {
  int m = (l + r) / 2;
  if (!check(m))
    l = m;
  else
    r = m;
}
return r
```

## Dynamic Programming

```cpp
int knapsack(int n, int maxW, vector<int> &w, vector<int> &p,
             vector<int> &result) {
  vector<vector<int>> dp(n, vector<int>(maxW + 1));
  for (int v = 0; v <= maxW; ++v) dp[0][v] = (v >= w[0] ? p[0] : 0);
  for (int i = 1; i < n; ++i) {
    for (int v = 0; v <= maxW; ++v) {
      dp[i][v] = dp[i - 1][v];
      if (v >= w[i]) dp[i][v] = max(dp[i][v], p[i] + dp[i - 1][v - w[i]]);
    }
  }
  int remainingW = maxW;
  for (int i = n - 1; i > 0; --i) {
    if (dp[i - 1][remainingW] == dp[i][remainingW]) {
    } else {
      result.push_back(i);
      remainingW -= w[i];
    }
  }
  if (remainingW >= w[0]) {
    result.push_back(0);
  }
  return dp[n - 1][maxW];
}

int lis(vector<int>& in) {
  vector<int> st;
  for (auto& x : in) {
    auto it = lower_bound(st.begin(), st.end(), x);
    if (it == st.end())
      st.push_back(x);
    else
      *it = x;
  }
  return st.size();
}
```

## Graph

### BFS/DFS

```cpp
void bfs(int start) {
  vector<bool> visited(V, false);
  queue<int> Q;           // FIFO
  stack<int> S;
  Q.push(start);          //
  visited[start] = true;  // start vertex
  while (!Q.empty()) {
    int v = Q.front();
    Q.pop();
    for (int u : adj[v])
      if (!visited[u]) {
        Q.push(u);
        visited[u] = true;
      }
  }
}

vector<bool> visited (V, false );
void dfs(int v) {
  visited[v] = true;
  for (int u : adj[v])
    if (!visited[u]) dfs(u);
}

// count SCC
int c = 0;
vector<bool> visited(V, false);
for (int i = 0; i < V; i++) {
  if (!visited[i]) {
    dfs(i);
    c++;
  }
}

// Toposort
int V; // number of nodes
vector<bool> visited(V, false);
deque<int> ts;     // the final topological sort
void dfs(int v) {  // modified dfs for toposort 
visited [v] = true;
  for (int u : adj[v])
    if (!visited[u]) dfs(u);
  ts.push_front(v);
}

// Shortest Path on DAG
vector<int> dist(V, INF);
dist[start] = 0;  // initialize all source vertices
// assume given toposort
vector<int> ts = compute_toposort();
for (auto u : ts)
	for (auto p : adj[u]) {
		int v = p.first;
		int w = p.second;  // relax (u,v)
		dist[v] = min(dist[v], dist[u] + w);
	}

// bipartite
int V;                      // number of nodes
vector<vector<int>> adj;    // the graph
vector<int> colors(V, -1);  // -1 means unvisited

void dfs(int v) {
  for (auto u : adj[v])
    if (colors[u] == -1) {
      colors[u] = 1 - colors[v];
      dfs(u);
    } else if (colors[u] == colors[v]) {
      cout << "Impossible" << endl;
      exit(0);
    }
}

void is_bipartite(int start) {
  colors[start] = 0;  // colors are 0, 1
  // assume adj to be connected
  dfs(start);
  cout << "Possible" << endl;
}

//find cycles
int V;                    // number of nodes
vector<vector<int>> adj;  // the graph

int UNVISITED = 0, EXPLORED = 1, VISITED = 2;
vector<int> visited(V, UNVISITED);

void dfs(int v) {
  visited[v] = EXPLORED;
  for (auto u : adj[v])
    if (visited[u] == UNVISITED) {
      dfs(u);
    } else {  // not part of dfs tree
      if (visited[u] == EXPLORED) {
        cout << "Cycle found" << endl;
        exit(0);
      }
    }
  visited[v] = VISITED;
}

void find_cycles(int start) {
  dfs(start);
  cout << "Graph is acyclic" << endl;
}
```

### Shortest Path

```cpp
// Dijkstra
vector<int> dist(V, INF);
dist[start] = 0;
priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>>
    pq;
pq.push({0, start});  // <distance , vertex>
while (!pq.empty()) {
  auto front = pq.top();
  pq.pop();
  int d = front.first, v = front.second;
  if (d > dist[v]) continue;  // lazy deletion
  for (auto p : adj[v]) {     // <target , weight>
    int u = p.first;
    int w = p.second;
    if (dist[v] + w < dist[u]) {
      dist[u] = dist[v] + w;
      pq.push({dist[u], u});  // can push
    }                         // duplicate vertices
  }
}

// Bellman-Ford/SPFA
vector<int> dist(V, INF);
queue<int> Q;
vector<bool> inQ(V, false);
dist[start] = 0;
Q.push(start);
inQ[start] = true;
while (!Q.empty()) {
	int v = Q.front();
	Q.pop();
	inQ[v] = false;
	for (auto p : adj[v]) {
		int u = p.first;
		int w = p.second;
		if (dist[u] > dist[v] + w) {
			dist[u] = dist[v] + w;
			if (!inQ[u]) {
				Q.push(u);
				inQ[u] = true;
			}
		}
	}
}

// Floyd-Warshall
  vector<vector<int>> adj;  // adjacency matrix !
  for (int k = 0; k < V; k++)
    for (int i = 0; i < V; i++)
      for (int j = 0; j < V; j++)
        adj[i][j] = min(adj[i][j], adj[i][k] + adj[k][j]);

```

### Tarjan

```cpp
// articulation points and bridges
int V;                    // number of nodes
vector<vector<int>> adj;  // the graph

int dfs_counter = 0;
const int UNVISITED = -1;
int dfsRoot, rootChildren;

vector<int> dfs_num(V, UNVISITED);
vector<int> dfs_min(V, UNVISITED);
vector<int> dfs_parent(V, -1);

void dfs(int u) {
  dfs_min[u] = dfs_num[u] = dfs_counter++;
  for (auto v : adj[u]) {
    if (dfs_num[v] == UNVISITED) {  // Tree Edge
      dfs_parent[v] = u;
      if (u == dfsRoot) rootChildren++;

      dfs(v);

      if (dfs_num[u] <= dfs_min[v] && u != dfsRoot)
        cout << u << " is AP" << endl;
      if (dfs_num[u] < dfs_min[v])
        cout << u << "-" << v << " is Bridge" << endl;
      dfs_min[u] = min(dfs_min[u], dfs_min[v]);
    } else if (v != dfs_parent[u])  // Back Edge
      dfs_min[u] = min(dfs_min[u], dfs_num[v]);
  }
}

void articulation_points_and_bridges() {
  for (int i = 0; i < V; i++)
    if (dfs_num[i] == UNVISITED) {
      dfsRoot = i;
      rootChildren = 0;
      dfs(i);  // code on next slide
      if (rootChildren > 1) cout << i << " is AP" << endl;
    }
}

```

```cpp
// scc
int V;                    // number of nodes
vector<vector<int>> adj;  // the graph

stack<int> S;  // stack
int dfs_counter = 0;
const int UNVISITED = -1;

vector<int> dfs_num(V, UNVISITED);
vector<int> dfs_min(V, UNVISITED);
vector<bool> on_stack(V, false);

void dfs(int u) {
  dfs_min[u] = dfs_num[u] = dfs_counter++;
  S.push(u);
  on_stack[u] = true;
  for (auto v : adj[u]) {
    if (dfs_num[v] == UNVISITED) {
      dfs(v);
      dfs_min[u] = min(dfs_min[u], dfs_min[v]);
    } else if (on_stack[v])  // only on_stack can use back edge
      dfs_min[u] = min(dfs_min[u], dfs_num[v]);
  }
  if (dfs_min[u] == dfs_num[u]) {  // output result
    cout << "SCC: ";
    int v = -1;
    while (v != u) {  // output SCC starting in u
      v = S.top();
      S.pop();
      on_stack[v] = false;
      cout << v << " ";
    }
    cout << endl;
  }
}

void scc() {
  for (int i = 0; i < V; i++) {
    if (dfs_num[i] == UNVISITED) dfs(i);  // on next slide
  }
}
```

### MaxFlow

```cpp
class EK {
 private:
  vector<vector<int>> capacity;  // adjacency matrix
  vector<vector<int>> adj;       // adjacency list
  vector<int> parent;            // used in BFS
  static constexpr int inf = numeric_limits<int>::max() / 2;
  void bfs(int s) ,u, v{
    parent.assign(adj.size(), -1);
    parent[s] = -2;  // s is visited
    queue<int> Q;
    Q.push(s);
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (int v : adj[u])
        if (parent[v] == -1 and capacity[u][v] > 0) {
          Q.push(v);
          parent[v] = u;
        }
    }
  }

 public:
  EK(int n_) {
    adj.resize(n_);
    capacity.assign(n_, vector<int>(n_, 0));
  }
  void addEdge(int u, int v, int c) {
    adj[u].push_back(v);
    adj[v].push_back(u);
    capacity[u][v] += c;
    capacity[v][u] += c; // for undirected graph
  }
  int maxFlow(int s, int t) {
    int totalflow = 0, u;
    while (true) {
      bfs(s);                      // build bfs tree
      if (parent[t] == -1) break;  // unreachable
      int bottleneck = inf;
      u = t;  // find bottleneck capacity
      while (u != s) {
        int v = parent[u];
        bottleneck = min(bottleneck, capacity[v][u]);
        u = v;
      }
      u = t;  // update capacities along path
      while (u != s) {
        int v = parent[u];
        capacity[v][u] -= bottleneck;
        capacity[u][v] += bottleneck;
        u = v;
      }
      totalflow += bottleneck;
    }
    return totalflow;
  }
};

class Dinic {
 private:
  static constexpr int inf = numeric_limits<int>::max() / 2;
  int n, s, t;
  vector<vector<int>> adj;             // adjacency list
  vector<vector<int>> capacity, flow;  // adjacency matrix
  vector<int> d, cur;
  vector<bool> vis;

  bool bfs() {
    vis.assign(n, false);
    queue<int> Q;
    Q.push(s);
    d[s] = 0;
    vis[s] = true;
    while (!Q.empty()) {
      int u = Q.front();
      Q.pop();
      for (int v : adj[u]) {
        if (!vis[v] && capacity[u][v] > flow[u][v]) {
          vis[v] = true;
          d[v] = d[u] + 1;
          Q.push(v);
        }
      }
    }
    return vis[t];
  }
  int dfs(int u, int a) {
    if (u == t || a == 0) return a;
    int fl = 0, f;
    for (auto v : adj[u]) {
      if (d[u] + 1 == d[v] &&
          (f = dfs(v, min(a, capacity[u][v] - flow[u][v]))) > 0) {
        flow[u][v] += f;
        flow[v][u] -= f;
        fl += f;
        a -= f;
        if (a == 0) break;
      }
    }
    return fl;
  }

 public:
  Dinic(int n) {
    this->n = n;
    adj.resize(n);
    capacity.assign(n, vector<int>(n, 0)), flow.assign(n, vector<int>(n, 0));
    d.resize(n), cur.resize(n), vis.resize(n);
  }
  void addEdge(int u, int v, int c) {
    adj[u].push_back(v);
    adj[v].push_back(u);
    capacity[u][v] += c;
    capacity[v][u] += c;
  }
  int maxFlow(int s, int t) {
    this->s = s;
    this->t = t;
    int flow = 0;
    while (bfs()) {
      cur.assign(n, 0);
      flow += dfs(s, inf);
    }
    return flow;
  }
};
```

### Bipartite Match

```cpp
class Hungarian {
 private:
  vector<vector<int>> adj;
  vector<int> pre;
  vector<bool> vis;
  int n;
  bool dfs(int u) {
    for (auto v : adj[u]) {
      if (vis[v] == false) {
        vis[v] = true;
        if (pre[v] == -1 or dfs(pre[v])) {
          pre[v] = u;
          return true;
        }
      }
    }
    return false;
  }

 public:
  Hungarian(int n_) : n(n_) {
    adj.assign(n, vector<int>(0));
    pre.assign(n, -1);
    vis.assign(n, false);
  }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  int match() {
    int sum = 0;
    for (int i = 0; i < n; i++) {
      vis.assign(n, false);
      if (dfs(i)) sum++;
    }
    return sum;
  }
};
```


### Euler
```cpp
int V;                    // number of nodes
int E;                    // number of edges
vector<vector<int>> adj;  // the graph
vector<int> indegree;     // store indegree of each vertex
deque<int> cycle;

void find_cycle(int u) {
  while (adj[u].size()) {
    int v = adj[u].back();
    adj[u].pop_back();
    find_cycle(v);
  }
  cycle.push_front(u);
}

void euler_cycle() {
  // test if solution can exist
  for (int i = 0; i < V; i++)
    if (indegree[i] != adj[i].size()) {
      cout << "IMPOSSIBLE" << endl;
      exit(0);
    }

  // start anywhere
  find_cycle(0);  // populate cycle
  // test against disconnected graphs
  if (cycle.size() != E + 1) {
    cout << "IMPOSSIBLE" << endl;
    exit(0);
  }
  for (auto v : cycle) cout << v << " ";
  cout << endl;
}
```

### LCA
```cpp
#include <iostream>
#include <vector>

using namespace std;

class UnionFind {
 private:
  vector<int> parent, rank;

 public:
  UnionFind(int N) {
    rank.assign(N, 0);
    parent.assign(N, 0);
    for (int i = 0; i < N; i++) parent[i] = i;
  }
  int findSet(int i) {
    if (parent[i] == i)
      return i;
    else  // path compression
      return parent[i] = findSet(parent[i]);
  }
  bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }
  void unionSet(int i, int j) {
    if (!isSameSet(i, j)) {
      int x = findSet(i), y = findSet(j);
      if (rank[x] > rank[y])
        parent[y] = x;
      else {
        parent[x] = y;
        if (rank[x] == rank[y]) rank[y]++;
      }
    }
  }
};

void dfs(int u, vector<vector<int>>& adj, UnionFind& UF, vector<int>& parent,
         vector<vector<int>>& queries, vector<bool>& visited) {
  for (auto v : adj[u]) {
    dfs(v, adj, UF, parent, queries, visited);
    UF.unionSet(u, v);
    parent[UF.findSet(u)] = u;
  }
  visited[u] = true;
  for (auto v : queries[u])
    if (visited[v])
      cout << "LCA of " << u << " and " << v << " is " << parent[UF.findSet(v)]
           << endl;
}

int main() {
  int V, E, root, Q;
  cin >> V >> E >> root >> Q;
  auto UF = UnionFind(V);
  vector<int> parent(V);
  vector<vector<int>> adj(V);
  vector<vector<int>> queries(V);
  vector<bool> visited(V, false);
  int u, v;
  for (int i = 0; i < E; i++) {
    cin >> u >> v;
    adj[u].push_back(v);
  }
  for (int i = 0; i < Q; i++) {
    cin >> u >> v;
    queries[u].push_back(v);
    queries[v].push_back(u);
  }
  for (int i = 0; i < V; i++) parent[i] = i;

  dfs(root, adj, UF, parent, queries, visited);
}

```

## String

### KMP

```cpp
  string s, t;
  int n, m;

  cin >> s >> t;
  n = s.size();
  m = t.size();

  vector<int> lsp(m, 0);
  for (int i = 1, prev = 0; i < m;) {
    if (t[i] == t[prev]) {
      prev++;
      lsp[i] = prev;
      i++;
    } else if (prev == 0) {
      lsp[i] = 0;
      i++;
    } else {
      prev = lsp[prev - 1];
    }
  }

  int start = 0, len = 0;
  while (start + len < n) {
    while (len >= m || s[start + len] != t[len]) {
      if (len == 0) {
        start++;
        len = -1;
        break;
      }
      int skip = len - lsp[len - 1];
      start += skip;
      len -= skip;
    }
    len++;
    if (len == m) cout << "t matches s at " << start << "\n";
  }
```

### Trie

```cpp
struct trie {
  bool isEndOfString = false;
  map<char, trie *> edges;
  void insert(string &s, int i = 0) {
    if (i == s.length()) {
      isEndOfString = true;
      return;
    }
    if (edges.count(s[i]) == 0) edges[s[i]] = new trie;
    edges[s[i]]->insert(s, i + 1);
  }
  bool contains(string &s, int i = 0) {
    if (i == s.length()) return isEndOfString;
    return edges.count(s[i]) > 0 && edges[s[i]]->contains(s, i + 1);
  }
};
```

### Ahocorasick

```cpp
struct ACTrie {
  map<char, ACTrie *> edges;
  ACTrie *lsp = nullptr;
  int cnt = 0;

  void insert(string &t, int i = 0) {
    if (i == t.length()) {
      cnt++;
      return;
    }
    if (edges.count(t[i]) == 0) edges[t[i]] = new ACTrie;
    edges[t[i]]->insert(t, i + 1);
  }
  void search(string &s, int i = 0, bool print = true) {
    if (print) cout << cnt << " matches ending at " << i - 1 << "\n";
    if (i == s.length()) return;
    if (edges.count(s[i]) == 0) {
      if (lsp == nullptr)
        search(s, i + 1, true);
      else
        lsp->search(s, i, false);
    } else {
      edges[s[i]]->search(s, i + 1, true);
    }
  }
};

void preprocess(ACTrie &root) {
  queue<ACTrie *> q;
  root.lsp = nullptr;
  q.push(&root);
  while (!q.empty()) {
    ACTrie *u = q.front();
    q.pop();
    for (auto it : u->edges) {
      char c = it.first;
      ACTrie *v = it.second;
      ACTrie *l = u->lsp;
      while (l != nullptr && l->edges.count(c) == 0) l = l->lsp;
      if (l == nullptr) {
        v->lsp = &root;
      } else {
        v->lsp = l->edges[c];
        v->cnt += v->lsp->cnt;
      }
      q.push(v);
    }
  }
}
```

### Math

```cpp
int64_t fexp(int64_t m, int64_t n, int64_t p) {
  if (n == 0)
    return 1;
  else if (n % 2 == 1)
    return (m * fexp(m, n - 1, p)) % p;
  else {  // n is even
    int64_t r = fexp(m, n / 2, p);
    return (r * r) % p;
  }
}

int64_t comb(int64_t n, int64_t k) {
  const int64_t p = 998244353;

  int64_t ni = 1, ki = 1;
  for (int64_t i = n; i > n - k; i--) {
    ni = (ni * i) % p;
  }

  for (int64_t i = k; i > 0; i--) {
    ki = (ki * i) % p;
  }
  int kni = fexp(ki, p - 2, p);

  return (ni * kni) % p;
}

```

### Parser

```cpp
struct tree {
  vector<shared_ptr<tree>> children;
  char op = 0;
  long long value = -1;

  tree(vector<shared_ptr<tree>> children, char op)
      : children(children), op(op) {}
  explicit tree(long long value) : value(value) {}
};

string expression;
int pos;

shared_ptr<tree> parse_stree();

shared_ptr<tree> parse_tree() {
  vector<shared_ptr<tree>> children;
  char op = expression[pos];
  ++pos;
  while (expression[pos] != ')') {
    ++pos;
    children.push_back(parse_stree());
  }
  ++pos;
  return make_shared<tree>(children, op);
}

shared_ptr<tree> parse_int() {
  int result = 0;
  while (expression[pos] >= '0' && expression[pos] <= '9') {
    result = 10 * result + expression[pos] - '0';
    ++pos;
  }
  return make_shared<tree>(result);
}

shared_ptr<tree> parse_stree() {
  char c = expression[pos];
  if (c == '+' || c == '*') {
    return parse_tree();
  } else {
    return parse_int();
  }
}

int main() {
  ios_base::sync_with_stdio(false);

  getline(cin, expression);
  shared_ptr<tree> t = parse_stree();

  // Do something with t
}
```