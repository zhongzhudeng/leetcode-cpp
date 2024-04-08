#include <vector>

#include <catch2/catch_all.hpp>
using namespace std;

class Solution {
 private:
  static constexpr int INF = 10000 * 101 + 1;

 public:
  int findCheapestPrice(int n, vector<vector<int>>& flights, int src, int dst,
                        int k) {
    vector<int> f(n, INF);

    f[src] = 0;
    int ans = INF;
    for (int t = 1; t <= k + 1; ++t) {
      vector<int> g(n, INF);
      for (auto&& flight : flights) {
        int j = flight[0], i = flight[1], cost = flight[2];
        g[i] = min(g[i], f[j] + cost);
      }
      f = std::move(g);
      ans = min(ans, f[dst]);
    }
    return (ans == INF ? -1 : ans);
  }
};

TEST_CASE("787. Cheapest Flights Within K Stops") {
  Solution s;
  int n = 3, src = 0, dst = 2, k = 1, output = 200;
  vector<vector<int>> flights = {{0, 1, 100}, {1, 2, 100}, {0, 2, 500}};
  REQUIRE(s.findCheapestPrice(n, flights, src, dst, k) == output);
  k = 0, output = 500;
  REQUIRE(s.findCheapestPrice(n, flights, src, dst, k) == output);
}