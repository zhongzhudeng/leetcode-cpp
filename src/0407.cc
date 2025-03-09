
#include <algorithm>
#include <queue>
#include <vector>
using namespace std;

typedef pair<int, int> pii;

class Solution {
public:
  int trapRainWater(vector<vector<int>> &heightMap) {
    size_t m = heightMap.size(), n = heightMap[0].size();
    if (m <= 2 || n <= 2)
      return 0;
    int ans = 0;
    std::vector<pii> pqv;
    vector<vector<bool>> visit(m, vector<bool>(n, false));
    for (int i = 0; i < m; i++) {
      pqv.push_back({heightMap[i][0], i * n});
      visit[i][0] = true;
      pqv.push_back({heightMap[i][n - 1], i * n + n - 1});
      visit[i][n - 1] = true;
    }
    for (int j = 1; j < n - 1; j++) {
      pqv.push_back({heightMap[0][j], j});
      visit[0][j] = true;
      pqv.push_back({heightMap[m - 1][j], (m - 1) * n + j});
      visit[m - 1][j] = true;
    }
    int res = 0;
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq(
        pqv.begin(), pqv.end());
    std::array<int, 5> dirs = {-1, 0, 1, 0, -1};
    while (not pq.empty()) {
      pii curr = pq.top();
      pq.pop();
      for (int k = 0; k < 4; k++) {
        int nx = curr.second / n + dirs[k];
        int ny = curr.second % n + dirs[k + 1];
        if (nx >= 0 and nx < m and ny >= 0 and ny < n and !visit[nx][ny]) {
          if (heightMap[nx][ny] < curr.first)
            res += curr.first - heightMap[nx][ny];
          visit[nx][ny] = true;
          pq.push({std::max(heightMap[nx][ny], curr.first), nx * n + ny});
        }
      }
    }
    return res;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("407. Trapping Rain Water 2", "[0407]") {
  Solution s;
  vector<vector<int>> in = {
      {1, 4, 3, 1, 3, 2}, {3, 2, 1, 3, 2, 4}, {2, 3, 3, 2, 3, 1}};
  int ans = 4;
  REQUIRE(s.trapRainWater(in) == ans);

  in = {{3, 3, 3, 3, 3},
        {3, 2, 2, 2, 3},
        {3, 2, 1, 2, 3},
        {3, 2, 2, 2, 3},
        {3, 3, 3, 3, 3}};
  ans = 10;
  REQUIRE(s.trapRainWater(in) == ans);
}