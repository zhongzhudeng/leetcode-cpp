#include <algorithm>
#include <unordered_map>
#include <vector>
using namespace std;

class Solution {
public:
private:
  static constexpr int MAXX = 1 << 10;
  static constexpr int INFTY = INT_MAX / 2;

public:
  int minChanges(vector<int> &nums, int k) {
    int n = nums.size();
    vector<int> f(MAXX, INFTY);
    f[0] = 0;

    for (int i = 0; i < k; ++i) {
      unordered_map<int, int> cnt;
      int size = 0;
      for (int j = i; j < n; j += k) {
        ++cnt[nums[j]];
        ++size;
      }

      int t2min = *min_element(f.begin(), f.end());

      vector<int> g(MAXX, t2min);
      for (int mask = 0; mask < MAXX; ++mask)
        for (auto [x, countx] : cnt)
          g[mask] = min(g[mask], f[mask ^ x] - countx);

      for_each(g.begin(), g.end(), [&](int &val) { val += size; });
      f = std::move(g);
    }

    return f[0];
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1787. Make the XOR of All Segments Equal to Zero", "[1787]") {
  Solution s;
  vector<int> nums = {1, 2, 0, 3, 0};
  int k = 1;
  REQUIRE(s.minChanges(nums, k) == 3);
  nums = {3, 4, 5, 2, 1, 7, 3, 4, 7};
  k = 3;
  REQUIRE(s.minChanges(nums, k) == 3);
  nums = {1, 2, 4, 1, 2, 5, 1, 2, 6};
  k = 3;
  REQUIRE(s.minChanges(nums, k) == 3);
  nums = {26, 19, 19, 28, 13, 14, 6, 25, 28, 19, 0, 15, 25, 11};
  k = 11;
  REQUIRE(s.minChanges(nums, k) == 3);
}