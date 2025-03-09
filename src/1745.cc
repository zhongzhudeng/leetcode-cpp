#include <string>
#include <vector>
using namespace std;

class Solution {
public:
  bool checkPartitioning(string s) {
    size_t n = s.length();
    vector<vector<bool>> dp(n, vector<bool>(n));
    for (size_t len = 1; len < n; len++) {
      for (size_t start = 0; start <= n - len; start++) {
        size_t end = start + len - 1;
        if (len == 1)
          dp[start][end] = true;
        else if (len == 2)
          dp[start][end] = (s[start] == s[end]);
        else
          dp[start][end] = (s[start] == s[end] && dp[start + 1][end - 1]);
      }
    }
    for (size_t start = 1; start < n - 1; start++) {
      if (!dp[0][start - 1])
        continue;
      for (size_t end = start; end < n - 1; end++)
        if (dp[start][end] && dp[end + 1][n - 1])
          return true;
    }
    return false;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1745. Palindrome Partitioning IV", "[1745]") {
  Solution s;
  string in = "abcbdd";
  REQUIRE(s.checkPartitioning(in) == true);
  in = "bcbddxy";
  REQUIRE(s.checkPartitioning(in) == false);
}
