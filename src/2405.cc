#include <string>
#include <string_view>
using namespace std;

class Solution {
public:
  int partitionString(string s) {
    int cnt = 0;
    size_t i = 0, j = 1, sz = s.length();
    while (j < sz) {
      string_view sv(s.begin() + i, s.begin() + j);
      char c = s[j];
      int fd = sv.find(s[j]);
      if (fd == std::string::npos)
        j++;
      else
        i = j, j++, cnt++;
    }
    cnt++;
    return cnt;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("2405. Optimal Partition of String", "[2405]") {
  Solution s;
  REQUIRE(s.partitionString("abacaba") == 4);
  REQUIRE(s.partitionString("ssssss") == 6);
}