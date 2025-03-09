#include <algorithm>
#include <string>
using namespace std;

class Solution {
public:
  string maximumBinaryString(string binary) {
    size_t n = binary.size(), i = binary.find('0');
    if (i == std::string::npos)
      return binary;

    size_t zeros = count(binary.begin(), binary.end(), '0');
    std::string s(n, '1');
    s[i + zeros - 1] = '0';
    return s;
  }
};

#include <catch2/catch_test_macros.hpp>
TEST_CASE("1702. Maximum Binary String After Change", "[1702]") {
  Solution s;
  std::string in = "000110";
  std::string ans = "111011";
  REQUIRE(s.maximumBinaryString(in) == ans);

  in = "01";
  ans = "01";
  REQUIRE(s.maximumBinaryString(in) == ans);
}