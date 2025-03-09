#include <list>
#include <vector>
using namespace std;

class MyHashMap {
  std::vector<std::list<std::pair<int, int>>> data;
  static const int base = 769;
  int hash(int key) { return key % base; }

public:
  MyHashMap() : data(base) {}

  void put(int key, int value) {
    int h = hash(key);
    for (auto it = data[h].begin(); it != data[h].end(); it++) {
      if (it->first == key) {
        it->second = value;
        return;
      }
    }
    data[h].emplace_back(key, value);
  }

  int get(int key) {
    int h = hash(key);
    for (auto it = data[h].begin(); it != data[h].end(); it++)
      if (it->first == key)
        return it->second;
    return -1;
  }

  void remove(int key) {
    int h = hash(key);
    for (auto it = data[h].begin(); it != data[h].end(); it++)
      if (it->first == key) {
        it = data[h].erase(it);
        return;
      }
  }
};

std::vector<std::string> run(std::vector<std::string> &cmd,
                             std::vector<std::vector<int>> &v) {
  std::vector<std::string> ans;
  ans.reserve(cmd.size());
  MyHashMap *mhm;
  char buf[64];
  for (int i = 0; i < cmd.size(); i++) {
    if (cmd[i] == "put") {
      mhm->put(v[i][0], v[i][1]);
      ans.push_back("null");
    } else if (cmd[i] == "get") {
      snprintf(buf, 64, "%d", mhm->get(v[i][0]));
      ans.push_back(buf);
    } else if (cmd[i] == "remove") {
      mhm->remove(v[i][0]);
      ans.push_back("null");
    } else {
      mhm = new MyHashMap;
      ans.push_back("null");
    }
  }
  delete mhm;
  return ans;
}

#include <catch2/catch_test_macros.hpp>
TEST_CASE("0706. Design HashMap", "[0706]") {
  std::vector<std::string> cmd = {"MyHashMap", "put", "put",    "get", "get",
                                  "put",       "get", "remove", "get"};
  std::vector<std::vector<int>> v = {{},     {1, 1}, {2, 2}, {1}, {3},
                                     {2, 1}, {2},    {2},    {2}};
  std::vector<std::string> ans = {"null", "null", "null", "1", "-1",
                                  "null", "1",    "null", "-1"};
  REQUIRE(run(cmd, v) == ans);
}