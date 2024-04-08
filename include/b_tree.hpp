#pragma once

#include <initializer_list>
#include <queue>

namespace LeetCode {

struct TreeNode {
  int val;
  TreeNode *left;
  TreeNode *right;
  TreeNode(int x = 0) : val(x), left(nullptr), right(nullptr) {}
};

void print_tree(TreeNode *root);
bool is_same_tree(TreeNode *root1, TreeNode *root2);

class Tree {
public:
  TreeNode *root;

private:
  std::queue<int> q;

public:
  Tree(std::initializer_list<int> il);
  ~Tree();
  void print();
  bool operator==(const Tree &t2) const;

private:
  TreeNode *deserialize(std::queue<int> &in);
};
} // namespace LeetCode