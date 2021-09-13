# C++ 笔记

## class 的 Rule of Five

当类成员变量中不包含指针时，我们直接让编译器默认生成就可以
当类成员变量时存在指针时，就需要我们注意 浅拷贝 和 深拷贝

C++ 和构造函数相关的五个函数,当你显示的定义了其中的任意一个,你应该同时实现另外几个

* destructor（析构）
* copy constructor（复制构造）
* copy assignment operator（复制赋值运算符）
* move constructor(移动构造)
* move assignment operator（移动赋值）
  
语句

```cpp
Widget w(w0);
Widget w = w0;
Widget w; w = w0;
Widget w(std::move(w0));
Widget w; w = std::move(w0);
```

分别调用的是

```cpp
Widget(const Widget& w);
Widget(const Widget& w);
Widget &operator=(const Widget& w); 
Widget(Widget&& w); 
Widget &operator=(Widget&&);
```