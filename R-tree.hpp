#include <iostream>
#include <cstdio>

using namespace std;
const int DIMENSION = 3; // R树维护的空间维数
const int M = 10; //最大分支数
class RTree {
private:
  class Rec { //n维空间矩形
  public:
    double bound[2 * DIMENSION];
    Rec() {
      memset(bound, 0, sizeof bound);
    }
    void init() {
      memset(bound, 0, sizeof bound);
    }
  };
  class RtreeNode { //Rtree结点
  public:
    int count;
    int level; //0 -> leaf
    RtreeBranch branch[M];
    RtreeNode() {
      count = 0;
      level = 0;
      for (int i = 0; i < M; i++) {
        brach[i].init();
      }
    }
  };
  class RtreeBranch { //结点中的分支
  public:
    Rec mbr;
    RtreeNode *son;
    RtreeBranch() {}
    void init() {
      mbr.init();
      son = NULL;
    }
  };
  class RtreeRoot {
  public:
    RtreeNode *rnode;
    RtreeBranch branch[M + 1];
    int count;
    Rec overall;
    double overallArea;
  };

  typedef RtreeNode* Node;
  typedef RtreeRoot* Root;

  Root R_root;

  Rtree() {
    R_root = NULL;
  }

  Rtree(const Rtree &other) {

  }

  ~Rtree() {

  }

  Rec NullRec();

  double RecArea(Rec *mbr); //n维矩形最小边界球体积

  Rec CombineRec(Rec *rc1, Rec *rc2); //返回包含rc1和rc2的最小矩形

  bool RecOverlap(Rec *rc1, Rec *rc2); //返回rc1和rc2是否overlap

  void SplitNode(Root root, Node node, RtreeBranch *br, Node *new_node); //将node分裂成node和new_node

  Rec CoverRec(Node node); //返回包含node中所有branch的最小矩形

  int ChooseBranch(Rec *mbr, Node node); //返回node的所有branch中加入mbr之后覆盖矩形增量最小的branch编号

  int RTreeAddBranch(Root root, RtreeBranch *br, Node node, Node *new_node);

  void CutNode(Node node, int i);

  void Destroy(Node node); //释放空间,析构用

  int RTreeSearch(Root root, Rec *mbr, Rec *target);

  int RTreeInsert(Root root, Rec *data_mbr, Rec *data, int level);

  int RTreeDelete(Root root, Rec *data_mbr, Rec *data);

  //..split 过程中还需要一些函数，未加入


public:
  void insert(vector<double> rec) { //插入一个区域(元素)

  }

  void delete(vector<double> rec) { //删除一个区域(元素)

  }

  int search(vector<double> rec) { // 查询一个区域内有多少元素

  }
};
