#include <iostream>
#include <cstdio>

using namespace std;
const int DIMENSION = 3; // R树维护的空间维数
const int M = 10; //最大分支数
const double eps = 1e-8;
const double xpi[] = {
 0.000000,  /* dimension   0 */
 2.000000,  /* dimension   1 */
 3.141593,  /* dimension   2 */
 4.188790,  /* dimension   3 */
 4.934802,  /* dimension   4 */
 5.263789,  /* dimension   5 */
 5.167713,  /* dimension   6 */
 4.724766,  /* dimension   7 */
 4.058712,  /* dimension   8 */
 3.298509,  /* dimension   9 */
 2.550164,  /* dimension  10 */
 1.884104,  /* dimension  11 */
 1.335263,  /* dimension  12 */
 0.910629,  /* dimension  13 */
 0.599265,  /* dimension  14 */
 0.381443,  /* dimension  15 */
 0.235331,  /* dimension  16 */
 0.140981,  /* dimension  17 */
 0.082146,  /* dimension  18 */
 0.046622,  /* dimension  19 */
 0.025807,  /* dimension  20 */
};
int sign(double x) {
  return x < -eps ? -1 : x > eps;
}
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
        branch[i].init();
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
    Destroy(R_root->rnode);
    delete R_root;
    R_root = NULL;
  }
  int sons(Node node) { //返回一个结点的儿子数， 叶子和非叶子会有区别
  }

  Rec NullRec();

  double RecArea(Rec *mbr); //n维矩形最小边界球体积 OK

  Rec CombineRec(Rec *rc1, Rec *rc2); //返回包含rc1和rc2的最小矩形 OK

  bool RecOverlap(Rec *rc1, Rec *rc2); //返回rc1和rc2是否overlap OK

  void SplitNode(Root root, Node node, RtreeBranch *br, Node &new_node); //将node分裂成node和new_node 

  Rec CoverRec(Node node); //返回包含node中所有branch的最小矩形 OK

  int ChooseBranch(Rec *mbr, Node node); //返回node的所有branch中加入mbr之后覆盖矩形增量最小的branch编号 OK

  int RTreeAddBranch(Root root, RtreeBranch *br, Node node, Node &new_node);

  void CutNode(Node node, int i); //

  void Destroy(Node node); //释放空间,析构用 OK

  int RtreeSearch(Root root, Rec *target);

  int RtreeInsert(Root root, Rec *data, int level);

  int RtreeDelete(Root root, Rec *data);

  //..split 过程中还需要一些函数，未加入


  Rec NullRec() {
    Rec ret;
    ret.bound[0] = 1;
    ret.bound[DIMENSION] = -1;
    for (int i = 1; i < DIMENSION; i++) {
      ret.bound[i] = ret.bound[i + DIMENSION] = 0;
    }
    return ret;
  }

  double RecArea(Rec *mbr) {
    double ret = 0;
    if (invalid_rec(mbr)) {
      return 0;
    }
    double sumr = 0;
    for (int i = 0; i < DIMENSION; i++) {
      double len = (buond[i + DIMENSION] - bound[i]) / 2;
      sumr += len * len;
    }
    double radius = sqrt(sumr);
    ret = pow(radius, DIMENSION) * xpi[DIMENSION];
    return ret;
  }

  Rec CombineRec(Rec *rc1, Rec *rc2) {
    Rec ret;
    if (invalid_rec(rc1)) {
      return *rc2;
    }
    if (invalid_rec(rc2)) {
      return *rc1;
    }
    for (int i = 0; i < DIMENSION; i++) {
      ret.bound[i] = min(rc1->bound[i], rc2->bound[i]);
    }
    for (int i = DIMENSION; i < 2 * DIMENSION; i++) {
      ret.bound[i] = max(rc1->bound[i], rc2->bound[i]);
    }
    return ret;
  }

  bool RecOverlap(Rec *rc1, Rec *rc2) {
    if (invalid_rec(rc2) || invalid_rec(rc1)) {
      return 0;
    }
    Rec tmp;
    for (int i = 0; i < DIMENSION; i++) {
      tmp.bound[i] = max(rc1->bound[i], rc2->bound[i]);
      tmp.bound[i + DIMENSION] = min(rc1->bound[i + DIMENSION], rc2->bound[i + DIMENSION]);
    }
    for (int i = 0; i < DIMENSION; i++) {
      if (tmp.bound[i] > tmp.bound[i + DIMENSION]) {
        return false;
      }
    }
    return true;
  }
  
  Rec CoverRec(Node node) {
    Rec ret;
    bool flg = 1;
    for (int i = 0; i < sons(node); i++) {
      if (flg) {
        ret = (node->branch[i]).mbr;
        flg = 0;
      } else {
        ret = CombineRec(&ret, &((node->branch[i]).mbr));
      }
    }
    return ret;
  }

  bool isbetter(double former, double increment, double bestinc, double bestfor) {
    if (sign(bestinc) < 0) {
      return true;
    }
    if (sign(increment - bestinc) < 0) {
      return true;
    }
    if (sign(increment - bestinc) == 0) {
      return sign(former - bestfor) < 0;
    }
    return 0;
  }

  int ChooseBranch(Rec *mbr, Node node) {
    int best = -1;
    double bestinc = -1;
    double bestfor = -1;
    for (int i = 0; i < sons(node); i++) {
      Rec tmp = CombineRec(mbr, &((node->branch[i]).mbr));
      double former = RecArea(&((node->branch[i]).mbr));
      double increment = RecArea(&tmp) - former; 
      if (best == -1 || isbetter(former, increment, bestinc, bestfor)) {
        best = i;
        bestinc = increment;
        bestfor = former;
      } 
    }
    return best;
  }

  void Destroy(Node node) {
    if (node == NULL) {
      return;
    }
    for (int i = 0; i < sons(node); i++) {
      Destroy(node->branch[i].son);
    }
    delete node;
    node = NULL;
  }

  int RtreeSearch(Node node, Rec *target) {
    int ret = 0;
    for (int i = 0; i < sons(node); i++) {
      if (overlap(&((node->branch[i]).mbr), target)) {
        if (node->level > 0) {
          ret += search(node->branch[i].son, target);
        } else {
          ret++;
        }
      }
    }
    return ret;
  }

  int RtreeSearch(Root root, Rec *target) {
    return RtreeSearch(root->rnode, target);
  }

  

public:
  void insert(vector<double> rec) { //插入一个区域(元素)
  }

  void delete(vector<double> rec) { //删除一个区域(元素)
  }

  int search(vector<double> rec) { // 查询一个区域内有多少元素
    Rec data;
    if (rec.size() != DIMENSION * 2) {
      throw invalid_input();
    }
    for (int i = 0; i < 2 * DIMENSION) {
      data.bound[i] = rec[i];
    }
    return RtreeSearch(root, &data);
  }
};
