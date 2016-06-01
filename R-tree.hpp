#ifndef SJTU_RTREE
#define SJTU_RTREE

#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>

const int DIMENSION = 3; // R树维护的空间维数
const int M = 10; //最大分支数



class RTree
{
protected:
	struct Rec;
	struct RtreeNode;
	struct RtreeBranch;
	struct RtreeRoot;
public:
	RTree()
	{
		R_root = 0;
	}

	RTree(const RTree &other)
	{

	}

	~RTree()
	{
		Destroy(R_root->rnode);
		delete R_root;
		R_root = 0;
	}

	void insert(const std::vector<double> &rec)
	{ //插入一个区域(元素)
	}

	void remove(const std::vector<double> &rec)
	{ //删除一个区域(元素)
		if(rec.size() != DIMENSION * 2)
			throw std::invalid_argument("The number of input doesn't match dimension");
		Rec target(rec);
		std::vector<Rec> vRebuild;
		remove(R_root->rnode, target, vRebuild);
		assert(R_root->rnode->count >= 1);
		if (R_root->rnode->count == 1)
		{
			RtreeNode *tmp = R_root->rnode;
			R_root->rnode = R_root->rnode->branch[0].child;
			delete tmp;
		}
		//TODO
	}

	int search(const std::vector<double> &rec)
	{ // 查询一个区域内有多少元素
		Rec data;
		if (rec.size() != DIMENSION * 2)
		{
			throw std::invalid_argument("The number of input doesn't match dimension");
		}
		for (int i = 0; i < 2 * DIMENSION;i++)
		{
			data.bound[i] = rec[i];
		}
		return RtreeSearch(R_root, &data);
	}


protected:
	struct Rec
	{ //n维空间矩形
		double bound[2 * DIMENSION];
		Rec()
		{
			memset(bound, 0, sizeof(bound));
		}
		Rec(const std::vector<double> &other)
		{
			assert(other.size() == 2 * DIMENSION);
			for (int i = 0; i < 2 * DIMENSION; i++)
				bound[i] = other[i];
		}
		void init()
		{
			memset(bound, 0, sizeof(bound));
		}
		bool is_valid() const
		{
			for (int i = 0; i < DIMENSION; i++)
				if (bound[i] > bound[i + DIMENSION])
					return false;
			return true;
		}
		bool operator==(const Rec &other) const
		{
			for (int i = 0; i < 2 * DIMENSION; i++)
			{
				if (sign(bound[i] - other.bound[i]) != 0)
					return false;
			}
			return true;
		}
	};
	struct RtreeBranch
	{ //结点中的分支
		Rec mbr;
		RtreeNode *child;
		RtreeBranch() {}
		void init()
		{
			mbr.init();
			child = 0;
		}
	};
	struct RtreeNode
	{ //Rtree结点
		int count;
		int level; //0 -> leaf
		RtreeBranch branch[M];
		RtreeNode()
		{
			count = 0;
			level = -1;
			for (int i = 0; i < M; i++)
			{
				branch[i].init();
			}
		}
	};
	struct RtreeRoot
	{
		RtreeNode *rnode;
	};

	

	typedef RtreeNode *Node;
	typedef RtreeRoot *Root;

	RtreeRoot *R_root;

	
	/*int sons(Node node)
	{ //返回一个结点的儿子数， 叶子和非叶子会有区别
	}*/

	//Rec NullRec(); //OK

	//double RecArea(Rec *mbr); //n维矩形最小边界球体积 OK

	//Rec CombineRec(Rec *rc1, Rec *rc2); //返回包含rc1和rc2的最小矩形 OK

	//bool RecOverlap(Rec *rc1, Rec *rc2); //返回rc1和rc2是否overlap OK

	void SplitNode(Root root, Node node, RtreeBranch *br, Node *new_node); //将node分裂成node和new_node 

	//Rec CoverRec(Node node); //返回包含node中所有branch的最小矩形 OK

	//int ChooseBranch(Rec *mbr, Node node); //返回node的所有branch中加入mbr之后覆盖矩形增量最小的branch编号 OK

	//Node ChooseLeaf(Rec *mbr, Root root);

	//bool AddBranch(Root root, RtreeBranch *br, Node node, Node *new_node); // OK

	//void DeletBranch(Node node, int i); // OK

	void Destroy(Node node); //释放空间,析构用 OK

	int RtreeSearch(Root root, Rec *target); // OK

	//int RtreeInsert(Root root, Rec *data, int level); // OK

	//bool RtreeInsert(Root root, Rec *mbr, RtreeNode *node, RtreeNode* &new_node, int level);

	//int RtreeDelete(Root root, Rec *data);

	//..split 过程中还需要一些函数，未加入
	//然而我现在全写在了split里面了 by Steiner

	Rec NullRec()
	{
		Rec ret;
		ret.bound[0] = 1;
		ret.bound[DIMENSION] = -1;
		for (int i = 1; i < DIMENSION; i++)
		{
			ret.bound[i] = ret.bound[i + DIMENSION] = 0;
		}
		return ret;
	}

	double RecArea(const Rec &mbr)
	{
		if (!mbr.is_valid())
			return 0;
		static const double xpi = vratio(DIMENSION);
		double ret = 0;
		double sumr = 0;
		for (int i = 0; i < DIMENSION; i++)
		{
			double len = (mbr.bound[i + DIMENSION] - mbr.bound[i]) / 2;
			sumr += len * len;
		}
		double radius = sqrt(sumr);
		ret = pow(radius, DIMENSION) * xpi;
		return ret;
	}

	Rec CombineRec(const Rec &rc1, const Rec &rc2)
	{
		if (!rc1.is_valid())
			return rc2;
		if (!rc2.is_valid())
			return rc1;
		Rec ret;
		for (int i = 0; i < DIMENSION; i++)
		{
			ret.bound[i] = std::min(rc1.bound[i], rc2.bound[i]);
		}
		for (int i = DIMENSION; i < 2 * DIMENSION; i++)
		{
			ret.bound[i] = std::max(rc1.bound[i], rc2.bound[i]);
		}
		return ret;
	}

	static bool RecOverlap(const Rec *rc1, const Rec *rc2)
	{
		assert(rc1->is_valid() && rc2->is_valid());
		Rec tmp;
		for (int i = 0; i < DIMENSION; i++)
		{
			tmp.bound[i] = std::max(rc1->bound[i], rc2->bound[i]);
			tmp.bound[i + DIMENSION] = std::min(rc1->bound[i + DIMENSION], rc2->bound[i + DIMENSION]);
		}
		for (int i = 0; i < DIMENSION; i++)
		{
			if (tmp.bound[i] > tmp.bound[i + DIMENSION])
			{
				return false;
			}
		}
		return true;
	}

	Rec CoverRec(const RtreeNode *node)
	{
		Rec ret;
		bool flag = 1;
		for (int i = 0; i < node->count; i++)
		{
			if (flag)
			{
				ret = (node->branch[i]).mbr;
				flag = 0;
			}
			else
			{
				ret = CombineRec(ret, node->branch[i].mbr);
			}
		}
		return ret;
	}



	void Destroy(RtreeNode *node)
	{
		if (!node)
		{
			return;
		}
		for (int i = 0; i < M; i++)
		{
			Destroy(node->branch[i].son);
		}
		delete node;
		node = 0;
	}



	int ChooseBranch(Rec *mbr, Node node)
	{
		int best = -1;
		double bestinc = -1;
		double bestfor = -1;
		for (int i = 0; i < node->count; i++)
		{
			Rec tmp = CombineRec(*mbr, ((node->branch[i]).mbr));
			double former = RecArea((node->branch[i]).mbr);
			double increment = RecArea(tmp) - former;
			if (best == -1 || isbetter(former, increment, bestinc, bestfor))
			{
				best = i;
				bestinc = increment;
				bestfor = former;
			}
		}
		return best;
	}

	bool AddBranch(RtreeBranch *br, RtreeNode *node, Node *new_node)
	{
		if (node->count < M)
		{
			node->branch[node->count++] = *br;
			return false;
		}
		SplitNode(node, br, new_node);
		return true;
	}

	void DeleteBranch(RtreeNode *node, int i)
	{
		assert(node->count >= 1);
		node->branch[i].child = nullptr;
		std::swap(node->branch[i], node->branch[node->count-1]);
		node->count--;
	}

	bool RtreeInsert(Rec *mbr, RtreeNode *node, RtreeNode* &new_node, int level)
	{
		assert(level == node->level);
		if (!level)
		{
			new_node = new RtreeNode();
			new_node->branch[new_node->count++].mbr = *mbr;
			return true;
		}
		int chosen = ChooseBranch(mbr, node);
		bool res = RtreeInsert(mbr, node->branch[chosen].child, new_node, level - 1);
		node->branch[chosen].mbr = CoverRec(node->branch[chosen].child);
		if (res)
		{
			RtreeBranch nbra;
			nbra.mbr = CoverRec(new_node);
			nbra.child = new_node;
			return AddBranch(&nbra, node, &new_node);
		}
		return false;
	}

	int RtreeSearch(const RtreeNode *node, const Rec *target)
	{
		int ret = 0;
		for (int i = 0; i < node->count; i++)
		{
			if (RecOverlap(&((node->branch[i]).mbr), target))
			{
				if (node->level > 0)
				{
					ret += RtreeSearch(node->branch[i].child, target);
				}
				else
				{
					ret++;
				}
			}
		}
		return ret;
	}

	int RtreeSearch(Root root, Rec *target)
	{
		return RtreeSearch(root->rnode, target);
	}

	void SplitNode(Node node, RtreeBranch *br, Node *new_node)
	{
		static Rec tmp[M + 1];
		static int taken[M + 1];
		static double tmpArea[M + 1];

		int level = node->level;
		for (int i = 0; i < M; i++)
		{
			tmp[i] = node->branch[i].mbr;
		}
		tmp[M] = br->mbr;
		Node nn = new RtreeNode;
		Rec uni = tmp[0];
		for (int i = 1; i < M + 1; i++)
		{
			uni = CombineRec(uni, tmp[i]);
		}
		for (int i = 0; i < M + 1; i++)
		{
			tmpArea[i] = RecArea(tmp[i]);
		}
		double area = RecArea(uni);

		int seed0 = -1, seed1 = -1;
		for (int i = 0; i < M; i++)
		{
			for (int j = i + 1; j < M + 1; j++)
			{
				if (area - tmpArea[i] - tmpArea[j] > worst)
				{
					worst = area - tmpArea[i] - tmpArea[j];
					seed0 = i;
					seed1 = j;
				}
			}
		}
		for (int i = 0; i < M + 1; i++)
		{
			taken[i] = -1;
		}

		taken[seed0] = 0;
		taken[seed1] = 1;

		Rec group[2];
		double garea[2];
		group[0] = tmp[seed0];
		group[1] = tmp[seed1];
		garea[0] = tmpArea[seed0];
		garea[1] = tmpArea[seed1];

		int count[2];
		count[0] = count[1] = 1;
		int rest = M - 2;
		int bigger = 0;
		double best = -1;
		std::pair<int, int> which = std::make_pair(-1, -1);

		for (; rest > 0 && count[bigger] < (M + 1) / 2; )
		{
			for (int i = 0; i < M + 1; i++)
			{
				if (taken[i] > -1)
				{
					continue;
				}
				double d0 = RecArea(CombineRec(group[0], tmp[i])) - garea[0];
				double d1 = RecArea(CombineRec(group[1], tmp[i])) - garea[1];
				if (fabs(d0 - d1) > best)
				{
					best = fabs(d0 - d1);
					which = std::make_pair((sign(d0 - d1) < 0 ? 0 : 1), i);
				}
			}
			taken[which.second] = which.first;
			rest--;
			count[which.first]++;
			if (count[which.first] > count[bigger])
			{
				bigger = which.first;
			}
			group[which.first] = CombineRec(group[which.first], tmp[which.second]);
			garea[which.first] = RecArea(group[which.first]);
		}
		if (rest > 0)
		{
			for (int i = 0; i < M + 1; i++)
			{
				if (taken[i] == -1)
				{
					group[bigger ^ 1] = CombineRec(group[bigger ^ 1], tmp[i]);
					taken[i] = (bigger ^ 1);
				}
			}
		}
		*new_node = new RtreeNode;
		(*new_node)->level = level;
		node->count = 0;
		for (int i = 0; i < M; i++)
		{
			node->branch[i].init();
		}
		for (int i = 0; i < M + 1; i++)
		{
			if (taken[i] == 0)
			{
				node->branch[(node->count)++].mbr = tmp[i];
			}
			else
			{
				(*new_node)->branch[(*new_node)->count++].mbr = tmp[i];
			}
		}
		
	}

	//vRebuild: the rects (leaf node) the should be reinsert
	//Return if the node o should be deleted
	bool remove(RtreeNode *o, const Rec &rec, std::vector<Rec> &vRebuild)
	{
		assert(o);
		if (o->level == 0)
		{
			assert(o->count == 1);
			if (rec == o->branch[0].mbr)
			{
				Destroy(o);
				return true;
			}
			return false;
		}
		int i = 0;
		while (i < o->count)
		{
			if (RecOverlap(&rec, &o->branch[i].mbr))
			{
				++i;
				continue;
			}
			if (!remove(o->branch[i].child, rec, vRebuild))
			{
				++i;
				continue;
			}
			DeleteBranch(o, i);
		}
		if (o->count < M / 2)
		{
			enumLeaf(o, vRebuild);
			Destroy(o);
			return true;
		}
		return false;
	}

private:
	static constexpr double eps = 1e-8;
	static constexpr double Pi = 3.141592653589793238;
	static constexpr double factorial(int n)
	{
		return n == 0 ? 1 : n * factorial(n - 1);
	}
	static constexpr double ipow(double a, int n)
	{
		return n == 0 ? 1 : a * ipow(a, n - 1);
	}
	static constexpr double vratio(int dim)
	{
		return dim % 2 ? 2 * factorial(dim / 2)*ipow(4 * Pi, dim / 2) / factorial(dim) : ipow(Pi, dim / 2) / factorial(dim / 2);
	}
	static constexpr int sign(double x)
	{
		return x < -eps ? -1 : x > eps;
	}

	static bool isbetter(double former, double increment, double bestinc, double bestfor)
	{
		if (sign(bestinc) < 0)
		{
			return true;
		}
		if (sign(increment - bestinc) < 0)
		{
			return true;
		}
		if (sign(increment - bestinc) == 0)
		{
			return sign(former - bestfor) < 0;
		}
		return 0;
	}

	//add o's all leaf nodes to vec
	void enumLeaf(RtreeNode *o, std::vector<Rec> &vec)
	{
		if (o->level == 0)
		{
			vec.push_back(o->branch[0].mbr);
			return;
		}
		for (int i = 0; i < o->count; i++)
			enumLeaf(o->branch[i].child, vec);
	}
};

#endif	//SJTU_RTREE