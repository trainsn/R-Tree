#include "R-Tree.hpp"
#include <iostream>
#include <cstdio>
#include <algorithm>

int main()
{
	RTree rt;
	for (int i = 0; i < 100; i++)
	{
		double f = i;
		rt.insert({0,0,0,f,f,f});
	}
	for (int i = 10; i < 60; i++)
	{
		double f = i;
		rt.remove({0,0,0,f,f,f});
	}
	for (int i = 0; i < 100; i++)
	{
		double f = i;
		printf("%d\n", rt.search({f,f,f,100,100,100}));
	}
	return 0;
}
