#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>  // memcpy
#include <search.h>  //  qsort
#include <math.h>
#include <time.h>    // time
#include <assert.h>  // assert

typedef unsigned char uchar;

// 向量数据结构
#define VECT_DIM 3    // 向量的维度
struct Vect {
	double arr[VECT_DIM];   // 向量的分量
	//double distance;        // 当前向量与聚类中心的距离
	size_t tpy;           // 用无符号整数表示不同的分类
};

void k_means_plusplus(Vect * observations, size_t nSetSize, size_t nClasses, size_t nDim, double epsilon = 0.01);
void k_means_plusplus(Vect * observations, size_t nSetSize,	Vect * centroids, size_t nClasses, size_t nDim, double epsilon = 0.01);
void printVects(Vect *vects, size_t nSetSize);