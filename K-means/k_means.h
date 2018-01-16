#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>  // memcpy
#include <search.h>  //  qsort
#include <math.h>
#include <time.h>    // time
#include <assert.h>  // assert


// �������ݽṹ
#define VECT_DIM 2    // ������ά��
struct Vect {
	double arr[VECT_DIM];   // �����ķ���
	//double distance;        // ��ǰ������������ĵľ���
	size_t tpy;           // ���޷���������ʾ��ͬ�ķ���
};

void k_means(Vect * observations, size_t nSetSize, size_t nClasses, size_t nDim, double epsilon);
void printVects(Vect *vects, size_t nSetSize);