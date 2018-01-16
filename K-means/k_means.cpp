#include "k_means.h"


// 距离计算函数，计算矢量a与矢量b的坐标距离
double calcDistance(double *a, double *b, size_t dim) {
	double sum = 0.0;
	double diff;
	for (size_t i = 0; i < dim; i++) {
		diff = a[i] - b[i];
		sum += (diff * diff);
	}
	return sqrt(sum);
}

// 找到数组a中最小值，返回其在a中的下标
// 接受类型为double, int, unsigned等支持比较运算符的基类
template <typename T>
size_t findMin(const T *a, size_t nLen) {
	size_t indexOfMin = 0;  // 最小值所在的下标
	for (size_t i = 1; i < nLen; i++) {
		if (a[i] < a[indexOfMin]) {
			indexOfMin = i;
		}
	}
	return indexOfMin;
}

// 找到数组a中最大值，返回其在a中的下标
// 接受类型为double, int, unsigned等支持比较运算符的基类
template <typename T>
size_t findMax(const T *a, size_t nLen) {
	size_t indexOfMax = 0;  // 最小值所在的下标
	for (size_t i = 1; i < nLen; i++) {
		if (a[i] > a[indexOfMax]) {
			indexOfMax = i;
		}
	}
	return indexOfMax;
}

// bsearch的回调函数（在markType中用到）
int cmpVects_arr(Vect *key, Vect *datum) {
	int diff = memcmp(key->arr, datum->arr, VECT_DIM * sizeof(double));
	return diff;
}

// bsearch的回调函数（在k_means中用到）
int cmpInt(int *key, int *datum) {
	int diff = *key - *datum;
	return diff;
}

// qsort的回调函数
int cmpVects_typ(Vect *elem1, Vect *elem2) {
	int diff = (int)(elem1->tpy) - (int)(elem2->tpy);
	return diff;
}


// 为observations各元素标记类型，并修改其typ分量
// observations为待分类样本集
// nSetSize为样本集中的元素个数
// nClasses为期望分类的分类总数
// centroidsInd是表示当前聚类中心的一维数组，其大小为nClasses。centroidsInd各元素为聚类中心在数组observations中的数组下标
void markType(Vect *observations, size_t nSetSize, Vect *centroids, size_t nClasses) {
	// 保证centroids是按照centroids[].tpy升序排序的
	qsort(centroids, nClasses, sizeof(Vect), (int (*)(const void*, const void*))cmpVects_typ);

	double *distancesToCentroids = new double[nClasses];
	for (size_t i = 0; i < nSetSize; i++) {
		//Vect *isDup = (Vect *)bsearch(observations + i, centroids, nClasses, sizeof(Vect), (int(*)(const void*, const void*))cmpVects_arr);
		////int diff = memcmp(observations + i, centroids + i, sizeof(Vect));
		//if (isDup) {  // 当前的observation和centroids是同一个点
		//	continue;
		//}
		for (size_t j = 0; j < nClasses; j++) {
			distancesToCentroids[j] = calcDistance(observations[i].arr, centroids[j].arr, VECT_DIM);
		}
		observations[i].tpy = findMin(distancesToCentroids, nClasses);
	}
	delete[] distancesToCentroids;
}

// 计算新的聚类点，并填入centroids（centroids必须是按centroids[].tpy升序排序的）
// nDim是centroids的维度（计算距离要使用）
// 返回新centroids与原来传入的centroids距离差的最大值
double calcNewCentroids(Vect *observations, size_t nSetSize, Vect *centroids, size_t nClasses, size_t nDim) {
	printf("--> calcNewCentroids\n");
	double *sum = new double[nDim];
	Vect *oldCentroids = new Vect[nClasses];
	memcpy(oldCentroids, centroids, nClasses * sizeof(Vect));
	int numThisTyp = 0;

	for (size_t thisTpy = 0; thisTpy < nClasses; thisTpy++) {
		for (size_t i = 0; i < nDim; i++) {
			sum[i] = 0.0;
		}
		numThisTyp = 0;    // 统计当前属于thisTpy类型的样本个数

		// 遍历observations，找出类型是thisTpy的
		for (size_t i = 0; i < nSetSize; i++) {
			if (observations[i].tpy == thisTpy) {
				++numThisTyp;
				for (size_t thisDim = 0; thisDim < nDim; thisDim++) {
					sum[thisDim] += observations[i].arr[thisDim];
				}
			}
		}
		for (size_t thisDim = 0; thisDim < nDim; thisDim++) {
			if (numThisTyp <= 0) {
				fprintf(stderr, "  !! numThisTyp is %u\n", numThisTyp);
			}
			assert(numThisTyp > 0);
			centroids[thisTpy].arr[thisDim] = (sum[thisDim] / (double)numThisTyp);
		}
	}
	delete[] sum;

	double *deltaDistances = new double[nClasses];    // 反映新旧centroids之间的距离
	for (size_t i = 0; i < nClasses; i++) {
		deltaDistances[i] = calcDistance(oldCentroids[i].arr, centroids[i].arr, nDim);
		printf("  deltaDistances[%u]: %f\n", i, deltaDistances[i]);
	}
	printf("  oldCentroids: \n");
	printVects(oldCentroids, nClasses);
	printf("  centroids: \n");
	printVects(centroids, nClasses);
	double maxDeltaDistance = deltaDistances[findMax(deltaDistances, nClasses)];  // deltaDistances中最大元素的值
	delete[] deltaDistances;
	printf("<-- calcNewCentroids\n");
	return maxDeltaDistance;
}

/**************************************************************************/

// K-means入口函数
// observations为待分类样本集
// nSetSize为样本集中的元素个数
// nClasses为期望分类的分类总数默认值为2
// nDim为向量的维度
// epsilon为判决阈值。当两次迭代之间的聚类中心最大值小于epsilon时，迭代停止
void k_means(Vect * observations, size_t nSetSize, size_t nClasses, size_t nDim, double epsilon = 0.01) {
	printf("--> k_means\n");
	if (nClasses < 2) {
		printf("nClasses < 2 \? I will change it to 2.\n");
		nClasses = 2;
	}

	Vect *centroids = new Vect[nClasses];
	// 随机选取observations中的nClasses个元素作为初始的聚类中心
/* 问题往往就是出在选取初始聚类点上。如果某两个初始聚类点靠的太近，那么极易造成错误的分类*/
	size_t *candidateIndex = new size_t[nClasses];
	candidateIndex[0] = (size_t)(rand() / (double)RAND_MAX *(nClasses - 0.0) + 0.0);
	for (size_t i = 1; i < nClasses; i++) {
		size_t nRand;
		size_t *isDup;
		do {
			nRand = (size_t)(rand() / (double)RAND_MAX *(double)(nSetSize - 0.0) + 0.0);
			isDup = (size_t *)bsearch(&nRand, candidateIndex, i, sizeof(size_t),(int (*)(const void*, const void*)) cmpInt);
		} while (isDup);
		candidateIndex[i] = nRand;
	}
	printf("  candidateIndex[] is \n");
	for (size_t i = 0; i < nClasses; i++) {
		printf("[%u]\t%u\n", i, candidateIndex[i]);
	}
	for (size_t thisType = 0; thisType < nClasses; thisType++) {
		observations[candidateIndex[thisType]].tpy = thisType;
		memcpy(centroids + thisType, observations + candidateIndex[thisType], sizeof(Vect));
	}
	delete[] candidateIndex;
	printf("  Selected centroids are:\n");
	printVects(centroids, nClasses);

	markType(observations, nSetSize, centroids, nClasses);
	printf("  Initial classification of observations:\n");
	printVects(observations, nSetSize);
	// 进行迭代
	double maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
	printf("  Initial maxDeltaDistance is %f\n", maxDeltaDistance);
	while (maxDeltaDistance > epsilon) {
		size_t cnt = 0;
		maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
		printf("  maxDeltaDistance is %f\n", maxDeltaDistance);
		markType(observations, nSetSize, centroids, nClasses);
		printf("  Classification of observations of the %u(th) iteration is:\n", ++cnt);
		printVects(observations, nSetSize);
	}

	delete[] centroids;
	printf("<-- k_means\n");
}

// 打印vects数组的全部信息
void printVects(Vect *vects, size_t nSetSize) {
	for (size_t i = 0; i < nSetSize; i++) {
		printf("[%u]\t[", i);
		for (size_t j = 0; j < VECT_DIM; j++) {
			printf("%f", vects[i].arr[j]);
			if (j < VECT_DIM - 1) {
				printf(", ");
			}
		}
		printf("]\t%u\n", vects[i].tpy);
	}
}