#include "k_means_plusplus.h"


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


//// 泛型最大最小查找，返回最大或最小元素的下标
//size_t findMaxOrMin(void* base, size_t numOfElements, size_t sizeOfElements, 
//                  int (*comparator)(const void*, const void*) ) {
//	size_t ind = 0;
//	for (size_t i = 1; i < numOfElements; i++) {
//		if (comparator( (uchar*)base + (i * sizeOfElements),
//			            (uchar*)base + (ind * sizeOfElements)) ) {
//			ind = i;
//		}
//	}
//	return ind;
//}

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

// K-means++独有。根据输入的第一个初始中心生成全部nClasses个聚类中心
// observations为待分类样本集
// nSetSize为样本集中的元素个数
// initialCentoridIndex为初始聚类中心在observations中的下标
// centroidsIndex为输出，保存所生成的所有nClasses个聚类中心在observations中的下标
// nClasses为期望分类的分类总数
// nDim为向量的维度
//! 此函数没写完
void genCentroids(Vect *observations, size_t nSetSize, 
	              size_t initialCentoridIndex, size_t *centroidsIndex, size_t nClasses,
                  size_t nDim) {
	//printf("--> genCentroids\n");
	centroidsIndex[0] = initialCentoridIndex;
	double *distanceCandidates = new double[nClasses];  // 保存当前数据点到各个聚类点的距离。由findMin函数找到最小值所在的下标
	double *distanceToNearestCentroid = new double[nSetSize]; // 保存当前数据点到最近聚类点的距离
	size_t nCentroids = 1;  // 当前已找到的聚类点数目
	while (nCentroids < nClasses) {
		for (size_t thisObsInd = 0; thisObsInd < nSetSize; thisObsInd++) {
			for (size_t j = 0; j < nCentroids; j++) {
				distanceCandidates[j] = calcDistance(observations[thisObsInd].arr,
					                                 observations[centroidsIndex[j]].arr,
					                                 nDim);
			}
			distanceToNearestCentroid[thisObsInd] = distanceCandidates[findMin(distanceCandidates, nCentroids)];
		}
		centroidsIndex[nCentroids] = findMax(distanceToNearestCentroid, nSetSize);
		++nCentroids;
	}
	delete[] distanceCandidates;
	delete[] distanceToNearestCentroid;
	//printf("<-- genCentroids\n");
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
	//printf("--> calcNewCentroids\n");
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
		//printf("  deltaDistances[%u]: %f\n", i, deltaDistances[i]);
	}
	//printf("  oldCentroids: \n");
	//printVects(oldCentroids, nClasses);
	//printf("  centroids: \n");
	//printVects(centroids, nClasses);
	double maxDeltaDistance = deltaDistances[findMax(deltaDistances, nClasses)];  // deltaDistances中最大元素的值
	delete[] deltaDistances;
	//printf("<-- calcNewCentroids\n");
	return maxDeltaDistance;
}

/**************************************************************************/

// K-means++入口函数 (1/2)
// observations为待分类样本集
// nSetSize为样本集中的元素个数
// nClasses为期望分类的分类总数默认值为2
// nDim为向量的维度
// epsilon为判决阈值。当两次迭代之间的聚类中心最大值小于epsilon时，迭代停止
void k_means_plusplus(Vect * observations, size_t nSetSize, size_t nClasses, size_t nDim, double epsilon) {
	//printf("--> k_means\n");
	if (nClasses < 2) {
		printf("nClasses < 2 \? I will change it to 2.\n");
		nClasses = 2;
	}

	Vect *centroids = new Vect[nClasses];

	// 利用K-means++的方法生成初始聚类中心
	size_t *centroidsIndex = new size_t[nClasses];
	centroidsIndex[0] = (size_t)(rand() / (double)RAND_MAX *(nClasses - 0.0) + 0.0);
	genCentroids(observations, nSetSize, centroidsIndex[0], centroidsIndex, nClasses, nDim);
	printf("  centroidsIndex[] is \n");
	for (size_t i = 0; i < nClasses; i++) {
		printf("[%u]\t%u\n", i, centroidsIndex[i]);
	}
	// 从0~(nClasses-1)标记类型
	for (size_t thisType = 0; thisType < nClasses; thisType++) {
		observations[centroidsIndex[thisType]].tpy = thisType;
		memcpy(centroids + thisType, observations + centroidsIndex[thisType], sizeof(Vect));
	}
	delete[] centroidsIndex;

	printf("  Selected centroids are:\n");
	printVects(centroids, nClasses);
	markType(observations, nSetSize, centroids, nClasses);
	//printf("  Initial classification of observations:\n");
	//printVects(observations, nSetSize);
	// 进行迭代
	double maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
	//printf("  Initial maxDeltaDistance is %f\n", maxDeltaDistance);
	while (maxDeltaDistance > epsilon) {
		//size_t cnt = 0;
		maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
		//printf("  maxDeltaDistance is %f\n", maxDeltaDistance);
		markType(observations, nSetSize, centroids, nClasses);
		//printf("  Classification of observations of the %u(th) iteration is:\n", ++cnt);
		//printVects(observations, nSetSize);
	}

	delete[] centroids;
	//printf("<-- k_means\n");
}

// K-means++入口函数 (2/2)
// observations为待分类样本集
// nSetSize为样本集中的元素个数
// centroids为回传的聚类中心
// nClasses为期望分类的分类总数默认值为2
// nDim为向量的维度
// epsilon为判决阈值。当两次迭代之间的聚类中心最大值小于epsilon时，迭代停止
void k_means_plusplus(Vect * observations, size_t nSetSize, 
	                  Vect * centroids,size_t nClasses, 
	                  size_t nDim, double epsilon) {
	//printf("--> k_means\n");
	if (nClasses < 2) {
		printf("nClasses < 2 \? I will change it to 2.\n");
		nClasses = 2;
	}

	// 利用K-means++的方法生成初始聚类中心
	size_t *centroidsIndex = new size_t[nClasses];
	centroidsIndex[0] = (size_t)(rand() / (double)RAND_MAX *(nClasses - 0.0) + 0.0);
	genCentroids(observations, nSetSize, centroidsIndex[0], centroidsIndex, nClasses, nDim);
	printf("  centroidsIndex[] is \n");
	for (size_t i = 0; i < nClasses; i++) {
		printf("[%u]\t%u\n", i, centroidsIndex[i]);
	}
	// 从0~(nClasses-1)标记类型
	for (size_t thisType = 0; thisType < nClasses; thisType++) {
		observations[centroidsIndex[thisType]].tpy = thisType;
		memcpy(centroids + thisType, observations + centroidsIndex[thisType], sizeof(Vect));
	}
	delete[] centroidsIndex;

	printf("  Selected centroids are:\n");
	printVects(centroids, nClasses);
	markType(observations, nSetSize, centroids, nClasses);
	//printf("  Initial classification of observations:\n");
	//printVects(observations, nSetSize);
	// 进行迭代
	double maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
	//printf("  Initial maxDeltaDistance is %f\n", maxDeltaDistance);
	while (maxDeltaDistance > epsilon) {
		size_t cnt = 0;
		maxDeltaDistance = calcNewCentroids(observations, nSetSize, centroids, nClasses, nDim);
		//printf("  maxDeltaDistance is %f\n", maxDeltaDistance);
		markType(observations, nSetSize, centroids, nClasses);
		//printf("  Classification of observations of the %u(th) iteration is:\n", ++cnt);
		//printVects(observations, nSetSize);
	}

	//printf("<-- k_means\n");
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