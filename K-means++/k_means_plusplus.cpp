#include "k_means_plusplus.h"


// ������㺯��������ʸ��a��ʸ��b���������
double calcDistance(double *a, double *b, size_t dim) {
	double sum = 0.0;
	double diff;
	for (size_t i = 0; i < dim; i++) {
		diff = a[i] - b[i];
		sum += (diff * diff);
	}
	return sqrt(sum);
}

// �ҵ�����a����Сֵ����������a�е��±�
// ��������Ϊdouble, int, unsigned��֧�ֱȽ�������Ļ���
template <typename T>
size_t findMin(const T *a, size_t nLen) {
	size_t indexOfMin = 0;  // ��Сֵ���ڵ��±�
	for (size_t i = 1; i < nLen; i++) {
		if (a[i] < a[indexOfMin]) {
			indexOfMin = i;
		}
	}
	return indexOfMin;
}


//// ���������С���ң�����������СԪ�ص��±�
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

// �ҵ�����a�����ֵ����������a�е��±�
// ��������Ϊdouble, int, unsigned��֧�ֱȽ�������Ļ���
template <typename T>
size_t findMax(const T *a, size_t nLen) {
	size_t indexOfMax = 0;  // ��Сֵ���ڵ��±�
	for (size_t i = 1; i < nLen; i++) {
		if (a[i] > a[indexOfMax]) {
			indexOfMax = i;
		}
	}
	return indexOfMax;
}

// bsearch�Ļص���������markType���õ���
int cmpVects_arr(Vect *key, Vect *datum) {
	int diff = memcmp(key->arr, datum->arr, VECT_DIM * sizeof(double));
	return diff;
}

// bsearch�Ļص���������k_means���õ���
int cmpInt(int *key, int *datum) {
	int diff = *key - *datum;
	return diff;
}

// qsort�Ļص�����
int cmpVects_typ(Vect *elem1, Vect *elem2) {
	int diff = (int)(elem1->tpy) - (int)(elem2->tpy);
	return diff;
}

// K-means++���С���������ĵ�һ����ʼ��������ȫ��nClasses����������
// observationsΪ������������
// nSetSizeΪ�������е�Ԫ�ظ���
// initialCentoridIndexΪ��ʼ����������observations�е��±�
// centroidsIndexΪ��������������ɵ�����nClasses������������observations�е��±�
// nClassesΪ��������ķ�������
// nDimΪ������ά��
//! �˺���ûд��
void genCentroids(Vect *observations, size_t nSetSize, 
	              size_t initialCentoridIndex, size_t *centroidsIndex, size_t nClasses,
                  size_t nDim) {
	//printf("--> genCentroids\n");
	centroidsIndex[0] = initialCentoridIndex;
	double *distanceCandidates = new double[nClasses];  // ���浱ǰ���ݵ㵽���������ľ��롣��findMin�����ҵ���Сֵ���ڵ��±�
	double *distanceToNearestCentroid = new double[nSetSize]; // ���浱ǰ���ݵ㵽��������ľ���
	size_t nCentroids = 1;  // ��ǰ���ҵ��ľ������Ŀ
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

// Ϊobservations��Ԫ�ر�����ͣ����޸���typ����
// observationsΪ������������
// nSetSizeΪ�������е�Ԫ�ظ���
// nClassesΪ��������ķ�������
// centroidsInd�Ǳ�ʾ��ǰ�������ĵ�һά���飬���СΪnClasses��centroidsInd��Ԫ��Ϊ��������������observations�е������±�
void markType(Vect *observations, size_t nSetSize, Vect *centroids, size_t nClasses) {
	// ��֤centroids�ǰ���centroids[].tpy���������
	qsort(centroids, nClasses, sizeof(Vect), (int (*)(const void*, const void*))cmpVects_typ);

	double *distancesToCentroids = new double[nClasses];
	for (size_t i = 0; i < nSetSize; i++) {
		//Vect *isDup = (Vect *)bsearch(observations + i, centroids, nClasses, sizeof(Vect), (int(*)(const void*, const void*))cmpVects_arr);
		////int diff = memcmp(observations + i, centroids + i, sizeof(Vect));
		//if (isDup) {  // ��ǰ��observation��centroids��ͬһ����
		//	continue;
		//}
		for (size_t j = 0; j < nClasses; j++) {
			distancesToCentroids[j] = calcDistance(observations[i].arr, centroids[j].arr, VECT_DIM);
		}
		observations[i].tpy = findMin(distancesToCentroids, nClasses);
	}
	delete[] distancesToCentroids;
}

// �����µľ���㣬������centroids��centroids�����ǰ�centroids[].tpy��������ģ�
// nDim��centroids��ά�ȣ��������Ҫʹ�ã�
// ������centroids��ԭ�������centroids���������ֵ
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
		numThisTyp = 0;    // ͳ�Ƶ�ǰ����thisTpy���͵���������

		// ����observations���ҳ�������thisTpy��
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

	double *deltaDistances = new double[nClasses];    // ��ӳ�¾�centroids֮��ľ���
	for (size_t i = 0; i < nClasses; i++) {
		deltaDistances[i] = calcDistance(oldCentroids[i].arr, centroids[i].arr, nDim);
		//printf("  deltaDistances[%u]: %f\n", i, deltaDistances[i]);
	}
	//printf("  oldCentroids: \n");
	//printVects(oldCentroids, nClasses);
	//printf("  centroids: \n");
	//printVects(centroids, nClasses);
	double maxDeltaDistance = deltaDistances[findMax(deltaDistances, nClasses)];  // deltaDistances�����Ԫ�ص�ֵ
	delete[] deltaDistances;
	//printf("<-- calcNewCentroids\n");
	return maxDeltaDistance;
}

/**************************************************************************/

// K-means++��ں��� (1/2)
// observationsΪ������������
// nSetSizeΪ�������е�Ԫ�ظ���
// nClassesΪ��������ķ�������Ĭ��ֵΪ2
// nDimΪ������ά��
// epsilonΪ�о���ֵ�������ε���֮��ľ����������ֵС��epsilonʱ������ֹͣ
void k_means_plusplus(Vect * observations, size_t nSetSize, size_t nClasses, size_t nDim, double epsilon) {
	//printf("--> k_means\n");
	if (nClasses < 2) {
		printf("nClasses < 2 \? I will change it to 2.\n");
		nClasses = 2;
	}

	Vect *centroids = new Vect[nClasses];

	// ����K-means++�ķ������ɳ�ʼ��������
	size_t *centroidsIndex = new size_t[nClasses];
	centroidsIndex[0] = (size_t)(rand() / (double)RAND_MAX *(nClasses - 0.0) + 0.0);
	genCentroids(observations, nSetSize, centroidsIndex[0], centroidsIndex, nClasses, nDim);
	printf("  centroidsIndex[] is \n");
	for (size_t i = 0; i < nClasses; i++) {
		printf("[%u]\t%u\n", i, centroidsIndex[i]);
	}
	// ��0~(nClasses-1)�������
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
	// ���е���
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

// K-means++��ں��� (2/2)
// observationsΪ������������
// nSetSizeΪ�������е�Ԫ�ظ���
// centroidsΪ�ش��ľ�������
// nClassesΪ��������ķ�������Ĭ��ֵΪ2
// nDimΪ������ά��
// epsilonΪ�о���ֵ�������ε���֮��ľ����������ֵС��epsilonʱ������ֹͣ
void k_means_plusplus(Vect * observations, size_t nSetSize, 
	                  Vect * centroids,size_t nClasses, 
	                  size_t nDim, double epsilon) {
	//printf("--> k_means\n");
	if (nClasses < 2) {
		printf("nClasses < 2 \? I will change it to 2.\n");
		nClasses = 2;
	}

	// ����K-means++�ķ������ɳ�ʼ��������
	size_t *centroidsIndex = new size_t[nClasses];
	centroidsIndex[0] = (size_t)(rand() / (double)RAND_MAX *(nClasses - 0.0) + 0.0);
	genCentroids(observations, nSetSize, centroidsIndex[0], centroidsIndex, nClasses, nDim);
	printf("  centroidsIndex[] is \n");
	for (size_t i = 0; i < nClasses; i++) {
		printf("[%u]\t%u\n", i, centroidsIndex[i]);
	}
	// ��0~(nClasses-1)�������
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
	// ���е���
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

// ��ӡvects�����ȫ����Ϣ
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