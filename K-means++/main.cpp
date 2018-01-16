#include "k_means_plusplus.h"

int main() {
#define N 18  // ��������С
#define K 3    // �ƻ�����
	Vect observations[N];
	Vect centroids[K];
	// ��ʼ������������
	srand(time(0));   // ����뿴���������һ�㣬��ʹ��srand(time(0));
	for (size_t i = 0; i < N / 3; i++) {
		observations[i].arr[0] = static_cast<double> (5.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
		observations[i].arr[1] = static_cast<double> (5.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
	}
	for (size_t i = N / 3; i < N / 3 * 2; i++) {
		observations[i].arr[0] = static_cast<double> (10.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
		observations[i].arr[1] = static_cast<double> (10.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
	}
	for (size_t i = N / 3 * 2; i < N; i++) {
		observations[i].arr[0] = static_cast<double> (12.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
		observations[i].arr[1] = static_cast<double> (12.0 + (double)rand() / RAND_MAX *(1.0 - 0.1) + 0.1);
	}

	//for (size_t i = 0; i < N; i++) {
	//	observations[i].arr[0] = static_cast<double> ((double)rand() / RAND_MAX *(1.0 - 0.0) + 0.0);
	//	observations[i].arr[1] = static_cast<double> ((double)rand() / RAND_MAX *(1.0 - 0.0) + 0.0);
	//}

	printf("observations = \n");
	for (size_t i = 0; i < N; i++) {
		printf("[%u]\t[%f, %f]\n", i, observations[i].arr[0], observations[i].arr[1]);
	}
	//printVects(observations, N);

	//k_means_plusplus(observations, N, K, 2, 0.01);
	k_means_plusplus(observations, N, centroids, K, VECT_DIM, 0.01);
	printf("\nobservations types are (see the last column)\n");
	//for (size_t i = 0; i < N; i++) {
	//	printf("[%u]\t%u\n", i, observations[i].tpy);
	//}
	printVects(observations, N);
	return 0;
}