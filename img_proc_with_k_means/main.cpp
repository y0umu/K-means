#include <opencv2/opencv.hpp>
#include <iostream>
#include "k_means_plusplus.h"

using namespace cv;

#define PROGRAM_NAME "img_proc_with_k_means"

static void help() {
	fprintf(stderr, "Usage:\n%s <input image file> [classes]\n", PROGRAM_NAME);
}

size_t parseCmd(int argc, char**argv) {
	if ((argc != 2) && (argc != 3)) {
		help();
		exit(-1);
	}
	if (argc == 3) {
		return atoi(argv[2]);
	}
	else {
		return 2;
	}
}

int main(int argc, char** argv) {
	size_t nClasses = parseCmd(argc, argv);
	if (nClasses < 2) {
		nClasses = 2;
	}
	Mat matImg = imread(argv[1], IMREAD_COLOR);
	if (matImg.empty()) {
		fprintf(stderr, "Bad read. Does the image really exist or is it really an image?\n");
		return 1;
	}
	if (matImg.isContinuous()) {
		size_t nImgPixels = matImg.rows * matImg.cols;
		std::cout << "matImg is continuous\n";
		std::cout << "There are " << matImg.rows << " rows, "
			<< matImg.cols << "cols, "
			<< nImgPixels << " Pixels\n";
		imshow("Original image", matImg);

		std::cout << "It has type of " << matImg.type() << std::endl;
		std::cout << "CV_32F is " << CV_32F << std::endl;
// My implementation of K-means++
		Vec3b *imgPixels = matImg.ptr<Vec3b>(0);
		
		Vect *vImg = new Vect[nImgPixels];
		Vect *vCentroids = new Vect[nClasses];
		for (size_t i = 0; i < nImgPixels; i++) {
			vImg[i].arr[0] = imgPixels[i][0];
			vImg[i].arr[1] = imgPixels[i][1];
			vImg[i].arr[2] = imgPixels[i][2];
		}
		k_means_plusplus(vImg, nImgPixels, vCentroids, nClasses, 3, 0.01);
		Mat matSegImg(matImg.rows, matImg.cols, matImg.type());
		Vec3b *ptrSegImgPixels = matSegImg.ptr<Vec3b>(0);
		for (size_t i = 0; i < nImgPixels; i++) {
			ptrSegImgPixels[i][0] = vCentroids[vImg[i].tpy].arr[0];
			ptrSegImgPixels[i][1] = vCentroids[vImg[i].tpy].arr[1];
			ptrSegImgPixels[i][2] = vCentroids[vImg[i].tpy].arr[2];
		}
		delete[] vImg;
		delete[] vCentroids;

		imwrite("my_kmeans.jpg", matSegImg);
		imshow("Segmented Image", matSegImg);

// OpenCV implemetation of K-means++
		matImg.convertTo(matImg, CV_32F);
		Mat oneRow = matImg.reshape(0, 1);
		Mat matBestLabels;
		Mat matCentroids;
		TermCriteria k_means_criterial(TermCriteria::EPS, 10, 0.1);
		kmeans(oneRow,
		   	   nClasses, 
			   matBestLabels, 
			   k_means_criterial, 
			   1, 
			   KMEANS_PP_CENTERS, 
			   matCentroids);
		//std::cout << matBestLabels << std::endl;
		std::cout << "centroids selected by opencv are:\n";
		std::cout << matCentroids << std::endl;
		for (size_t thisPixel = 0; thisPixel < nImgPixels; thisPixel++) {
			Vec3f *ptr = oneRow.ptr<Vec3f>(0);
			ptr[thisPixel] = *(matCentroids.ptr<Vec3f>( matBestLabels.at<int>(thisPixel)));
		}
		Mat matSegImg2 = oneRow.reshape(0, matImg.rows);
		matSegImg2.convertTo(matSegImg2, CV_8U);
		//std::cout << matSegImg2 << std::endl;
		imwrite("opencv_kmeans.jpg", matSegImg2);
		imshow("Segmented using OpenCV", matSegImg2);
		waitKey(0);
	}
	else {
		printf("Sorry, matImg is NOT continuous\n");
	}

	return 0;
}