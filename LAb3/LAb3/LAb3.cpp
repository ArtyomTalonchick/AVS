#include "pch.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>

using namespace std;

void multiplyVectors() {
	const int n = 1e+6;
	vector<int> vector1(n, 1), vector2(n, 2);

	auto start = clock();
	int result = 0;
	for (int i = 0; i < vector1.size(); i++) {
		result += vector1[i] * vector2[i];
	}
	auto time = clock() - start;
	cout << "1. result:" << result << "\ntime: " << time << " ms" << endl;

	auto start2 = clock();
	result = 0;
#pragma omp parallel for shared(vector1, vector2) reduction(+:result) schedule(static)
	for (int i = 0; i < vector1.size(); i++) {
		result += vector1[i] * vector2[i];
	}
	auto time2 = clock() - start2;
	cout << "2. result:" << result << "\ntime: " << time2 << " ms" << endl;
}


void multiplyMatrices() {
	const int a = 1; //заполняем matrix1
	const int b = 2; //заполняем matrix2
	const int l = 60;
	const int m = 70;
	const int n = 80;
	vector<vector<int>> matrix1(l), matrix2(m), result(l); // (l,m), (m,n), (l,n)

	for (int i = 0; i < l; ++i) {
		matrix1[i].resize(m);
		for (int j = 0; j < m; ++j)
			matrix1[i][j] = a;
	}
	for (int i = 0; i < m; ++i) {
		matrix2[i].resize(n);
		for (int j = 0; j < n; ++j)
			matrix2[i][j] = b;
	}
	for (int i = 0; i < l; ++i) {
		result[i].resize(n);
	}

	auto start = clock();
	for (int i = 0; i < l; i++) {
		for (int j = 0; j < n; j++) {
			result[i][j] = 0;
			for (int k = 0; k < m; k++) {
				result[i][j] += matrix1[i][k] * matrix2[k][j];
			}
			//cout << result[i][j] << "\t"; //m*a*b
		}
		//cout << endl;
	}
	auto time = clock() - start;
	cout << "time1: " << time << " ms" << endl;

	auto start2 = clock();
#pragma omp parallel for shared(matrix1, matrix2, result) schedule(static)
	for (int i = 0; i < l; i++) {
		for (int j = 0; j < n; j++) {
			int temp_value = 0;
			for (int k = 0; k < m; k++) {
				temp_value += matrix1[i][k] * matrix2[k][j];
			}
			result[i][j] = temp_value;
		}
//#pragma omp critical
//		{
//			cout << result[i][1] << " ";
//		}
	}
	auto time2 = clock() - start2;
	cout << "time2: " << time2 << " ms" << endl;
	cout << (double)time / time2 << endl;
}

int main()
{
	//multiplyVectors();

	multiplyMatrices();
}
