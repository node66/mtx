#include "Matrix.h"
#include<iostream>
#include <iomanip>
#include<time.h>

namespace algebra {

	using namespace std;


	Matrix::Matrix() :rows_(1), cols_(1) {
		allocSpace();
		for (int i = 0; i < rows_; i++) {
			for (int j = 0; j < cols_; j++) {
				data[i*cols_ + j] = 0;
			}
		}
	}

	Matrix::Matrix(const Matrix& m) : rows_(m.rows_), cols_(m.cols_) {
		allocSpace();
		for (int i = 0; i < rows_; i++) {
			for (int j = 0; j < cols_; j++) {
				data[i*cols_ + j] = m.data[i*cols_ + j];
			}
		}
	}

	Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
		allocSpace();
		for (int i = 0; i < rows_; i++) {
			for (int j = 0; j < cols_; j++) {
				data[i*cols_ + j] = 0;
			}
		}
	}

	Matrix::~Matrix()
	{
		delete data;
	}

	double & Matrix::operator()(int i, int j) {
		return data[i*(cols_)+j];
	}


	Matrix & Matrix::operator=(const Matrix& m) {
		if (this == &m) {
			return *this;
		}

		if (rows_ != m.rows_ || cols_ != m.cols_) {
			delete[] data;

			rows_ = m.rows_;
			cols_ = m.cols_;
			allocSpace();
		}

		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data[i*cols_ + j] = m.data[i*cols_ + j];
			}
		}
		return *this;
	}

	Matrix & Matrix::operator=(const double value)
	{
		for (int i = 0; i < rows_; i++) {
			for (int j = 0; j < cols_; j++) {
				data[i*cols_ + j] = value;
			}
		}
		return *this;
	}

	Matrix& Matrix::operator+=(const Matrix& m) {
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data[i*cols_ + j] += m.data[i*cols_ + j];
			}
		}
		return *this;
	}

	Matrix& Matrix::operator-=(const Matrix& m) {
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data[i*cols_ + j] -= m.data[i*cols_ + j];
			}
		}
		return *this;
	}

	Matrix & Matrix::operator-=(double value)
	{
		Matrix temp(rows_, cols_);
		for (int i = 0; i < rows_; ++i)
			for (int j = 0; j < cols_; ++j)
				data[i*cols_ + j] -= value;
		return *this;
	}

	Matrix& Matrix::operator*=(const Matrix& m) {
		Matrix temp(rows_, m.cols_);
		for (int i = 0; i < temp.rows_; ++i) {
			for (int j = 0; j < temp.cols_; ++j) {
				for (int k = 0; k < cols_; ++k) {
					temp.data[i*cols_ + j] += (data[i*cols_ + k] * m.data[k*cols_ + j]);
				}
			}
		}
		return (*this = temp);
	}

	Matrix& Matrix::operator*=(double value) {
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data[i*cols_ + j] *= value;
			}
		}
		return *this;
	}

	Matrix& Matrix::operator/=(double value) {
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				data[i*cols_ + j] /= value;
			}
		}
		return *this;
	}

	Matrix Matrix::operator-(const Matrix & right_arg)
	{
		Matrix temp(right_arg.rows_, right_arg.cols_);
		for (int j = 0; j < right_arg.rows_; j++)
			for (int i = 0; i < right_arg.cols_; i++)
				temp.data[j * right_arg.cols_ + i] = data[j * right_arg.cols_ + i] - right_arg.data[j * right_arg.cols_ + i];
		return temp;
	}

	Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
		Matrix temp(lhs);
		return (temp += rhs);
	}

	Matrix operator-(const Matrix& lhs, const Matrix& rhs) {
		Matrix temp(lhs);
		return (temp -= rhs);
	}

	Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
		Matrix temp(lhs);
		return (temp *= rhs);
	}

	Matrix operator*(const Matrix& lhs, double value) {
		Matrix temp(lhs);
		return (temp *= value);
	}

	Matrix operator-(const Matrix &lhs, double value)
	{
		Matrix temp(lhs);
		return (temp -= value);
	}

	Matrix operator*(double value, const Matrix& m) {
		return (m * value);
	}

	Matrix operator/(const Matrix& m, double value) {
		Matrix temp(m);
		return (temp /= value);
	}

	std::ostream & operator<<(std::ostream & os, const Matrix & m)
	{
		for (int i = 0; i < m.rows_; i++)
		{
			for (int j = 0; j < m.cols_; j++)
			{
				os << setw(8) << setprecision(3) << setiosflags(ios::fixed | ios::showpoint) << m.data[i*m.cols_ + j] << " ";
			}
			os << "\n";
		}
		return os;
	}


	std::istream & operator>>(std::istream &is, Matrix &m)
	{
		for (int i = 0; i < m.rows_; i++)
		{
			for (int j = 0; j < m.cols_; j++)
			{
				is >> m.data[i, j];
			}
		}
		return is;
	}

	void Matrix::fillrand() {
		srand(time(NULL));
		for (int i = 0; i < rows_; i++)
			for (int j = 0; j < cols_; j++)
				data[i*cols_ + j] = rand() % 6 - 5;
	}
	/*!
		LU разложение с дальнейшим решением СЛАУ
	*/
	void Matrix::LU(Matrix& L, Matrix& U, Vector& b)
	{
		int n = cols_;
		U = *this;

		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++)
				L(j, i) = U(j, i) / U(i, i);

		for (int k = 1; k < n; k++)
		{
			for (int i = k - 1; i < n; i++)
				for (int j = i; j < n; j++)
					L(j, i) = U(j, i) / U(i, i);

			for (int i = k; i < n; i++)
				for (int j = k - 1; j < n; j++)
					U(i, j) = U(i, j) - L(i, k - 1) * U(k - 1, j);
		}
		cout << "Матрица А " << endl;
		cout << *this << endl;
		cout << "Матрица L " << endl;
		cout << L << endl;
		cout << "Матрица U " << endl;
		cout << U << endl;
		cout << "Матрица L * U " << endl;
		cout << L * U << endl;

		/*Решение СЛАУ */
		/*     L*y=b    */
		/*  Прямая подстановка  */
		Vector y(n);
		y(0) = b(0);
		for (int i = 1; i < n; i++) {
			double s = b(i);
			for (int j = 0; j < i + 1; j++) {
				s -= L(i, j)*y(j);
			}
			y(i) = s;
		}
		/*    U*x=y   */
		/*   Обратная подстановка  */
		Vector x(n);
		x(n - 1) = y(n - 1) / U(n - 1, n - 1);
		for (int i = n - 2; i >= 0; i--) {
			double s = y(i);
			for (int j = n - 1; j >= i; j--) {
				s -= U(i, j)*x(j);
			}
			x(i) = s / U(i, i);
		}
		cout << " Решение СЛАУ " << endl;
		cout << " L*y=b    Прямая подстановка " << endl;
		cout << " y = " << y << endl;
		cout << " U*x=y   Обратная подстановка " << endl;
		cout << " x = " << x << endl;
	}

	/*!
		Просто LU разложение
	*/
	void Matrix::LU(Matrix& L, Matrix& U)
	{
		int n = cols_;
		U = *this;

		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++)
				L(j, i) = U(j, i) / U(i, i);

		for (int k = 1; k < n; k++)
		{
			for (int i = k - 1; i < n; i++)
				for (int j = i; j < n; j++)
					L(j, i) = U(j, i) / U(i, i);

			for (int i = k; i < n; i++)
				for (int j = k - 1; j < n; j++)
					U(i, j) = U(i, j) - L(i, k - 1) * U(k - 1, j);
		}

	}

	/*!
	LU обращение матрицы
	*/
	Matrix Matrix::LU_Inversion()
	{
		unsigned int n = this->cols_;
		Matrix E(n, n), Answer(n, n), L(n, n), U(n, n);
		Vector e(n), y(n), x(n);
		this->LU(L, U);

		for (unsigned int i = 0; i < n; i++)
			for (unsigned int j = 0; j < n; j++)
				i == j ? E(i, j) = 1 : E(i, j) = 0;

		for (unsigned int i = 0; i < n; i++)
		{
			for (unsigned int j = 0; j < n; j++)
				e(j) = E(i, j);

			y.fill(0);
			y(0) = e(0);
			for (unsigned int i = 1; i < n; i++)
			{
				double s = e(i);
				for (unsigned int j = 0; j < i + 1; j++)
					s -= L(i, j)*y(j);

				y(i) = s;
			}
			x.fill(0);
			x(n - 1) = y(n - 1) / U(n - 1, n - 1);
			for (int i = n - 2; i >= 0; i--)
			{
				double s = y(i);
				for (int j = n - 1; j >= i; j--)
					s -= U(i, j)*x(j);

				x(i) = s / U(i, i);
			}
			for (unsigned int j = 0; j < n; j++)
				Answer(j, i) = x(j);
		}
		return Answer;
	}


	Matrix Matrix::LDL_T_Inversion()
	{
		unsigned int n = this->cols_;
		Matrix E(n, n), Answer(n, n), L(n, n), D(n, n), U(n, n);
		Vector e(n), y(n), x(n);
		this->LDL_T(L, D);
		U = D * L.transpose();

		for (unsigned int i = 0; i < n; i++)
			E(i, i) = 1;

		for (unsigned int i = 0; i < n; i++)
		{
			for (unsigned int j = 0; j < n; j++)
				e(j) = E(i, j);

			y.fill(0);
			y(0) = e(0);
			for (unsigned int i = 1; i < n; i++)
			{
				double s = e(i);
				for (unsigned int j = 0; j < i + 1; j++)
					s -= L(i, j)*y(j);

				y(i) = s;
			}
			x.fill(0);
			x(n - 1) = y(n - 1) / U(n - 1, n - 1);
			for (int i = n - 2; i >= 0; i--)
			{
				double s = y(i);
				for (int j = n - 1; j >= i; j--)
					s -= U(i, j)*x(j);

				x(i) = s / U(i, i);
			}
			for (unsigned int j = 0; j < n; j++)
				Answer(j, i) = x(j);
		}

		return Answer;
	}
	void Matrix::QR(Matrix & Q, Matrix & R)
	{
		unsigned int length = this->cols_;
		Matrix A = *this;
		Vector q(length), a(length);

		for (unsigned int z = 0; z < length; z++)
		{
			for (unsigned int j = 0; j < length; j++)
			{
				a(j) = A(j, z);
				q(j) = a(j);
			}
			for (unsigned int j = 0; j < length; j++)
			{
				if (j > 1)
				{
					for (unsigned int i = 0; i < j - 1; i++)
					{
						R(i, j) = q(i) *  a(j);
						q(j) = q(j) - R(i, j)*q(j);
					}
				}
				R(j, j) = q.norm();

				if (R(j, j) == 0)
					cerr << "ai линейно зависит от a1,..., ai−1";

				q(j) = q(j) / R(j, j);
			}
			cout << q;
			for (unsigned int j = 0; j < length; j++) 
				Q(z, j) = q(z);
		}
		cout << *this << endl << R << endl << Q << endl << R * Q;
	}

	void Matrix::LDL_T(Matrix &L, Matrix &D)
	{
		int n = this->cols_;
		Vector d(n);
		double sum = 0;
		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++)
			{
				sum = data[j*n + i];
				for (int k = 0; k < i; k++)
					sum -= L(i, k) *d(k) * L(j, k);
				if (i == j)
				{
					if (sum <= 0) new exception;
					{
						d(i) = sum;
						L(i, i) = 1;
					}
				}
				else L(j, i) = sum / d(i);
			}
		for (int i = 0; i < n; i++)
			D(i, i) = d(i);
	}
	

	void Matrix::LDL_T(Matrix &L, Matrix &D, Vector &b)
	{
		int n = cols_;
		Vector d(n);	
		double sum = 0;
		for (int i = 0; i < n; i++)
			for (int j = i; j < n; j++)
			{
				sum = data[j*n + i];
				for (int k = 0; k < i; k++)
					sum -= L(i, k) *d(k) * L(j, k);
				if (i == j)
				{
					if (sum <= 0) new exception;
					{
						d(i) = sum;
						L(i, i) = 1;
					}
				}
				else L(j, i) = sum / d(i);
			}
		for (int i = 0; i < n; i++)
			D(i, i) = d(i);

		cout << "Матрица А " << endl;
		cout << *this << endl;
		cout << "Матрица L " << endl;
		cout << L << endl;
		cout << "Матрица D " << endl;
		cout << D << endl;
		cout << "Матрица L^t " << endl;
		cout << L.transpose() << endl;
		cout << "Матрица L*D*L^t " << endl;
		cout << L * D*L.transpose() << endl;
	
		/*Решение СЛАУ */
		/*  Прямая подстановка  */
		Vector y(n);
		y(0) = b(0);
		for (int i = 1; i < n; i++) {
			double s = b(i);
			for (int j = 0; j < i + 1; j++) {
				s -= L(i, j)*y(j);
			}
			y(i) = s;
		}
		/*   Обратная подстановка  */
		Vector x(n);
		Matrix U(n, n);
		U = D * L.transpose();
		x(n - 1) = y(n - 1) / U(n - 1, n - 1);
		for (int i = n - 2; i >= 0; i--) {
			double s = y(i);
			for (int j = n - 1; j >= i; j--) {
				s -= U(i, j)*x(j);
			}
			x(i) = s / U(i, i);
		}
		cout << " Решение СЛАУ " << endl;
		cout << " L*y=b    Прямая подстановка " << endl;
		cout << " y = " << y << endl;
		cout << " D*L^T*x=y   Обратная подстановка " << endl;
		cout << " x = " << x << endl;
	}

	Matrix Matrix::transpose()
	{
		Matrix temp(cols_, rows_);
		for (int i = 0; i < rows_; ++i) {
			for (int j = 0; j < cols_; ++j) {
				temp.data[j*cols_+i] = data[i*cols_+j];
			}
		}
		return temp;
	}

	void Matrix::allocSpace() {
		data = new double[rows_ * cols_];
	}
}