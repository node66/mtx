#include "Vector.h"
#include "Matrix.h"
#include<iostream>
#include <iomanip>
#include <cmath>

namespace algebra
{
	using namespace std;

	Vector::Vector() : cols_(1) {
		allocSpace();
		data[0] = 0;
	}

	Vector::~Vector() {
		delete[] data;
	}

	double & Vector::operator()( int i) {
		return data[i];
	}

	double Vector::operator*(Vector & v)
	{
		double result = 0.0;
		for ( int i = 0; i < v.cols_; i++)
		{
			result += (*this)(i) * v(i);
		}
		return result;
	}

	Vector::Vector(int cols) :cols_(cols) {
		allocSpace();
		for (int i = 0; i < cols_; i++)
			data[i] = 0;
	}

	Vector::Vector(const Vector& v) : cols_(v.cols_) {
		allocSpace();
		for (int i = 0; i < cols_; i++)
			data[i] = v.data[i];
	}

	Vector & Vector::operator=(const double & val)
	{
		for (int i = 0; i < this->cols_; i++)
		{
			this->data[i] = val;
		}
		return *this;
	}

	Vector& Vector::operator=(const Vector&v) {
		if (this == &v) {
			return *this;
		}
		if (cols_ != v.cols_) {
			cols_ = v.cols_;
			allocSpace();
		}
		for (int i = 0; i < cols_; i++) {
			data[i] = v.data[i];
		}
		return *this;
	}

	Vector& Vector::operator+=(const Vector&v) {
		for (int i = 0; i < cols_; i++) {
			data[i] += v.data[i];
		}
		return *this;
	}

	Vector& Vector::operator-=(const Vector&v) {
		for (int i = 0; i < cols_; i++) {
			data[i] -= v.data[i];
		}
		return *this;
	}

	Vector& Vector::operator*=(const Vector&v) {
		Vector temp(cols_);
		for (int i = 0; i < cols_; i++) {
			temp.data[i] += data[i] * v.data[i];
		}
		return (*this = temp);
	}


	Vector& Vector::operator*=(double value) {
		for (int i = 0; i < cols_; i++) {
			data[i] *= value;
		}
		return *this;
	}


	Vector& Vector::operator/=(double value) {
		for (int i = 0; i < cols_; i++) {
			data[i] /= value;
		}
		return *this;
	}

	Vector & Vector::operator/(double value )
	{
		unsigned int n = this->cols_;
		Vector temp(n);
		temp = *this;
		for (unsigned i = 0; i < n; i++)
		{
			temp(i) /= value;
		}
		return (*this = temp);
	}

	Vector & Vector::operator*(Matrix & m)
	{
		Vector result(cols_);
		for (int i = 0; i < cols_; i++)
			for (int j = 0; j < cols_; j++)
				result(i) += data[j] * m(j, i);
		
		return (*this = result);
	}

	Vector operator*(const Vector &v, Matrix &m)
	{
		Vector result(v.cols_);
		for (int i = 0; i < v.cols_; i++)
			for (int j = 0; j < v.cols_; j++)
				result(i) += v.data[j] * m(i,j);

		return  result;
	}

	Vector operator*(Matrix &m,  Vector &v)
	{
		int n = v.cols_;
		Vector res(n);
		for (int i = 0; i < n; i++)
		{
			double temp = 0.0;
			for (int j = 0; j < n; j++)
			{
				temp += m(i,j) * v(j);
			}
			res(i) = temp;
		}
		return res;
	}

	Vector operator+(const Vector& lhs, const Vector& rhs) {
		Vector temp(lhs);
		return (temp += rhs);
	}

	Vector operator-(const Vector& lhs, const Vector& rhs) {
		Vector temp(lhs);
		return (temp -= rhs);
	}

	Vector operator*(const Vector& lhs, const Vector& rhs) {
		Vector temp(lhs);
		return (temp *= rhs);
	}

	Vector operator*(const Vector& lhs, double value) {
		Vector temp(lhs);
		return (temp *= value);
	}

	std::ostream & operator<<(std::ostream & os, const Vector & v)
	{
		for (int i = 0; i < v.cols_; i++)
		{
			os << setw(8) << setprecision(3) << setiosflags(ios::fixed | ios::showpoint) << v.data[i] << " ";
		}
		cout << endl;
		return os;
	}

	std::istream & operator>>(std::istream &is, Vector &v)
	{
		for (int i = 0; i < v.cols_; i++)
		{
			is >> v.data[i];
		}
		return is;
	}

	Vector operator*(double value, const Vector& rhs) {
		return (rhs * value);
	}

	Vector operator/(const Vector& lhs, double value) {
		Vector temp(lhs);
		return (temp /= value);
	}

	void Vector::allocSpace() {
		data = new double[cols_];
	}

	double Vector::norm()
	{
		Vector norm = * this;
		double temp;
		double s = 0.0;
		for ( int i = 0; i < this->cols_; i++)
		{
			temp = pow(abs(norm(i)), 2.0);
			s += temp;
		}
		return sqrt(s) ;
	}


	void Vector::fill(double value) {
		for (int i = 0; i < cols_; i++) {
			data[i] = value;
		}
	}

	void Vector::fillrand() {
		for (int i = 0; i < cols_; i++) {
			data[i] = rand() % 6 - 5;
		}
	}

}