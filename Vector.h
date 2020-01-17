#pragma once
#include"Matrix.h"
#include <iostream>


namespace algebra
{
	class Matrix;
	class Vector
	{
	public:
		Vector();
		~Vector();
		Vector(int);
		Vector(const Vector&);
    
		double operator*(Vector &);
	
		double & operator()( int);	
		Vector & operator=(const double&);
		Vector & operator=(const Vector&);
		Vector & operator+=(const Vector&);
		Vector & operator-=(const Vector&);
		Vector & operator*=(const Vector&);
		Vector & operator*=(double);
		Vector & operator/=(double);
		Vector & operator/(double);
		Vector & operator*(Matrix&);
	
		friend	Vector operator*(const Vector&,  Matrix&);
		friend	Vector operator*( Matrix&,  Vector&);
		friend	Vector operator+(const Vector&, const Vector&);
		friend	Vector operator-(const Vector&, const Vector&);
		friend	Vector operator*(const Vector&, const Vector&);
		friend	Vector operator*(const Vector&, double);
		friend std::ostream& operator<<(std::ostream&, const Vector&);
		friend std::istream& operator>>(std::istream&, Vector&);

		void fill(double);
		void fillrand();
		void allocSpace();
		double norm();

	private:
		int cols_;
		double * data;
	};
}