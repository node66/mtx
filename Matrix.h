#pragma once
#include"Vector.h"
#include <iostream>

namespace algebra
{
	class Vector;

	class Matrix
	{
	public:
		Matrix();
		~Matrix();
		Matrix(int, int);
		Matrix(const Matrix&);

		int GetRows() { return rows_; }
		int GetCols() { return cols_; }
		void SetRows(int value) { value = rows_; }
		void SetCols(int value) { value = cols_; }

		double&operator()(int, int);
		Matrix&operator=(const Matrix&);
		Matrix&operator=(const double);
		Matrix&operator+=(const Matrix&);
		Matrix&operator-=(const Matrix&);
		Matrix&operator-=(double);
		Matrix&operator*=(const Matrix&);
		Matrix&operator*=(double);
		Matrix&operator/=(double);

		Matrix  operator-(const Matrix&);
	
		friend Matrix operator+(const Matrix&, const Matrix&);
		friend Matrix operator-(const Matrix&, const Matrix&);
		friend Matrix operator*(const Matrix&, const Matrix&);
		friend Matrix operator*(const Matrix&, double);
		friend Matrix operator-(const Matrix&, double);
		friend Matrix operator*(double, const Matrix&);
		friend Matrix operator/(const Matrix&, double);
		friend std::ostream& operator<<(std::ostream&, const Matrix&);
		friend std::istream& operator>>(std::istream&, Matrix&);

		void allocSpace();
		void fillrand();
		Matrix transpose();
		
		void LU(Matrix& L, Matrix& U);//LU ����������
		void LU(Matrix& L, Matrix&U, Vector&b);//LU ���������� c �������� ����
		void LDL_T(Matrix & L, Matrix &D);//LDL_T ����������
		void LDL_T(Matrix &L, Matrix &D, Vector&b); //LDL_T ���������� c �������� ����
		Matrix LU_Inversion(); //LU ��������� �������
		Matrix LDL_T_Inversion();//LDL_T ��������� �������
		void QR(Matrix& Q, Matrix& R);//QR ����������

	private:
		int rows_, cols_;
		double * data;
	};
}