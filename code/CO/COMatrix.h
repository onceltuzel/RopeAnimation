//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

/////////////////////////////////////////////////////////////////////////////

#ifndef _CO_MATRIX_H
#define _CO_MATRIX_H

#include "COVector.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "CODefaults.h"

#define CO_MATRIX_ASCII 0
#define CO_MATRIX_BIN64 1
#define CO_MATRIX_BIN32 2


class COMatrix
{
public:

   static double Determinant(double*, int);
   static bool Inverse(double*, int, double*);
   static void Transpose(double*, int, int, double*);
   static void MultMatVec(double*, int, int, double*, double*);
   static void MultMatMat(double*, int, int, double*, int, double*);
   static void CopyMat(double*, int, int, double*);
   static void OutProduct(double*, double* , int, double*);
   static void ApproximateRotation(double*, double*);

   double ApproximateRotation(void);
  
   enum COMatrixType {CODiagonal,
      COTranslationX, COTranslationY, COTranslationZ,
      CORotationX, CORotationY, CORotationZ};
   
   // 
   //   Constructors
   //

   // default constructor
   COMatrix(void);

   // create and allocate: rows, cols
   COMatrix(int, int);
   
   // create, allocate and init matrix: rows, cols, value
   COMatrix(int, int, double);
   
   // create and allocate: filename
   COMatrix(const char*);
   
   // create, allocate and init projective 4x4 matrix: type, value (diag, translation, angle)
   COMatrix(COMatrixType, double);
   
   // copy constructor: matrix
   COMatrix(const COMatrix&);

   // create, alocate and init matrix: rows, cols, pointer to values
   COMatrix(int, int, double*);


   //
   // Destructors
   //
   
   // default destructor
   ~COMatrix(void);

   //
   // operators
   //
   
   // copy: matrix
   const COMatrix& operator=(const COMatrix&);
   
   // test equal: matrix
   int operator==(const COMatrix&) const;
   
   // test notequal: matrix
   int operator!=(const COMatrix&) const;
   
   // add value: value
   COMatrix operator+(double) const;
   
   // add matrix: matrix
   COMatrix operator+(const COMatrix&) const;
   
   // add matrix += : matrix
   COMatrix operator+=(const COMatrix&);
   
   // add value += : value
   COMatrix operator+=(double);
   
   // substract value: value
   COMatrix operator-(double) const;
   
   // substract matrix: matrix
   COMatrix operator-(const COMatrix&) const;
   
   // negate
   COMatrix operator-(void) const;
   
   // substract matrix -= : matrix
   COMatrix operator-=(const COMatrix&);
   
   // substract value -= : value
   COMatrix operator-=(double);
   
   // multiply value: value
   COMatrix operator*(double) const;
   
   // multiply matrix: matrix
   COMatrix operator*(const COMatrix&) const;
   
   // multiply vector: vector
   COVector operator*(const COVector&) const;

   // multiply vector *= : vector
   COMatrix operator*=(const COVector&);
   
   // multiply matrix *= : matrix
   COMatrix operator*=(const COMatrix&);
   
   // multiply value *= : value
   COMatrix operator*=(double);
   
   // divide value: value
   COMatrix operator/(double) const; 
   
   // divide value /= : value
   COMatrix operator/=(double);
   
   //
   // Element access
   //
   
   // value of element (r,c): row, column
   inline double operator()(int, int) const;
   
   // reference of element (r,c): row, column
   inline double& operator()(int, int);
   
   //
   // Memory functions
   //
   
   // test if matrix is allocated
   inline int IsAllocated(void) const;
   
   // return data pointer
   inline double* GetData(void);
   
   // return data pointer (const) 
   inline double* GetData(void) const;
   
   //
   // Matrix dimensions
   //
   
   // get number of rows
   inline int GetNRows(void) const;
   
   // get number of columns
   inline int GetNCols(void) const;
   
   // get number of elements
   inline int GetNElements(void) const;

   // get submatrix (fr-lr) to (fc-lc): first raw, last row, first column, last column
   COMatrix GetSubMatrix(int, int, int, int) const;
   void SetSubMatrix(int , int , COMatrix& );
   
   //
   // Unary functions
   //

   COVector GetColumn(int ) const;
   COVector GetRow(int ) const;
   void SetColumn(int, COVector& );
   void SetRow(int, COVector& );
   void SetDiagonal( double );
   void SetDiagonal( COVector& );
   void SetIdentity();

   // set matrix from data
   void SetSameMatrix(double*);

   // set matrix from data
   void GetSameMatrix(double*);
   
   // normalize matrix by the last element
   void Normalize(void);

   // normalize matrix by frobenious norm
   void NormalizeFrobenius(void);

   // normalize projection matrix
   void NormalizeProjection(void);
   
   // compute inverse
   COMatrix Inverse(void) const;

   // compute pseudo inverse
   COMatrix PInverse(void) const;
   COMatrix PInverse(int imprank) const;

   // compute transpose
   COMatrix Transpose(void) const;
  
   double NormFrobenius(void);

   double Trace(void);

   COVector VarianceCols() const;
   COVector VarianceCols(int) const;
   
   //
   // Structure modifications
   //
   
   // resize matrix to (nr,nc): rows, columns
   void Resize(int, int);
   
   // resize matrix to (nr, nc) and init with value: rows, columns, value
   void Resize(int, int, double);
   
   //
   // Linear algebra functions
   //
   
   // compute determinant
   double Determinant(void) const;
   
   // compute SVD m = u * d * Transpose(v): u, d, v, u type=full ("A", "N"), v type=full ("A", "N")
   void Svd(COMatrix&, COMatrix&, COMatrix&, char* uType="A", char* vType="A") const;
   void SvdSpec(int nr, int nc, COMatrix& d) const;

   // compute LU decompozition of the matrix (with pivoting)
   void LUDecomp(COMatrix& LU, long* permidx) const;
   void LUDecomp(COMatrix& L, COMatrix& U, COMatrix& P) const;

   // compute RQ factorization 
   void RQDecomp(COMatrix& RQ, COVector& Tau) const;
   void RQDecomp(COMatrix& R, COMatrix& Q) const;

   // decompose matrix such that M = L'*L (M must be symetric)
   COMatrix Decomp();

   // special rotation matrix
   void GetRotationAngles(double* alpha, double* beta, double* gamma);
   void SetRotationAngles(double alpha, double beta, double gamma);

   void GetRotationQuaternion(COVector&);
   void GetRotationQuaternion(double* q0, double* q1, double* q2, double* q3);
   void SetRotationQuaternion(COVector&);
   void SetRotationQuaternion(double , double , double , double );

   COMatrix CorrMat();


   // 
   // Input/Output
   //
   
   // write matrix on filedescriptor: filedescriptor
   void Write(FILE*, int mtype=CO_MATRIX_ASCII) const;
   //void LogFile() const;
   
   // write matrix on file: filename
   void Write(const char*, int mtype=CO_MATRIX_ASCII) const;
   
   // read matrix from filedescriptor: filedescriptor
   void Read(FILE*);
   
   // read matrix from file: filename
   void Read(const char*);

   void SetPersistentData(double* data, int nrows, int ncols);
   
private:
    
    // matrix is allocated
    int isAllocated_;

    // data is persistent
    int isPersistent_;
    
    // matrix data
    double* data_;
    
    // rows
    int nRows_;
    
    // columns
    int nCols_;
    
    // allocate storage: rows, columns
    void PrivateAllocateData(int, int);
    
    // resize data: rows, columns
    void PrivateResize(int, int);
    
    // replace this by matrix: matrix.
    void PrivateCopyToThis(const COMatrix&);
    
    // clear storage
    void CleanData(void);                       
    
};

//
// Non member functions
//

// value * matrix: value, matrix
extern COMatrix operator*(double, const COMatrix&);

// value + matrix: value, matrix
extern COMatrix operator+(double, const COMatrix&);

// value - matrix: value, matrix
extern COMatrix operator-(double, const COMatrix&);

// vector * matrix: vector, matrix
extern COMatrix operator*(const COVector&, const COMatrix&);

extern COMatrix OuterProduct(const COVector&, const COVector&);

inline int COMatrix::IsAllocated(void) const
{
   return isAllocated_;
}

inline int COMatrix::GetNRows(void) const
{
   assert(IsAllocated());
   return nRows_;
}

inline int COMatrix::GetNCols(void) const
{
   assert(IsAllocated());
   return nCols_;
}

inline int COMatrix::GetNElements(void) const
{
   assert(IsAllocated());
   return nRows_*nCols_;
}

inline double COMatrix::operator()(int r, int c) const
{
   assert(IsAllocated() && (r >= 0) && (r < GetNRows()) && (c >= 0) && (c < GetNCols()));
   return data_[r+c*nRows_];
}

inline double& COMatrix::operator()(int r, int c)
{
   assert(IsAllocated() && (r >= 0) && (r < GetNRows()) && (c >= 0) && (c < GetNCols()));
   return data_[r+c*nRows_];
}

inline double* COMatrix::GetData(void)
{
   assert(IsAllocated());
   return data_;
}

inline double* COMatrix::GetData(void) const
{
   assert(IsAllocated());
   return data_;
}

extern void COLogFile(COMatrix&);
extern void COLogFile(COVector&);

#endif
