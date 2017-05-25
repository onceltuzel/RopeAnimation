//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#ifndef _CO_VECTOR_H
#define _CO_VECTOR_H

#include "COVector.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


class COVector
{
public:

   // 
   //   Constructors
   //

   // default constructor
   COVector(void );

   // create and allocate: vector length
   COVector(int );
   
   // create, allocate and init: vector length, value
   COVector(int, double );
   
   // create and allocate: filename
   COVector(const char* );
      
   // copy constructor: vector
   COVector(const COVector&);

   // create, alocate and init: vector length, pointer to values
   COVector(int , double* );


   //
   // Destructors
   //
   
   // default destructor
   ~COVector(void);

   //
   // operators
   //
   
   // copy: vector
   const COVector& operator=(const COVector&);
   
   // test equal: vector
   int operator==(const COVector&) const;
   
   // test notequal: vector
   int operator!=(const COVector&) const;
   
   // add value: value
   COVector operator+(double) const;
   
   // add vector: vector
   COVector operator+(const COVector&) const;
   
   // add vector += : vector
   //COVector operator+=(const COVector&);
	COVector operator+=(const COVector&);
   
   // add value += : value
   COVector operator+=(double);
   
   // substract value: value
   COVector operator-(double) const;
   
   // substract vector: vector
   COVector operator-(const COVector&) const;
   
   // negate
   COVector operator-(void) const;
   
   // substract vector -= : vector
   COVector operator-=(const COVector&);
   
   // substract value -= : value
   COVector operator-=(double);
   
   // multiply value: value
   COVector operator*(double) const;
   
   // multiply vector: vector (element by element)
   COVector operator*(const COVector&) const;
      
   // multiply vector *= : vector (element by element)
   COVector operator*=(const COVector&);
   
   // multiply value *= : value
   COVector operator*=(double);
   
   // divide value: value
   COVector operator/(double) const;

   // divide vector: vector (element by element)
   COVector operator/(const COVector&) const;
   
   // divide value /= : value
   COVector operator/=(double);

   // divide vector /= : vector (element by element)
   COVector operator/=(const COVector&);
   
   //
   // Element access
   //
   
   // value of element (i): index
   inline double operator()(int ) const;
   
   // reference of element (i): index
   inline double& operator()(int );
   
   //
   // Memory functions
   //
   
   // test if vector is allocated
   inline int IsAllocated(void) const;
   
   // return data pointer
   inline double* GetData(void);
   
   // return data pointer (const)
   inline double* GetData(void) const;
   
   //
   // Vector functions
   //
   
   // get number of elements
   inline int GetLength(void) const;
   
   // get subvector (idxs, idxd): index start, index end (inclusive)
   COVector GetSubVector(int, int) const;

   // set subvector (idxs, idxd): index start, index end (inclusive)
   void SetSubVector(int , COVector& );
   
   //
   // Unary functions
   //

   // set vector from data
   void SetSameVector(double*);

   // set data from vector
   void GetSameVector(double*);
   
   // normalize vector (2-norm)
   void Normalize(void);
   void Normalize1Norm(void);

   // get norm
   double GetNorm(void);
   double GetNorm2(void);
   double GetNorm1(void);

   // multiply vector: vector (inner product)
   double InnerProd(const COVector&) const;
   // 2D Cross Product
   double Cross2D(const COVector&) const;

   // mean and variance
   double Mean() const;
   double Variance() const;
   double Variance(double) const;

   //
   // Structure modifications
   //
   
   // resize vector to (length): new length
   void Resize(int);
   
   // resize vector to (length) and init: vector length, value
   void Resize(int, double);
         
   // 
   // Input/Output
   //
   
   // write vector on filedescriptor: filedescriptor
   void Write(FILE*) const;
   
   // write vector on file: filename
   void Write(const char*) const;
   
   // read vector from filedescriptor: filedescriptor
   void Read(FILE*);
   
   // read vector from file: filename
   void Read(const char*);

	void Clear();
   
private:
    
    // vector is allocated
    int isAllocated_;
    
    // vector data
    double* data_;
    
    // length
    int length_;
        
    // allocate storage: length
    void PrivateAllocateData(int );
    
    // resize data: length
    void PrivateResize(int );
    
    // replace this by vector: vector.
    void PrivateCopyToThis(const COVector&);
    
    // clear storage
    void CleanData(void);                       
    
};

//
// Non member functions
//

// value * vector: value, vector
extern COVector operator*(double, const COVector&);

// value + vector: value, vector
extern COVector operator+(double, const COVector&);

// value - vector: value, vector
extern COVector operator-(double, const COVector&);

// vector * vector: vector, vector
extern COVector operator*(const COVector&, const COVector&);

inline int COVector::IsAllocated(void) const
{
   return isAllocated_;
}

inline int COVector::GetLength(void) const
{
   assert(IsAllocated());
   return length_;
}

inline double COVector::operator()(int idx) const
{
   assert(IsAllocated() && (idx >= 0) && (idx < GetLength()));
   return data_[idx];
}

inline double& COVector::operator()(int idx)
{
   assert(IsAllocated() && (idx >= 0) && (idx < GetLength()));
   return data_[idx];
}

inline double* COVector::GetData(void)
{
   assert(IsAllocated());
   return data_;
}

inline double* COVector::GetData(void) const
{
   assert(IsAllocated());
   return data_;
}

#endif
