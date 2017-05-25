//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//

#include "COMatrix.h"

#define MATRIX_ZERO_TRESH 0.0000000001
#define SING_TRESH 0.00000000001
#define SING_TRESH_ENERGY 0.99999

#define _CO_USELAPACK_H
#ifdef _CO_USELAPACK_H

#include "f2c.h"
extern "C" int dgesvd_(char *, char *, integer *, integer *, 
                       doublereal *, integer *, doublereal *, doublereal *, integer *, 
                       doublereal *, integer *, doublereal *, integer *, integer *);

extern "C" int dgetrf_(integer *, integer *, doublereal *, integer *,
                       integer *, integer *);
extern "C" int dgerqf_(integer *, integer *, doublereal *, integer *,
	                    doublereal *, doublereal *, integer *, integer *);

#else
typedef long int integer;
typedef double doublereal;
int dgesvd_(char *, char *, integer *, integer *, 
                       doublereal *, integer *, doublereal *, doublereal *, integer *, 
                       doublereal *, integer *, doublereal *, integer *, integer *)
{
   printf("Please compile with lapack.\n");
   exit(-1);
   return -1;
}
int dgetrf_(integer *, integer *, doublereal *, integer *,
                       integer *, integer *)
{
   printf("Please compile with lapack.\n");
   exit(-1);
   return -1;
}
int dgerqf_(integer *, integer *, doublereal *, integer *,
	         doublereal *, doublereal *, integer *, integer *)
{
   printf("Please compile with lapack.\n");
   exit(-1);
   return -1;
}

#endif



double COMatrix::Determinant(double* mat, int n)
{
   if (n==2)
   {
      return (mat[0]*mat[3]-mat[1]*mat[2]);
   }
   if (n==3)
   {
      return (mat[0]*mat[4]*mat[8]+mat[3]*mat[7]*mat[2]+mat[1]*mat[5]*mat[6]-
         mat[2]*mat[4]*mat[6]-mat[1]*mat[3]*mat[8]-mat[0]*mat[5]*mat[7]);
   }
   return 0;
}

bool COMatrix::Inverse(double* inmat, int n, double* outmat)
{
   double det=0;
   if (n==2 || n==3) {
      det = COMatrix::Determinant(inmat, n);
      if (fabs(det)<MATRIX_ZERO_TRESH) // near singular
         return false;
   }
   if (n==2)
   {
      outmat[0] = inmat[3]/det;
      outmat[2] = -inmat[2]/det;
      outmat[1] = -inmat[1]/det;
      outmat[3] = inmat[0]/det;
      return true;
   }
   if (n==3)
   {
      outmat[0] = (inmat[4]*inmat[8]-inmat[5]*inmat[7])/det;
      outmat[3] = -(inmat[3]*inmat[8]-inmat[5]*inmat[6])/det;
      outmat[6] = (inmat[3]*inmat[7]-inmat[4]*inmat[6])/det;
      outmat[1] = -(inmat[1]*inmat[8]-inmat[2]*inmat[7])/det;
      outmat[4] = (inmat[0]*inmat[8]-inmat[2]*inmat[6])/det;
      outmat[7] = -(inmat[0]*inmat[7]-inmat[1]*inmat[6])/det;
      outmat[2] = (inmat[1]*inmat[5]-inmat[2]*inmat[4])/det;
      outmat[5] = -(inmat[0]*inmat[5]-inmat[2]*inmat[3])/det;
      outmat[8] = (inmat[0]*inmat[4]-inmat[1]*inmat[3])/det;
      return true;
   }

   integer svdldu = n;
   integer svdlda = n;
   integer svdldvt = n;
   integer svdlwork = n*n;
   integer svdinfo;
   integer svdm = n;
   integer svdn = n;
   double* svdu = new double[n*n];
   double* svdvt = new double[n*n];
   double* svdwork = new double[n*n];
   double* svds = new double[n*n];
   double* svdinmat = new double[n*n];
   COMatrix::CopyMat(inmat, n, n, svdinmat);

   dgesvd_("A", "A", &svdm, &svdn, svdinmat, &svdlda, svds, svdu, &svdldu,
      svdvt, &svdldvt, svdwork, &svdlwork, &svdinfo);
   
   // make inverse diagonal svds
   int i;
   for (i=0; i<n*n; i++)
      svdwork[i] = 0;
   for (i=0; i<n; i++)
   {
      if (fabs(svds[i])>MATRIX_ZERO_TRESH)
         svdwork[i+i*n] = 1/svds[i];
   }

   COMatrix::Transpose(svdvt, n, n, svds);
   COMatrix::MultMatMat(svds, n, n, svdwork, n, svdvt);
   COMatrix::Transpose(svdu, n, n, svds);
   COMatrix::MultMatMat(svdvt, n, n, svds, n, outmat);

   delete svdinmat;
   delete svds;
   delete svdwork;
   delete svdvt;
   delete svdu;

   return false;

}

void COMatrix::Transpose(double* inmat, int m, int n, double* outmat)
{
   int c, r;
   for (c=0; c<n; c++)
   {
      for (r=0; r<m; r++)
      {
         outmat[c+r*n]=inmat[r+c*m];
      }
   } 
}

void COMatrix::MultMatVec(double* mat, int m, int n, double* vec, double* rez)
{
   int r, c;
   for (r=0; r<m; r++)
   {
      rez[r] = 0;
      for (c=0; c<n; c++)
      {
         rez[r] += mat[r+c*m]*vec[c];
      }
   }
}

void COMatrix::MultMatMat(double* mat1, int m, int n, double* mat2, int p, double* rez)
{
   int r, c, i;
   for (r=0; r<m; r++)
   {
      for (c=0; c<p; c++)
      {
         rez[r+c*m]=0;
         for (i=0; i<n; i++)
         {
            rez[r+c*m] += mat1[r+i*m]*mat2[i+c*n];
         }
      }
   }
}

void COMatrix::CopyMat(double* inmat, int m, int n, double* outmat)
{
   int i;
   for (i=0; i<m*n; i++)
      outmat[i] = inmat[i];
}

void COMatrix::OutProduct(double* v1, double* v2, int n, double* mat)
{
   int l, c;
   for (l=0; l<n; l++)
   {
      for (c=0; c<n; c++)
      {
         mat[l+c*n] = v1[l]*v2[c];
      }
   }
}

void COMatrix::ApproximateRotation(double* inmat, double* outmat)
{
   integer svdldu = 3;
   integer svdlda = 3;
   integer svdldvt = 3;
   integer svdlwork = 36;
   integer svdinfo;
   integer svdm = 3;
   integer svdn = 3;
   double svdu[9];
   double svdvt[9];
   double svdwork[36];
   double svds[9];
   double svdinmat[9];
   COMatrix::CopyMat(inmat, 3, 3, svdinmat);

   dgesvd_("A", "A", &svdm, &svdn, svdinmat, &svdlda, svds, svdu, &svdldu,
      svdvt, &svdldvt, svdwork, &svdlwork, &svdinfo);

   COMatrix::MultMatMat(svdu, 3, 3, svdvt, 3, svdinmat);
   double detuvt = COMatrix::Determinant(svdinmat, 3);
   for (int i=0; i<9; i++)
      svds[i] = 0;
   svds[0] = svds[4] = 1;
   svds[8] = detuvt;

   COMatrix::MultMatMat(svdu, 3, 3, svds, 3, svdinmat);
   COMatrix::MultMatMat(svdinmat, 3, 3, svdvt, 3, outmat);
}

double COMatrix::ApproximateRotation()
{
   assert((nCols_ == 3) && (nRows_ == 3));
   double retval;
   COMatrix u(3,3);
   COMatrix v(3,3);
   COMatrix d(3,3);
   Svd(u, d, v, "A", "A");
   retval = d(0,0);
   d(0,0) = 1.0;
   d(1,1) = 1.0;
   d(2,2) = 1.0/(u.Determinant()*v.Determinant());
   *this = u*d*v.Transpose();
   return retval;
}

/****************************************************************************/

void COMatrix::PrivateAllocateData(int nr, int nc)
{
  assert((nr > 0) && (nc > 0));

  nRows_ = nr;
  nCols_ = nc;
  data_ = new double[nRows_*nCols_];
  assert(data_ != 0);
  isAllocated_ = 1;
  isPersistent_ = 0;
}

void COMatrix::PrivateResize(int nr, int nc)
{
  if (!IsAllocated())
    PrivateAllocateData(nr, nc);
  //else if ((GetNRows() != nr) || (GetNCols() != nc))
  else if(GetNElements() != (nr*nc))
  {
    CleanData();
    PrivateAllocateData(nr, nc);
  } else
  {
     nRows_ = nr;
     nCols_ = nc;
  }
}

void COMatrix::PrivateCopyToThis(const COMatrix& m)
{
  PrivateResize(m.GetNRows(), m.GetNCols());

  int i;
  double* indata;
  indata = m.GetData();
  for (i=0; i<GetNElements(); i++)
     data_[i] = indata[i];
}

void COMatrix::CleanData(void)
{
   if (isPersistent_)
      data_ = 0;
   else
      delete [] data_;
   isAllocated_ = 0;
}

COMatrix::COMatrix(void)
{
   isPersistent_ = 0;
   isAllocated_ = 0;
   data_ = 0;
}

COMatrix::COMatrix(int nr, int nc)
{
   assert((nr > 0) && (nc > 0));
   PrivateAllocateData(nr, nc);
}

COMatrix::COMatrix(int nr, int nc, double val)
{
   assert((nr > 0) && (nc > 0));
   PrivateAllocateData(nr, nc);
   
   int i;
   for (i=0; i<nr*nc; i++)
      data_[i] = val;
}

COMatrix::COMatrix(int nr, int nc, double* ival)
{
   assert((nr > 0) && (nc > 0));
   PrivateAllocateData(nr, nc);
   int i;
   for (i=0; i<nr*nc; i++)
      data_[i] = ival[i];
}

void COMatrix::SetSameMatrix(double* ival)
{
   int i;
   for (i=0; i<nCols_*nRows_; i++)
      data_[i] = ival[i];
}

void COMatrix::GetSameMatrix(double* ival)
{
   int i;
   for (i=0; i<nCols_*nRows_; i++)
      ival[i] = data_[i];
}

COMatrix::COMatrix(const char* name)
{
   isAllocated_ = 0;
   Read(name);
}

COMatrix::COMatrix(COMatrixType type, double val)
{
   PrivateAllocateData(4, 4);
   
   switch(type)
   {
      
   case CODiagonal:
      (*this)(1, 1) = val; (*this)(1, 2) =   0; (*this)(1, 3) =   0; (*this)(1, 4) =   0;
      (*this)(2, 1) =   0; (*this)(2, 2) = val; (*this)(2, 3) =   0; (*this)(2, 4) =   0;
      (*this)(3, 1) =   0; (*this)(3, 2) =   0; (*this)(3, 3) = val; (*this)(3, 4) =   0;
      (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = val;
      break;
      
   case COTranslationX:
      (*this)(1, 1) = 1.0; (*this)(1, 2) =   0; (*this)(1, 3) =   0; (*this)(1, 4) = val;
      (*this)(2, 1) =   0; (*this)(2, 2) = 1.0; (*this)(2, 3) =   0; (*this)(2, 4) =   0;
      (*this)(3, 1) =   0; (*this)(3, 2) =   0; (*this)(3, 3) = 1.0; (*this)(3, 4) =   0;
      (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
      break;
      
   case COTranslationY:
      (*this)(1, 1) = 1.0; (*this)(1, 2) =   0; (*this)(1, 3) =   0; (*this)(1, 4) =   0;
      (*this)(2, 1) =   0; (*this)(2, 2) = 1.0; (*this)(2, 3) =   0; (*this)(2, 4) = val;
      (*this)(3, 1) =   0; (*this)(3, 2) =   0; (*this)(3, 3) = 1.0; (*this)(3, 4) =   0;
      (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
      break;
      
   case COTranslationZ:
      (*this)(1, 1) = 1.0; (*this)(1, 2) =   0; (*this)(1, 3) =   0; (*this)(1, 4) =   0;
      (*this)(2, 1) =   0; (*this)(2, 2) = 1.0; (*this)(2, 3) =   0; (*this)(2, 4) =   0;
      (*this)(3, 1) =   0; (*this)(3, 2) =   0; (*this)(3, 3) = 1.0; (*this)(3, 4) = val;
      (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
      break;
      
   case CORotationX:
      {
         double c = cos(val);
         double s = sin(val);
         (*this)(1, 1) = 1.0; (*this)(1, 2) =   0; (*this)(1, 3) =   0; (*this)(1, 4) =   0;
         (*this)(2, 1) =   0; (*this)(2, 2) =   c; (*this)(2, 3) = - s; (*this)(2, 4) =   0;
         (*this)(3, 1) =   0; (*this)(3, 2) =   s; (*this)(3, 3) =   c; (*this)(3, 4) =   0;
         (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
         break;
      }
      
   case CORotationY:
      {
         double c = cos(val);
         double s = sin(val);
         (*this)(1, 1) =   c; (*this)(1, 2) =   0; (*this)(1, 3) =   s; (*this)(1, 4) =   0;
         (*this)(2, 1) =   0; (*this)(2, 2) = 1.0; (*this)(2, 3) =   0; (*this)(2, 4) =   0;
         (*this)(3, 1) = - s; (*this)(3, 2) =   0; (*this)(3, 3) =   c; (*this)(3, 4) =   0;
         (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
         break;
      }
      
   case CORotationZ:
      {
         double c = cos(val);
         double s = sin(val);
         (*this)(1, 1) =   c; (*this)(1, 2) = - s; (*this)(1, 3) =   0; (*this)(1, 4) =   0;
         (*this)(2, 1) =   s; (*this)(2, 2) =   c; (*this)(2, 3) =   0; (*this)(2, 4) =   0;
         (*this)(3, 1) =   0; (*this)(3, 2) =   0; (*this)(3, 3) = 1.0; (*this)(3, 4) =   0;
         (*this)(4, 1) =   0; (*this)(4, 2) =   0; (*this)(4, 3) =   0; (*this)(4, 4) = 1.0;
         break;
      }
      
   default:
      assert(0);
      
   }
}

COMatrix::COMatrix(const COMatrix& m)
{
   isAllocated_ = 0;
   if (!m.IsAllocated())
      return;
   PrivateCopyToThis(m);
}

const COMatrix& COMatrix::operator=(const COMatrix& m)
{
   if (this == &m)
      return *this;
   
   if (!m.IsAllocated())
   {
      isAllocated_ = 0;
      return *this;
   }
   
   PrivateCopyToThis(m);
   return *this;
}

COMatrix::~COMatrix(void) 
{
   if (IsAllocated())
      CleanData();
}

void COMatrix::SetPersistentData(double* data, int nr, int nc)
{
   if (IsAllocated())
      CleanData();
   data_ = data;
   isAllocated_ = 1;
   nRows_ = nr;
   nCols_ = nc;
   isPersistent_ = 1;
}

int COMatrix::operator==(const COMatrix& m) const
{
   assert((IsAllocated()) && (m.IsAllocated()));
   
   if ((GetNRows() != m.GetNRows()) || (GetNCols() != m.GetNCols()))
      return 0;
   int i;
   double* datam;
   datam = m.GetData();
   for (i=0; i<GetNElements(); i++)
   {
      if (data_[i] != datam[i])
      {
         return 0;
      }
   }
   return 1;
}

int COMatrix::operator!=(const COMatrix& m) const
{
   return !(*this == m);
}

COMatrix COMatrix::operator+(double d) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNRows(), GetNCols());
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]+d;      
   return r;
}

COMatrix COMatrix::operator+(const COMatrix& m) const
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNRows() == m.GetNRows()) && (GetNCols() == m.GetNCols()));
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   double* datam;
   datar = r.GetData();
   datam = m.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]+datam[i];
   return r;
}

COMatrix operator+(double d, const COMatrix& m)
{
   assert(m.IsAllocated());
   
   COMatrix r(m.GetNRows(), m.GetNCols());
   
   int i;
   double* datar;
   double* datam;
   datar = r.GetData();
   datam = m.GetData();
   for (i=0; i<m.GetNElements(); i++)
      datar[i] = d+datam[i];      
   return r;
}

COMatrix COMatrix::operator+=(const COMatrix& m)
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNRows() == m.GetNRows()) && (GetNCols() == m.GetNCols()));
   
   int i;
   double* datam;
   datam = m.GetData();
   for (i=0; i<GetNElements(); i++)
      data_[i] += datam[i];      
   return *this;
}

COMatrix COMatrix::operator+=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] +=d;
   return *this;
}

COMatrix COMatrix::operator-(double d) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]-d;
      
   return r;
}

COMatrix COMatrix::operator-(const COMatrix& m) const
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNRows() == m.GetNRows()) && (GetNCols() == m.GetNCols()));
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   double* datam;
   datar = r.GetData();
   datam = m.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]-datam[i];
   return r;
}

COMatrix operator-(double d, const COMatrix& m)
{
   assert(m.IsAllocated());
   
   COMatrix r(m.GetNRows(), m.GetNCols());
   
   int i;
   double* datar;
   double* datam;
   datar = r.GetData();
   datam = m.GetData();
   for (i=0; i<m.GetNElements(); i++)
      datar[i] = d-datam[i];
   return r;
}

COMatrix COMatrix::operator-(void) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = -data_[i];
   return r;
}

COMatrix COMatrix::operator-=(const COMatrix& m)
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNRows() == m.GetNRows()) && (GetNCols() == m.GetNCols()));
   
   int i;
   double* datam;
   datam = m.GetData();
   for (i=0; i <GetNElements(); i++)
      data_[i] -= datam[i];
   return *this;
}

COMatrix COMatrix::operator-=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] -=d;
   return *this;
}

COMatrix COMatrix::operator*(double d) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]*d;
   return r;
}

COMatrix COMatrix::operator*(const COMatrix& m) const
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNCols() == m.GetNRows()));
   
   COMatrix r(GetNRows(), m.GetNCols());
   
   int i, j, k;
   double s;
   for (i=0; i<GetNRows(); i++)
   {
      for (j=0; j<m.GetNCols(); j++)
      {
         s = 0;
         for (k=0; k<GetNCols(); k++)
         {
            s += (*this)(i, k) * m(k, j);
         }
         r(i, j) = s;
      }
   }
   
   return r;
}

COVector COMatrix::operator*(const COVector& v) const
{
   assert(IsAllocated() && (v.IsAllocated()) && 
      (GetNCols() == v.GetLength()));

   COVector r(GetNRows(), 0.0);
   int i, j;
   for (i=0; i<nRows_; i++)
   {
      for (j=0; j<nCols_; j++)
      {
         r(i) += (*this)(i,j) * v(j);
      }
   }
   return r;
}

COMatrix operator*(double d, const COMatrix& m)
{
   assert(m.IsAllocated());
   
   COMatrix r(m.GetNRows(), m.GetNCols());
   
   int i;
   double* datar;
   double* datam;
   datar = r.GetData();
   datam = m.GetData();
   for (i=0; i<m.GetNElements(); i++)
      datar[i] = d*datam[i];
   return r;
}

COMatrix COMatrix::operator*=(const COMatrix& m)
{
   assert((IsAllocated()) && (m.IsAllocated()) &&
      (GetNCols() == m.GetNRows()));
   
   COMatrix r(GetNRows(), m.GetNCols());
   
   int i, j, k;
   double s;
   for (i=0; i<GetNRows();i++)
   {
      for (j=0; j<m.GetNCols();j++)
      {
         s = 0;
         for (k=0; k <GetNCols(); k++)
         {
            s += (*this)(i, k) * m(k, j);
         }
         r(i, j) = s;
      }
   }
   
   PrivateCopyToThis(r);
   return *this;
}

COMatrix COMatrix::operator*=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] *= d;
   return *this;
}

COMatrix COMatrix::operator/(double d) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNRows(), GetNCols());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetNElements(); i++)
      datar[i] = data_[i]/d;
   return r;
}

COMatrix COMatrix::operator/=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] /= d;
   return *this;
}

COMatrix COMatrix::GetSubMatrix(int fr, int lr, int fc, int lc) const
{
   assert((IsAllocated()) &&
      (fr >= 0) && (fr < GetNRows()) && (lr >= fr) && (lr < GetNRows()) &&
      (fc >= 0) && (fc < GetNCols()) && (lc >= fc) && (lc < GetNCols()));
   
   COMatrix r(lr - fr + 1, lc - fc + 1);
   
   int i, j;
   for (i=fr; i<=lr; i++)
      for (j=fc; j<=lc; j++)
         r(i-fr, j-fc) = (*this)(i, j);
      
   return r;
}

COVector COMatrix::GetColumn(int c) const
{
   assert((IsAllocated()) && (c>=0) && (c<nCols_));
   COVector v(nRows_);
   int i;
   double* vdata;
   vdata = v.GetData();
   for (i=0; i<nRows_; i++)
      vdata[i] = (*this)(i, c);

   return v;
}

COVector COMatrix::GetRow(int r) const
{
   assert((IsAllocated()) && (r>=0) && (r<nRows_));
   COVector v(nCols_);
   int i;
   double* vdata;
   vdata = v.GetData();
   for (i=0; i<nCols_; i++)
      vdata[i] = (*this)(r, i);

   return v;
}

void COMatrix::SetColumn(int c, COVector& v)
{
   assert((IsAllocated()) && (c>=0) && (c<nCols_) && (nRows_ == v.GetLength()));
   int i;
   double* vdata;
   vdata = v.GetData();
   for (i=0; i<nRows_; i++)
      (*this)(i, c) = vdata[i];
}

void COMatrix::SetRow(int r, COVector& v)
{
   assert((IsAllocated()) && (r>=0) && (r<nRows_) && (nCols_ == v.GetLength()));
   int i;
   double* vdata;
   vdata = v.GetData();
   for (i=0; i<nCols_; i++)
      (*this)(r, i) = vdata[i];
}

void COMatrix::SetDiagonal(double v)
{
   assert(IsAllocated());
   int i;
   memset(data_, 0, sizeof(double)*nCols_*nRows_);

   for (i=0; i<nCols_; i++)
      (*this)(i, i) = v;
}


void COMatrix::SetDiagonal(COVector& v)
{
   assert((IsAllocated())  && (nCols_ == v.GetLength()));
   int i;
   double* vdata;

   memset(data_, 0, sizeof(double)*nCols_*nRows_);

   vdata = v.GetData();
   for (i=0; i<nCols_; i++)
      (*this)(i, i) = vdata[i];
}


void COMatrix::SetIdentity()
{
   assert(IsAllocated());

   int i;
   for (i=0; i<nCols_; i++)
      (*this)(i, i) = 1.0;
}

void COMatrix::SetSubMatrix(int fr, int fc, COMatrix& dm)
{
   int lr, lc;
   lr = fr+dm.GetNRows()-1;
   lc = fc+dm.GetNCols()-1;
   assert((IsAllocated()) &&
      (fr >= 0) && (fr < GetNRows()) && (lr >= fr) && (lr < GetNRows()) &&
      (fc >= 0) && (fc < GetNCols()) && (lc >= fc) && (lc < GetNCols()));

   int i,j;
   for (i=fr; i<=lr; i++)
      for (j=fc; j<=lc; j++)
         (*this)(i, j) = dm(i-fr, j-fc);
}

void COMatrix::Normalize(void)
{
   assert(IsAllocated());
   
   double lastElement = data_[GetNElements()-1];
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] /= lastElement;
      
   data_[GetNElements()-1] = 1;
}

void COMatrix::NormalizeFrobenius(void)
{
   assert(IsAllocated());
   
   double normf = NormFrobenius();
   if (normf < MATRIX_ZERO_TRESH)
      return;
   int i;
   for (i=0; i<GetNElements(); i++)
      data_[i] /= normf;
}

void COMatrix::NormalizeProjection(void)
{
   assert(IsAllocated() && (nCols_ == 4) && (nRows_ == 3));

   double det,norm;
   int i;
   det= data_[0]*data_[4]*data_[8]+data_[2]*data_[3]*data_[7]+data_[1]*data_[5]*data_[6]-
      data_[2]*data_[4]*data_[6]-data_[1]*data_[3]*data_[8]-data_[0]*data_[5]*data_[7];
   norm = sqrt(data_[2]*data_[2]+data_[5]*data_[5]+data_[8]*data_[8]);
   if (det < 0)
      norm = -norm;
   if (fabs(norm) < 0.0000001)
      return;
   for (i=0; i<12; i++)
      data_[i] /= norm;
}

COMatrix COMatrix::Inverse(void) const
{
   assert((IsAllocated()) && (GetNRows() == GetNCols()));
   
   COMatrix u(GetNRows(), GetNRows());
   COMatrix d(GetNRows(), GetNCols(), (double) 0);
   COMatrix v(GetNCols(), GetNCols());
   Svd(u, d, v, "A", "A");
   
   int i;
   for (i=0; i<d.GetNRows(); i++)
   {
      assert(d(i, i) != 0);
      d(i, i) = 1/d(i, i);
   }
   
   return v * d * u.Transpose();
}

COMatrix COMatrix::PInverse(void) const
{
   assert((IsAllocated()));// && (GetNRows() == GetNCols()));
   
   COMatrix u(GetNRows(), GetNRows());
   COMatrix d(GetNRows(), GetNCols(), (double) 0);
   COMatrix v(GetNCols(), GetNCols());
   Svd(u, d, v, "A", "A");
   
   int i;
   double diag = 0;
   int mindim = COMin(d.GetNRows(), d.GetNCols());
   for (i=0; i<mindim; i++)
   {
      diag += d(i,i);
   }
   if (diag<SING_TRESH)
   {
      d.Resize(GetNRows(), GetNCols(), 0.0);
      return d;
   }
   for (i=0; i<mindim; i++)
   {
      if (fabs(d(i,i)/diag)<SING_TRESH)
         d(i,i) = 0;
      else
         d(i, i) = 1/d(i, i);
   }
   
   return v * d.Transpose() * u.Transpose();
}

/*
COMatrix COMatrix::PInverse(void) const
{
   assert((IsAllocated()));// && (GetNRows() == GetNCols()));
   
   COMatrix u(GetNRows(), GetNRows());
   COMatrix d(GetNRows(), GetNCols(), (double) 0);
   COMatrix v(GetNCols(), GetNCols());
   Svd(u, d, v, "A", "A");
   
   int i;
   double diag = 0;
   int mindim = COMin(d.GetNRows(), d.GetNCols());
   for (i=0; i<mindim; i++)
   {
      diag += d(i,i);
   }
   if (diag<SING_TRESH)
   {
      d.Resize(GetNRows(), GetNCols(), 0.0);
      return d;
   }
   double cdiagsum = 0;
   for (i=0; i<mindim; i++)
   {
      cdiagsum += d(i,i);
      if ((cdiagsum/diag) > SING_TRESH_ENERGY)
         d(i,i) = 0;
      else
         d(i, i) = 1/d(i, i);
   }
   
   return v * d.Transpose() * u.Transpose();
}
*/

COMatrix COMatrix::PInverse(int imprank) const
{
   assert((IsAllocated()) && (GetNRows() == GetNCols()));
   
   COMatrix u(GetNRows(), GetNRows());
   COMatrix d(GetNRows(), GetNCols(), (double) 0);
   COMatrix v(GetNCols(), GetNCols());
   Svd(u, d, v, "A", "A");
   
   int i;
   int mindim = COMin(d.GetNRows(), d.GetNCols());
   imprank = (imprank < mindim) ? imprank : mindim;
   for (i=0; i<imprank; i++) {
      assert(d(i, i) != 0);
      d(i, i) = 1/d(i, i);
   }
   for (i=imprank; i<mindim; i++)
      d(i, i) = 0;
   
   return v * d.Transpose() * u.Transpose();
}

double COMatrix::NormFrobenius(void)
{
   assert(IsAllocated());
   double nfro = 0;
   int i;
   for (i=0; i<(nCols_*nRows_); i++)
      nfro += data_[i]*data_[i];
   nfro = sqrt(nfro);
   return nfro;
}

double COMatrix::Trace() {
   assert(IsAllocated());
   int mnd = (nCols_ < nRows_) ? nCols_ : nRows_;
   double trc = 0;
   int i;
   for (i=0; i<mnd; i++)
      trc += (*this)(i, i);
   return trc;
}

COVector COMatrix::VarianceCols() const
{
   assert(IsAllocated());
   COVector v(GetNRows());
   double mean;
   double var;
   int r,c;
   for (r=0; r<nRows_; r++)
   {
      mean = 0;
      for (c=0; c<nCols_; c++)
         mean += data_[r+c*nRows_];
      mean /= nCols_;
      var = 0;
      for (c=0; c<nCols_; c++)
         var += (data_[r+c*nRows_]-mean)*(data_[r+c*nRows_]-mean);
      v(r) = var/nCols_;
   }
   return v;
}

COVector COMatrix::VarianceCols(int col) const
{
   assert(IsAllocated());
   COVector v(GetNRows());
   double mean;
   double var;
   int r,c;
   for (r=0; r<nRows_; r++)
   {
      mean = data_[r+col*nRows_];
      var = 0;
      for (c=0; c<nCols_; c++)
         var += (data_[r+c*nRows_]-mean)*(data_[r+c*nRows_]-mean);
      v(r) = var/nCols_;
   }
   return v;
}

COMatrix COMatrix::Transpose(void) const
{
   assert(IsAllocated());
   
   COMatrix r(GetNCols(), GetNRows());
   
   int i, j;
   for (i=0; i<GetNRows(); i++)
      for (j=0; j<GetNCols(); j++)
         r(j, i) = (*this)(i, j);
      
   return r;
}

void COMatrix::Resize(int nr, int nc)
{
   assert((nr > 0) && (nc > 0));
   PrivateResize(nr, nc);
}

void COMatrix::Resize(int nr, int nc, double val)
{
   assert((nr > 0) && (nc > 0));
   PrivateResize(nr, nc);
   int i;
   for (i=0; i<nr*nc; i++)
      data_[i] = val;
}

double COMatrix::Determinant(void) const
{
   assert((IsAllocated()) && (GetNRows() == GetNCols()));
/*   
   COMatrix d(GetNRows(), GetNCols());
   COMatrix u, v;
   Svd(u, d, v, "N", "N");
   
   double det = 1;
   for (int i=0; i<d.GetNRows(); i++)
      det *= d(i, i);
   
   return det;
*/
   integer* permidx;
   permidx = new integer[nRows_];
   COMatrix LU(nRows_, nCols_);
   LUDecomp(LU, permidx);
   double det = 1;
   int i;
   for (i=0; i<nRows_; i++)
      det *= LU(i, i);
   // set also the sign
   double sgn=1;
   for (i=0; i<nRows_; i++)
   {
      if (permidx[i]!=(i+1))
         sgn = -sgn;
   }
   det *= sgn;
   delete [] permidx;
   return det;
}

COMatrix COMatrix::Decomp()
{
   assert(IsAllocated() && (nCols_ == nRows_));
   COMatrix v(nCols_, nRows_);
   COMatrix d(nCols_, nRows_);
   COMatrix u;
   Svd(u, d, v, "N", "A");
   int i;
   for (i=0; i<nCols_; i++)
      d(i,i) = sqrt(d(i,i));

   v = d*v.Transpose();
   return v;
}

void COMatrix::LUDecomp(COMatrix& LU, long* permidx) const
{
   // set parameters
   integer m = GetNRows();
   integer n = GetNCols();
   
   doublereal *a;
   LU = (*this);
   a = LU.GetData();
   integer lda = m;
   integer *ipiv;
   ipiv = permidx;
   integer info;

   dgetrf_(&m, &n, a, &lda, ipiv, &info);
}

void COMatrix::RQDecomp(COMatrix& RQ, COVector& Tau) const
{
   integer m = GetNRows();
   integer n = GetNCols();
   doublereal *a;
   RQ = (*this);
   a = RQ.GetData();
   integer lda = m;
   doublereal *tau;
   tau = Tau.GetData();
   doublereal *work;
   integer lwork;
   integer info;
   lwork = m*3;
   work = new doublereal[lwork];

   dgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);

   delete [] work;
}

void COMatrix::RQDecomp(COMatrix& R, COMatrix& Q) const
{
   // assume orthogonal matrix
   if (nRows_ != nCols_)
   {
      COLog("rq no orthogonal\n");
      exit(0);
   }
   R.Resize(nRows_, nRows_);
   Q.Resize(nRows_, nCols_);
   COMatrix RQ(nRows_, nCols_);
   COVector Tau(nRows_);
   RQDecomp(RQ, Tau);

   // retrieve params
   int r, c;
   R.Resize(nRows_, nRows_, 0.0);
   for (r=0; r<nRows_; r++)
   {
      for (c=r; c<nRows_; c++)
      {
         R(r,c) = RQ(r,c);
      }
   }

   int i, j;
   COMatrix Eye(nRows_, nRows_, 0.0);
   for (i=0; i<nRows_; i++)
      Eye(i,i) = 1.0;
   Q = Eye;
   COVector V(nRows_);
   COMatrix Oprod(nRows_, nRows_);
   double* opData;
   opData = Oprod.GetData();

   for (i=0; i<nRows_; i++)
   {
      for (j=0; j<=(i-1); j++)
         V(j) = RQ(i, j);
      V(i) = 1.0;
      for (j=(i+1); j<nRows_; j++)
         V(j) = 0.0;
      COMatrix::OutProduct(V.GetData(), V.GetData(), nRows_, opData);
      Q = Q*(Eye-Tau(i)*Oprod);
   }
}


void COMatrix::LUDecomp(COMatrix& L, COMatrix& U, COMatrix& P) const
{
   L.Resize(nRows_, nCols_, 0.0);
   integer *permidx = new integer[COMin(nRows_, nCols_)];
   LUDecomp(L, permidx);
   int r, c;
   if (nRows_>=nCols_)
   {
      U.Resize(nCols_, nCols_, 0.0);
      for (r=0; r<nCols_; r++)
      {
         for (c=r; c<nCols_; c++)
         {
            U(r, c) = L(r, c);
            L(r, c) = 0;
         }
         L(r,r) = 1.0;
      }
   } else
   {
      U = L;
      L.Resize(nRows_, nRows_, 0.0);
      for (r=0; r<nRows_; r++)
      {
         for (c=(r+1); c<nRows_; c++)
         {
            L(r, c) = U(r, c);
            U(r, c) = 0;
         }
         L(r,r) = 1.0;
      }
   }
   // compute p
   int i;
   P.Resize(nRows_, nRows_, 0.0);
   for (i=0; i<COMin(nRows_, nCols_); i++)
      P(i, permidx[i]) = 1.0;

   delete [] permidx;
}

void COMatrix::SetRotationQuaternion(COVector& quaternion)
{
   assert(quaternion.GetLength() == 4);
   PrivateResize(3,3);
   quaternion.Normalize();
   double q0,q1,q2,q3;
   q0 = quaternion(1);
   q1 = quaternion(2);
   q2 = quaternion(3);
   q3 = -quaternion(0);

   double x2 = q0*q0;
   double y2 = q1*q1;
   double z2 = q2*q2;
   double r2 = q3*q3;

   (*this)(0,0) = r2 + x2 - y2 - z2;		// fill diagonal terms
   (*this)(1,1) = r2 - x2 + y2 - z2;
   (*this)(2,2) = r2 - x2 - y2 + z2;
   double xy = q0 * q1;
   double yz = q1 * q2;
   double zx = q2 * q0;
   double rx = q3 * q0;
   double ry = q3 * q1;
   double rz = q3 * q2;
   (*this)(0,1) = 2 * (xy + rz);			// fill off diagonal terms
   (*this)(0,2) = 2 * (zx - ry);
   (*this)(1,2) = 2 * (yz + rx);
   (*this)(1,0) = 2 * (xy - rz);
   (*this)(2,0) = 2 * (zx + ry);
   (*this)(2,1) = 2 * (yz - rx);

}

void COMatrix::SetRotationQuaternion(double q3, double q0, double q1, double q2)
{
   q3 = -q3;
   PrivateResize(3,3);

   double x2 = q0*q0+q1*q1+q2*q2+q3*q3;
   x2 = sqrt(x2);
   q0 /= x2;
   q1 /= x2;
   q2 /= x2;
   q3 /= x2;

   x2 = q0*q0;
   double y2 = q1*q1;
   double z2 = q2*q2;
   double r2 = q3*q3;

   (*this)(0,0) = r2 + x2 - y2 - z2;		// fill diagonal terms
   (*this)(1,1) = r2 - x2 + y2 - z2;
   (*this)(2,2) = r2 - x2 - y2 + z2;
   double xy = q0 * q1;
   double yz = q1 * q2;
   double zx = q2 * q0;
   double rx = q3 * q0;
   double ry = q3 * q1;
   double rz = q3 * q2;
   (*this)(0,1) = 2 * (xy + rz);			// fill off diagonal terms
   (*this)(0,2) = 2 * (zx - ry);
   (*this)(1,2) = 2 * (yz + rx);
   (*this)(1,0) = 2 * (xy - rz);
   (*this)(2,0) = 2 * (zx + ry);
   (*this)(2,1) = 2 * (yz - rx);

}


void COMatrix::GetRotationQuaternion(COVector& quaternion)
{
   assert((nRows_ == 3) && (nCols_==3));
   // orthogonal matrix ??
   
   double d0 = (*this)(0,0), d1 = (*this)(1,1), d2 = (*this)(2,2);
   double xx = 1.0 + d0 - d1 - d2;		// from the diagonal of (*this)ation
   double yy = 1.0 - d0 + d1 - d2;		// matrix, find the terms in
   double zz = 1.0 - d0 - d1 + d2;		// each Quaternion compoment
   double rr = 1.0 + d0 + d1 + d2;
   
   double max = rr;				// find the maximum of all
   if (xx > max) max = xx;				// diagonal terms.
   if (yy > max) max = yy;
   if (zz > max) max = zz;
   
   if (rr == max) {
      double r4 = sqrt(rr * 4.0);
      quaternion(1) = ((*this)(1,2) - (*this)(2,1)) / r4;	// find other components from
      quaternion(2) = ((*this)(2,0) - (*this)(0,2)) / r4;	// off diagonal terms of
      quaternion(3) = ((*this)(0,1) - (*this)(1,0)) / r4;	// (*this)ation matrix.
      quaternion(0) = -r4 / 4.0;
   } else if (xx == max) {
      double x4 = sqrt(xx * 4.0);
      quaternion(1) = x4 / 4.0;
      quaternion(2) = ((*this)(0,1) + (*this)(1,0)) / x4;
      quaternion(3) = ((*this)(0,2) + (*this)(2,0)) / x4;
      quaternion(0) = -((*this)(1,2) - (*this)(2,1)) / x4;
   } else if (yy == max) {
      double y4 = sqrt(yy * 4.0);
      quaternion(1) = ((*this)(0,1) + (*this)(1,0)) / y4;
      quaternion(2) =  y4 / 4.0;
      quaternion(3) = ((*this)(1,2) + (*this)(2,1)) / y4;
      quaternion(0) = -((*this)(2,0) - (*this)(0,2)) / y4;
   } else {
      double z4 = sqrt(zz * 4.0);
      quaternion(1) = ((*this)(0,2) + (*this)(2,0)) / z4;
      quaternion(2) = ((*this)(1,2) + (*this)(2,1)) / z4;
      quaternion(3) =  z4 / 4.0;
      quaternion(0) = -((*this)(0,1) - (*this)(1,0)) / z4;
   }
}

void COMatrix::GetRotationQuaternion(double* q0, double* q1, double* q2, double* q3)
{
   assert((nRows_ == 3) && (nCols_==3));
   // orthogonal matrix ??
   
   double d0 = (*this)(0,0), d1 = (*this)(1,1), d2 = (*this)(2,2);
   double xx = 1.0 + d0 - d1 - d2;		// from the diagonal of (*this)ation
   double yy = 1.0 - d0 + d1 - d2;		// matrix, find the terms in
   double zz = 1.0 - d0 - d1 + d2;		// each Quaternion compoment
   double rr = 1.0 + d0 + d1 + d2;
   
   double max = rr;				// find the maximum of all
   if (xx > max) max = xx;				// diagonal terms.
   if (yy > max) max = yy;
   if (zz > max) max = zz;
   
   if (rr == max) {
      double r4 = sqrt(rr * 4.0);
      *q1 = ((*this)(1,2) - (*this)(2,1)) / r4;	// find other components from
      *q2 = ((*this)(2,0) - (*this)(0,2)) / r4;	// off diagonal terms of
      *q3 = ((*this)(0,1) - (*this)(1,0)) / r4;	// (*this)ation matrix.
      *q0 = -r4 / 4.0;
   } else if (xx == max) {
      double x4 = sqrt(xx * 4.0);
      *q1 = x4 / 4.0;
      *q2 = ((*this)(0,1) + (*this)(1,0)) / x4;
      *q3 = ((*this)(0,2) + (*this)(2,0)) / x4;
      *q0 = -((*this)(1,2) - (*this)(2,1)) / x4;
   } else if (yy == max) {
      double y4 = sqrt(yy * 4.0);
      *q1 = ((*this)(0,1) + (*this)(1,0)) / y4;
      *q2 =  y4 / 4.0;
      *q3 = ((*this)(1,2) + (*this)(2,1)) / y4;
      *q0 = -((*this)(2,0) - (*this)(0,2)) / y4;
   } else {
      double z4 = sqrt(zz * 4.0);
      *q1 = ((*this)(0,2) + (*this)(2,0)) / z4;
      *q2 = ((*this)(1,2) + (*this)(2,1)) / z4;
      *q3 =  z4 / 4.0;
      *q0 = -((*this)(0,1) - (*this)(1,0)) / z4;
   }
//   max = (*q0)*(*q0) + (*q1)*(*q1) +(*q2)*(*q2) +(*q3)*(*q3);
//   max = sqrt(max);
//   *q0 /= max;
//   *q1 /= max;
//   *q2 /= max;
//   *q3 /= max;
}


void COMatrix::SetRotationAngles(double a, double b, double g)
{
   PrivateResize(3,3);
   data_[0] = cos(b)*cos(g);
   data_[1] = sin(a)*sin(b)*cos(g)+cos(a)*sin(g);
   data_[2] = -cos(a)*sin(b)*cos(g)+sin(a)*sin(g);
   data_[3] = -cos(b)*sin(g);
   data_[4] = -sin(a)*sin(b)*sin(g)+cos(a)*cos(g);
   data_[5] = cos(a)*sin(b)*sin(g)+sin(a)*cos(g);
   data_[6] = sin(b);
   data_[7] = -sin(a)*cos(b);
   data_[8] = cos(a)*cos(b);
}

void COMatrix::GetRotationAngles(double* alpha, double* beta, double* gamma)
{
   assert((nRows_ == 3) && (nCols_==3));

   double b1, b2, g, a1, a2;
   b1 = asin(data_[6]);
   if (b1>=0)
      b2 = PI - b1;
   else
      b2 = -PI - b1;
   COMatrix trot(3,3);
   double mn, tmn;

   double cb;
   cb = cos(b1);
   if (fabs(cb) > 0.0005)
   {
      // normal case
      
      // first for b1
      g = atan2(-data_[3]/cb, data_[0]/cb);
      a1 = asin(-data_[7]/cb);
      if (a1 >= 0)
         a2 = PI - a1;
      else
         a2 = -PI - a1;

      // test for a1,a2
      trot.SetRotationAngles(a1, b1, g);
      mn = ((*this)-trot).NormFrobenius();
      (*alpha) = a1;
      (*beta) = b1;
      (*gamma) = g;

      trot.SetRotationAngles(a2, b1, g);
      tmn = ((*this)-trot).NormFrobenius();
      if (tmn<mn)
      {
         mn=tmn;
         (*alpha) = a2;
         (*beta) = b1;
         (*gamma) = g;
      }

      // next for b2
      cb = cos(b2);
      g = atan2(-data_[3]/cb, data_[0]/cb);
      a1 = asin(-data_[7]/cb);
      if (a1 >= 0)
         a2 = a1+PI/2;
      else
         a2 = a1-PI/2;

      // test for a1,a2
      trot.SetRotationAngles(a1, b2, g);
      tmn = ((*this)-trot).NormFrobenius();
      if (tmn<mn)
      {
         mn=tmn;
         (*alpha) = a1;
         (*beta) = b2;
         (*gamma) = g;
      }

      trot.SetRotationAngles(a2, b2, g);
      tmn = ((*this)-trot).NormFrobenius();
      if (tmn<mn)
      {
         mn=tmn;
         (*alpha) = a2;
         (*beta) = b2;
         (*gamma) = g;
      }
   } else
   {
      // not normal case, assume b = +- pi/2 and one value (b2 ~= b1)
      // then gamma = 0;
      g = 0;
      a1 = atan2(data_[1]/data_[6], -data_[2]/data_[6]);
      (*alpha) = a1;
      (*beta) = b1;
      (*gamma) = g;
   }
}

void COMatrix::Svd(COMatrix& u, COMatrix& d, COMatrix& v, char* uType, char* vType) const
{
   assert(IsAllocated());
   int mindim;

   integer svdldu = (strcmp(uType,"A")==0) ? GetNRows() : 1;
   integer svdlda = GetNRows();
   integer svdldvt = (strcmp(vType,"A")==0) ? GetNCols() : 1;
   integer svdlwork = GetNElements()+25;
   integer svdinfo;
   integer svdm = GetNRows();
   integer svdn = GetNCols();
   double* svdu = (strcmp(uType,"A")==0) ? u.GetData() : 0;
   double* svdvt = (strcmp(vType,"A")==0) ? v.GetData() : 0;
   double* svdwork = new double[svdlwork];
   mindim = (svdm<svdn) ? svdm : svdn;//min(svdm,svdn);
   double* svdd = new double[mindim];
   COMatrix inmat(*this);
   double* svdinmat = inmat.GetData();

   dgesvd_(uType, vType, &svdm, &svdn, svdinmat, &svdlda, svdd, svdu, &svdldu,
      svdvt, &svdldvt, svdwork, &svdlwork, &svdinfo);

   double* dtmp = d.GetData();
   for (int i=0; i<d.GetNElements(); i++)
      dtmp[i] = 0;
   for (int i=0; i<mindim; i++)
      d(i,i) = svdd[i];
   if (strcmp(vType,"A")==0)
      v=v.Transpose();
   delete [] svdd;
   delete [] svdwork;
}

void COMatrix::SvdSpec(int nr, int nc, COMatrix& d) const
{
   assert(IsAllocated());
   int mindim;

   integer svdldu = 1;
   integer svdlda = nr;
   integer svdldvt = 1;
   integer svdlwork = nr*nc+25;
   integer svdinfo;
   integer svdm = nr;
   integer svdn = nc;
   double* svdu = 0;
   double* svdvt = 0;
   double* svdwork = new double[svdlwork];
   mindim = (svdm < svdn) ? svdm : svdn; //min(svdm,svdn);
   double* svdd = new double[mindim];
   COMatrix inmat(*this);
   double* svdinmat = inmat.GetData();

   dgesvd_("N", "N", &svdm, &svdn, svdinmat, &svdlda, svdd, svdu, &svdldu,
      svdvt, &svdldvt, svdwork, &svdlwork, &svdinfo);

   double* dtmp = d.GetData();
   for (int i=0; i<(nr*nc); i++)
      dtmp[i] = 0;
   for (int i=0; i<mindim; i++)
      d(i,i) = svdd[i];
   delete [] svdd;
   delete [] svdwork;
}

void COMatrix::Write(FILE* fd, int mtype) const
{
   assert(IsAllocated());

   if (mtype == CO_MATRIX_ASCII)
   {
      fprintf(fd,"COMatrix\n%d %d\nASCII\n", GetNRows(), GetNCols());
      int i, j;
      for (i=0; i<GetNRows(); i++)
      {
         for (j=0; j<GetNCols(); j++)
         {
            fprintf(fd,"%g ", (*this)(i,j));
         }
         fprintf(fd, "\n");
      }
   }

   if (mtype == CO_MATRIX_BIN64)
   {
      fprintf(fd,"COMatrix\n%d %d\nBIN64\n", GetNRows(), GetNCols());
      fwrite(data_, sizeof(double), GetNElements(), fd);
   }

   if (mtype == CO_MATRIX_BIN32)
   {
      fprintf(fd,"COMatrix\n%d %d\nBIN32\n", GetNRows(), GetNCols());
      int ne = GetNElements();
      float* temp;
      temp = new float[ne];
      int i;
      for (i=0; i<ne; i++)
         temp[i] = (float) data_[i];
      fwrite(temp, sizeof(float), ne, fd);
      delete [] temp;
   }
}

/*
void COMatrix::LogFile() const
{
   assert(IsAllocated());
   COLogFile("COMatrix\n%d %d\nASCII\n", GetNRows(), GetNCols());
   int i, j;
   for (i=0; i<GetNRows(); i++)
   {
      for (j=0; j<GetNCols(); j++)
      {
         COLogFile("%9.6f ", (*this)(i,j));
      }
      COLogFile("\n");
   }
}
*/
void COMatrix::Write(const char* name, int mtype) const
{
   FILE* fd;
   fd = fopen(name, "wb");
   if (fd == NULL)
   {
      printf("COMatrix::write: cannot write file %s\n",name);
      exit(-1);
   }
   Write(fd,mtype);
   fclose(fd);
}

void COMatrix::Read(FILE* fd)
{
   char s[100];
   fscanf(fd, "%s", s);
   if (strcmp(s, "COMatrix") != 0) {
      printf("COMatrix::Read: file doesn't start with MtMatrix.\n");
      exit(-1);
   }
   
   int nr, nc;
   fscanf(fd, "%d %d", &nr, &nc);
   assert((nr >= 0) && (nc >= 0));
   Resize(nr, nc);
   
   fscanf(fd, "%s", s);

   if (strcmp(s, "ASCII") == 0)
   {
      // read ascii   
      int i, j;
      float val;
      for (i=0; i<GetNRows(); i++)
      {
         for (j=0; j<GetNCols(); j++)
         {
            fscanf(fd,"%g", &val);
            data_[i+j*nRows_] = val;
         }
      }
   } else if(strcmp(s, "BIN64") == 0)
   {
      // read double directly
      unsigned char ct;
      fread(&ct, 1, 1, fd);
      while(ct != 10)
         fread(&ct, 1, 1, fd);
      fread(data_, sizeof(double), nr*nc, fd);
   } else if(strcmp(s, "BIN32") == 0)
   {
      unsigned char ct;
      fread(&ct, 1, 1, fd);
      while(ct != 10)
         fread(&ct, 1, 1, fd);

      // read from float
      int i;
      float* temp;
      temp = new float[nr*nc];
      fread(temp, sizeof(float), nr*nc, fd);
      for (i=0; i<nr*nc; i++)
         data_[i] = temp[i];
      delete [] temp;
   } else
   {
      printf("COMatrix::Read: unknown matrix format.\n");
      exit(-1);
   }

}

void COMatrix::Read(const char* name)
{
   FILE* fd;
   fd = fopen(name,"rb");
   if (fd == NULL)
   {
      printf("COMatrix::Read: cannot read file %s",name);
      exit(-1);
   }
   Read(fd);
   fclose(fd);
}

COMatrix OuterProduct(const COVector& v1, const COVector& v2)
{
   assert(v1.IsAllocated() && v2.IsAllocated());
   int nRows = v1.GetLength();
   int nCols = v2.GetLength();
   int r, c;
   COMatrix m(nRows, nCols);
   for (r = 0; r<nRows; r++)
   {
      for (c = 0; c<nCols; c++)
      {
         m(r, c) = v1(r) * v2(c);
      }
   }

   return m;
}

COMatrix COMatrix::CorrMat()
{
   int i,j,k;
   assert(IsAllocated());
   // compute M'*M
   COMatrix MtM(nCols_, nCols_);

   double *dki, *dkj, dtemp;
   for (i=0; i<nCols_; i++)
   {
      dki = data_+i*nRows_;
      for (j=i; j<nCols_; j++)
      {
         dkj = data_+j*nRows_;
         dtemp = 0;
         for (k=0; k<nRows_; k++)
         {
            dtemp += dki[k]*dkj[k];
         }
         MtM(i,j) = dtemp;
         MtM(j,i) = dtemp;
      }
   }    
   return MtM;
}
