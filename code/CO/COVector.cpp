//
//  Created by Oncel Tuzel.
//  Copyright © 2004 Oncel Tuzel. All rights reserved.
//


#include "COVector.h"



/****************************************************************************/

void COVector::Clear()
{
	memset(data_,0,sizeof(double)*length_);
}

void COVector::PrivateAllocateData(int length)
{
  assert(length > 0);

  length_ = length;
  data_ = new double[length_];
  assert(data_ != 0);
  isAllocated_ = 1;
}

void COVector::PrivateResize(int length)
{
  if (!IsAllocated())
    PrivateAllocateData(length);
  else if(length_ != length)
  {
    CleanData();
    PrivateAllocateData(length);
  }
}

void COVector::PrivateCopyToThis(const COVector& v)
{
  PrivateResize(v.GetLength());

  int i;
  double* indata;
  indata = v.GetData();
  for (i=0; i<length_; i++)
     data_[i] = indata[i];
}

void COVector::CleanData(void)
{
   delete [] data_;
   isAllocated_ = 0;
}

COVector::COVector(void)
{
   isAllocated_ = 0;
   data_ = 0;
}

COVector::COVector(int length)
{
   assert(length > 0);
   PrivateAllocateData(length);
}

COVector::COVector(int length, double val)
{
   assert(length > 0);
   PrivateAllocateData(length);
   
   int i;
   for (i=0; i<length_; i++)
      data_[i] = val;
}

COVector::COVector(int length, double* data)
{
   assert(length > 0);
   PrivateAllocateData(length);
   int i;
   for (i=0; i<length_; i++)
      data_[i] = data[i];
}

void COVector::SetSameVector(double* data)
{
   int i;
   for (i=0; i<length_; i++)
      data_[i] = data[i];
}

void COVector::GetSameVector(double* data)
{
   int i;
   for (i=0; i<length_; i++)
      data[i] = data_[i];
}

COVector::COVector(const char* name)
{
   isAllocated_ = 0;
   Read(name);
}

COVector::COVector(const COVector& v)
{
   isAllocated_ = 0;
   if (!v.IsAllocated())
      return;
   PrivateCopyToThis(v);
}

const COVector& COVector::operator=(const COVector& v)
{
   if (this == &v)
      return *this;
   
   if (!v.IsAllocated())
   {
      isAllocated_ = 0;
      return *this;
   }
   
   PrivateCopyToThis(v);
   return *this;
}

COVector::~COVector(void) 
{
   if (IsAllocated())
      CleanData();
}

int COVector::operator==(const COVector& v) const
{
   assert((IsAllocated()) && (v.IsAllocated()));
   
   if (GetLength() != v.GetLength())
      return 0;
   int i;
   double* datav;
   datav = v.GetData();
   for (i=0; i<length_; i++)
   {
      if (data_[i] != datav[i])
      {
         return 0;
      }
   }
   return 1;
}

int COVector::operator!=(const COVector& v) const
{
   return !(*this == v);
}

COVector COVector::operator+(double d) const
{
   assert(IsAllocated());
   
   COVector r(GetLength());
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]+d;      
   return r;
}

COVector COVector::operator+(const COVector& v) const
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]+datav[i];
   return r;
}

COVector operator+(double d, const COVector& v)
{
   assert(v.IsAllocated());
   
   COVector r(v.GetLength());
   
   int i;
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   for (i=0; i<v.GetLength(); i++)
      datar[i] = d+datav[i];      
   return r;
}

COVector COVector::operator+=(const COVector& v)
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   int i;
   double* datav;
   datav = v.GetData();
   for (i=0; i<GetLength(); i++)
      data_[i] += datav[i];      
   return *this;
}

COVector COVector::operator+=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] +=d;
   return *this;
}

COVector COVector::operator-(double d) const
{
   assert(IsAllocated());
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]-d;
      
   return r;
}

COVector COVector::operator-(const COVector& v) const
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]-datav[i];
   return r;
}

COVector operator-(double d, const COVector& v)
{
   assert(v.IsAllocated());
   
   COVector r(v.GetLength());
   
   int i;
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   for (i=0; i<v.GetLength(); i++)
      datar[i] = d-datav[i];
   return r;
}

COVector COVector::operator-(void) const
{
   assert(IsAllocated());
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = -data_[i];
   return r;
}

COVector COVector::operator-=(const COVector& v)
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   int i;
   double* datav;
   datav = v.GetData();
   for (i=0; i <GetLength(); i++)
      data_[i] -= datav[i];
   return *this;
}

COVector COVector::operator-=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] -=d;
   return *this;
}

COVector COVector::operator*(double d) const
{
   assert(IsAllocated());
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]*d;
   return r;
}

COVector COVector::operator*(const COVector& v) const
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   COVector r(GetLength());
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   
   int i;
   for (i=0; i<GetLength(); i++)
   {
      datar[i] = data_[i]*datav[i];
   }
   
   return r;
}


COVector operator*(double d, const COVector& v)
{
   assert(v.IsAllocated());
   
   COVector r(v.GetLength());
   
   int i;
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   for (i=0; i<v.GetLength(); i++)
      datar[i] = d*datav[i];
   return r;
}

COVector COVector::operator*=(const COVector& v)
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   double* datav;
   datav = v.GetData();

   int i;
   for (i=0; i<GetLength();i++)
   {
      data_[i] *= datav[i];
   }
   return *this;
}

COVector COVector::operator*=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] *= d;
   return *this;
}


COVector COVector::operator/(const COVector& v) const
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   COVector r(GetLength());
   double* datar;
   double* datav;
   datar = r.GetData();
   datav = v.GetData();
   
   int i;
   for (i=0; i<GetLength(); i++)
   {
      datar[i] = data_[i]/datav[i];
   }
   
   return r;
}

COVector COVector::operator/(double d) const
{
   assert(IsAllocated());
   
   COVector r(GetLength());
   
   int i;
   double* datar;
   datar = r.GetData();
   for (i=0; i<GetLength(); i++)
      datar[i] = data_[i]/d;
   return r;
}

COVector COVector::operator/=(const COVector& v)
{
   assert((IsAllocated()) && (v.IsAllocated()) &&
      (GetLength() == v.GetLength()));
   
   double* datav;
   datav = v.GetData();

   int i;
   for (i=0; i<GetLength();i++)
   {
      data_[i] /= datav[i];
   }
   return *this;
}

COVector COVector::operator/=(double d)
{
   assert(IsAllocated());
   
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] /= d;
   return *this;
}

COVector COVector::GetSubVector(int idxs, int idxe) const
{
   assert((IsAllocated()) &&
      (idxs >= 0) && (idxs < GetLength()) &&
      (idxe >= 0) && (idxe < GetLength()) && (idxe>=idxs));
   
   COVector r(idxe-idxs+1);
   double* datar;
   datar = r.GetData();
   
   int i;
   for (i=idxs; i<=idxe; i++)
      datar[i-idxs] = data_[i];
      
   return r;
}

void COVector::SetSubVector(int idxs, COVector& dv)
{
   assert(dv.IsAllocated());
   int idxe;
   idxe = idxs+dv.GetLength()-1;
   assert((IsAllocated()) &&
      (idxs >= 0) && (idxs < GetLength()) &&
      (idxe >= 0) && (idxe < GetLength()));

   double* datav;
   datav = dv.GetData();
   int i;
   for (i=idxs; i<=idxe; i++)
      data_[i] = datav[i-idxs];
}

void COVector::Normalize(void)
{
   assert(IsAllocated());
   
   double norm = GetNorm();
   if (norm < 0.00000000001)
      return;
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] /= norm;
      
}


void COVector::Normalize1Norm(void)
{
   assert(IsAllocated());
   
   double norm = GetNorm1();
   if (norm < 0.00000000001)
      return;
   int i;
   for (i=0; i<GetLength(); i++)
      data_[i] /= norm;
      
}

double COVector::GetNorm(void)
{
   assert(IsAllocated());
   double norm = 0;
   int i;
   for (i=0; i<length_; i++)
      norm += data_[i]*data_[i];
   return (double) sqrt(norm);
}

double COVector::GetNorm2(void)
{
   assert(IsAllocated());
   double norm = 0;
   int i;
   for (i=0; i<length_; i++)
      norm += data_[i]*data_[i];
   return norm;
}

double COVector::GetNorm1(void)
{
   assert(IsAllocated());
   double norm = 0;
   int i;
   for (i=0; i<length_; i++)
      norm += fabs(data_[i]);
   return norm;
}

double COVector::InnerProd(const COVector& v) const
{
   assert(IsAllocated() && v.IsAllocated() && (GetLength() == v.GetLength()));
   int i;
   double* datav;
   datav = v.GetData();
   double ip = 0;
   for (i=0; i<length_; i++)
      ip += data_[i]*datav[i];
   return ip;
}

double COVector::Cross2D(const COVector& v) const
{
   assert(IsAllocated() && v.IsAllocated() && (GetLength() == v.GetLength()));

   return data_[0] * v(1) - v(0)*data_[1];
}


double COVector::Mean() const
{
   int i;
   double mean = 0;
   for (i=0; i<length_; i++)
      mean += data_[i];
   mean /= length_;
   return mean;
}

double COVector::Variance() const
{
   int i;
   double mean = 0;
   for (i=0; i<length_; i++)
      mean += data_[i];
   mean /= length_;
   double variance = 0;
   for (i=0; i<length_; i++)
      variance += (data_[i]-mean)*(data_[i]-mean);
   variance /= length_;
   return variance;
}

double COVector::Variance(double mean) const
{
   int i;
   double variance = 0;
   for (i=0; i<length_; i++)
      variance += (data_[i]-mean)*(data_[i]-mean);
   variance /= length_;
   return variance;
}

void COVector::Resize(int length)
{
   assert(length > 0);
   PrivateResize(length);
}

void COVector::Resize(int length, double val)
{
   assert(length > 0);
   PrivateResize(length);
   int i;
   for (i=0; i<length_; i++)
      data_[i] = val;
}

void COVector::Write(FILE* fd) const
{
   assert(IsAllocated());
   fprintf(fd,"COVector\n%d\nASCII\n", GetLength());
   int i;
   for (i=0; i<GetLength(); i++)
   {
      fprintf(fd,"%g ", data_[i]);
   }
   fprintf(fd, "\n");
}

void COVector::Write(const char* name) const
{
   FILE* fd;
   fd = fopen(name, "w");
   if (fd == NULL)
   {
      printf("COVector::write: cannot write file %s\n",name);
      exit(-1);
   }
   Write(fd);
   fclose(fd);
}

void COVector::Read(FILE* fd)
{
   char s[100];
   fscanf(fd, "%s", s);
   if (strcmp(s, "COVector") != 0) {
      printf("COVector::Read: file doesn't start with COVector.\n");
      exit(-1);
   }
   
   int length;
   fscanf(fd, "%d", &length);
   assert(length > 0);
   Resize(length);
   
   fscanf(fd, "%s", s);
   if (strcmp(s, "ASCII") != 0) {
      printf("COVector::Read: file doesn't contain the ASCII keyword.\n");
      exit(-1);
   }
   
   int i;
   float idf;
   for (i=0; i<GetLength(); i++)
   {
      fscanf(fd,"%g", &idf);
      data_[i] = idf;
   }
}

void COVector::Read(const char* name)
{
   FILE* fd;
   fd = fopen(name,"rb");
   if (fd == NULL)
   {
      printf("COVector::Read: cannot read file %s",name);
      exit(-1);
   }
   Read(fd);
   fclose(fd);
}


