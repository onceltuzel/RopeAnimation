


#include <math.h>
#include <stdio.h>
#include <string.h>

// !!! only MSV
//#include <sys/timeb.h>
#include <time.h>
#include <stdarg.h>

#include "CODefaults.h"


void COLog(const char *varStr, ...)
{
//   return;

	//obtain argument list using ANSI standard...
	va_list	argList;
	va_start(argList, varStr);

	//print the output string to stderr using
	vfprintf(stdout, varStr, argList);
	va_end(argList);

	//done.
	return;
}


double COSign(double x)
{
   if (x>=0)
      return 1.0;
   return -1.0;
}

inline int COMin(int a_in, int b_in)
{
   if (a_in<b_in)
      return a_in;
   return b_in;
}

inline int COMax(int a_in, int b_in)
{
   if (a_in>b_in)
      return a_in;
   return b_in;
}

inline double COMin(double a_in, double b_in)
{
   if (a_in<b_in)
      return a_in;
   return b_in;
}

inline int COMin(double* v_in, int nv_in)
{
   int idx_in,i_in;
   i_in = 0;
   for (idx_in=1; idx_in<nv_in; idx_in++)
      if (v_in[idx_in]<v_in[i_in])
         i_in=idx_in;
   return i_in;
}

inline int COIMin(double a_in, double b_in, double c_in)
{
   int idx_in;
   if (a_in <= b_in) 
   {
      if (a_in <= c_in)
         idx_in = 1;
      else
         idx_in = 3;
   } else if (b_in <= c_in) idx_in = 2;
   else idx_in = 3;
   return idx_in;
}

inline double COMin3(double a_in, double b_in, double c_in)
{  
   double min_in;
   if (a_in <= b_in) 
   {
      if (a_in <= c_in)
         min_in = a_in;
      else
         min_in = c_in;
   } else if (b_in <= c_in) min_in =  b_in;
   else min_in = c_in;
   return min_in;
}

inline double COMax(double a_in, double b_in)
{
   if (a_in>b_in)
      return a_in;
   return b_in;
}


double COFactorial(double num)
{
   if (num==0 || num==1)
      return 1;
   
   return (num * COFactorial(num - 1));
} 

// return 0 - error
// return 1 - one real solution
// return 3 - three real solutions
int COSolveCubic(double a, double b, double c, double d, double& s1, double& s2, double& s3)
{
   double p, q;
   double r, s, t, z;

   // convert to canonical form
   r = b/a;
   s = c/a;
   t = d/a;
   p = s-(r*r)/3.0;
   q = (2*r*r*r)/27.0-(r*s)/3.0+t;

   double D, phi, R;
   R = COSign(q)*sqrt(fabs(p)/3.0); 
   D = pow(p/3.0,3)+pow(q/2,2);

   if (p<0)
   {
      if (D<=0)
      {
         phi = acos(q/(2*R*R*R));
         s1 = -2*R*cos(phi/3)-r/3;
         s2 = -2*R*cos(phi/3+2*PI/3)-r/3;
         s3 = -2*R*cos(phi/3+4*PI/3)-r/3;
         return 3;
      }
      else
      {
         z = q/(2*R*R*R);
         phi = log(z+sqrt(z*z-1));
         s1 = -2*R*cosh(phi/3)-r/3;
         return 1;
      }
   }
   else
   {
      z = q/(2*R*R*R);
      phi = log(z+sqrt(z*z+1));
      s1 = -2*R*sinh(phi/3)-r/3;
      return 1;
   }
   return 0;
}

inline int CORound(double in_x)
{
  return int(floor(in_x + 0.5));
}
inline int CORoundSign(double in_x)
{
   if (in_x>=0)
   {
      return ((int) (in_x + 0.5));
   }
   else
   {
      return ((int) (in_x - 0.5));
   }
}

inline double COEuclidianDistance2D(double* v1, double* v2)
{
   return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]));
}

inline double COEuclidianDistance2D2(double* in_v1, double* in_v2)
{
   return ((in_v1[0]-in_v2[0])*(in_v1[0]-in_v2[0])+(in_v1[1]-in_v2[1])*(in_v1[1]-in_v2[1]));
}


void COSort(double* ra, int nVec)
{
   unsigned long n, l, ir, i, j;
   n = nVec;
   double rra;
   
   if (n<2)
      return;
   l = (n>>1)+1;
   ir = n;
   for (;;)
   {
      if (l>1)
      {
         rra = ra[(--l)-1];
      }
      else
      {
         rra = ra[ir-1];
         ra[ir-1] = ra[1-1];
         if (--ir==1)
         {
            ra[1-1] = rra;
            break;
         }
      }
      i = l;
      j = l+l;
      while (j<=ir)
      {
         if (j<ir && ra[j-1]<ra[j+1-1])
            j++;
         if (rra<ra[j-1])
         {
            ra[i-1] = ra[j-1];
            i = j;
            j <<= 1;
         }
         else
            j = ir+1;
      }
      ra[i-1] = rra;
   }

}

void COISort(double* ra, int nVec, int* ira)
{
   unsigned long n, l, ir, i, j;
   n = nVec;
   double rra;
   int irra;
   
   if (n<2)
      return;
   l = (n>>1)+1;
   ir = n;
   for (;;)
   {
      if (l>1)
      {
         irra = ira[(--l)-1];
         rra = ra[l-1];
      }
      else
      {
         irra = ira[ir-1];
         rra = ra[ir-1];

         ira[ir-1] = ira[1-1];
         ra[ir-1] = ra[1-1];

         if (--ir==1)
         {
            ira[1-1] = irra;
            ra[1-1] = rra;
            break;
         }
      }
      i = l;
      j = l+l;
      while (j<=ir)
      {
         if (j<ir && ra[j-1]<ra[j+1-1])
            j++;
         if (rra<ra[j-1])
         {
            ira[i-1] = ira[j-1];
            ra[i-1] = ra[j-1];

            i = j;
            j <<= 1;
         }
         else
            j = ir+1;
      }
      ira[i-1] = irra;
      ra[i-1] = rra;
   }

}

// rank in 0-1 range, 0 min, 1 max
// inplace sort vec
double COMedian(double* vec, int nVec, double rank)
{

   COSort(vec, nVec);
   int krank = int(floor(rank*nVec));
   if (rank == 1)
      krank = nVec-1;
   return vec[krank];
}

double COMedianToSigmaGaussian(double med)
{
  return med * 1.482602219;
}

int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
    int cols, char *comment, int maxval)
{
   FILE *fp;

   /***************************************************************************
   * Open the output image file for writing if a filename was given. If no
   * filename was provided, set fp to write to standard output.
   ***************************************************************************/
   if(outfilename == NULL) fp = stdout;
   else{
      if((fp = fopen(outfilename, "wb")) == NULL){
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
         return(0);
      }
   }

   /***************************************************************************
   * Write the header information to the PGM file.
   ***************************************************************************/
   fprintf(fp, "P5\n%d %d\n", cols, rows);
   if(comment != NULL)
      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   /***************************************************************************
   * Write the image data to the file.
   ***************************************************************************/
   if((unsigned int) rows != fwrite(image, cols, rows, fp)){
      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
      if(fp != stdout) fclose(fp);
      return(0);
   }

   if(fp != stdout) fclose(fp);
   return(1);
}

int read_pgm_image(char *infilename, unsigned char *image, int *rows,
    int *cols)
{
   FILE *fp;
   char buf[71];

   /***************************************************************************
   * Open the input image file for reading if a filename was given. If no
   * filename was provided, set fp to read from standard input.
   ***************************************************************************/
   if(infilename == NULL) fp = stdin;
   else{
      if((fp = fopen(infilename, "rb")) == NULL){
         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
            infilename);
         return(0);
      }
   }

   /***************************************************************************
   * Verify that the image is in PGM format, read in the number of columns
   * and rows in the image and scan past all of the header information.
   ***************************************************************************/
   fgets(buf, 70, fp);
   if(strncmp(buf,"P5",2) != 0){
      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
      fprintf(stderr, "read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   sscanf(buf, "%d %d", cols, rows);
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

   /***************************************************************************
   * Allocate memory to store the image then read the image from the file.
   ***************************************************************************/
   /*if(((*image) = new unsigned char[(*rows)*(*cols)]) == NULL){
      fprintf(stderr, "Memory allocation failure in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
      }*/
   if((unsigned int) (*rows) != fread(image, (*cols), (*rows), fp)){
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }

   if(fp != stdin) fclose(fp);
   return(1);
}
/*
int read_ppm_image(char *infilename, unsigned char **image_red, 
    unsigned char **image_grn, unsigned char **image_blu, int *rows,
    int *cols)
{
   FILE *fp;
   char buf[71];
   int p, size;

   if(infilename == NULL) fp = stdin;
   else{
      if((fp = fopen(infilename, "r")) == NULL){
         fprintf(stderr, "Error reading the file %s in read_ppm_image().\n",
            infilename);
         return(0);
      }
   }

   fgets(buf, 70, fp);
   if(strncmp(buf,"P6",2) != 0){
      fprintf(stderr, "The file %s is not in PPM format in ", infilename);
      fprintf(stderr, "read_ppm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  // skip all comment lines
   sscanf(buf, "%d %d", cols, rows);
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  // skip all comment lines

   if(((*image_red) = (unsigned char *) malloc((*rows)*(*cols))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_ppm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   if(((*image_grn) = (unsigned char *) malloc((*rows)*(*cols))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_ppm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   if(((*image_blu) = (unsigned char *) malloc((*rows)*(*cols))) == NULL){
      fprintf(stderr, "Memory allocation failure in read_ppm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }

   size = (*rows)*(*cols);
   for(p=0;p<size;p++){
      (*image_red)[p] = (unsigned char)fgetc(fp);
      (*image_grn)[p] = (unsigned char)fgetc(fp);
      (*image_blu)[p] = (unsigned char)fgetc(fp);
   }

   if(fp != stdin) fclose(fp);
   return(1);
}
*/
int write_ppm_image(char *outfilename, unsigned char *image_red,
    unsigned char *image_grn, unsigned char *image_blu, int rows,
    int cols, char *comment, int maxval)
{
   FILE *fp;
   long size, p;

   if(outfilename == NULL) fp = stdout;
   else{
      if((fp = fopen(outfilename, "wb")) == NULL){
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
         return(0);
      }
   }

   fprintf(fp, "P6\n%d %d\n", cols, rows);
   if(comment != NULL)
      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   size = (long)rows * (long)cols;
   for(p=0;p<size;p++){      // Write the image in pixel interleaved format.
      fputc(image_red[p], fp);
      fputc(image_grn[p], fp);
      fputc(image_blu[p], fp);
   }

   if(fp != stdout) fclose(fp);
   return(1);
}

int write_ppm_image2(char *outfilename, unsigned char *image_all,
    int rows, int cols, char *comment, int maxval)
{
   FILE *fp;
   long size, p;

   if(outfilename == NULL) fp = stdout;
   else{
      if((fp = fopen(outfilename, "wb")) == NULL){
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
         return(0);
      }
   }

   fprintf(fp, "P6\n%d %d\n", cols, rows);
   if(comment != NULL)
      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   size = (long)rows * (long)cols;
   for(p=0;p<size;p++){      // Write the image in pixel interleaved format.
      fputc(image_all[3*p+0], fp);
      fputc(image_all[3*p+1], fp);
      fputc(image_all[3*p+2], fp);
   }

   if(fp != stdout) fclose(fp);
   return(1);
}

/*
static struct _timeb timestart;
static struct _timeb timeend;
void timer_start()
{
   _ftime( &timestart );
   COLog("timer start...\n");
}
void timer_stop()
{
   _ftime( &timeend );

   unsigned long seconds;
   unsigned long milliseconds;
   seconds = timeend.time - timestart.time;
   if (timeend.millitm > timestart.millitm)
   {
      milliseconds = timeend.millitm - timestart.millitm;
   }
   else
   {
      seconds--;
      milliseconds = (timeend.millitm + 1000) - timestart.millitm;
   }

   COLog("timer stop, elapsed %d.%d seconds.\n", seconds, milliseconds);
}

void timer_elapsed()
{
   _ftime( &timeend );

   unsigned long hours;
   unsigned long minutes;
   unsigned long seconds;
   unsigned long milliseconds;
   seconds = timeend.time - timestart.time;
   long msdif;
   msdif = timeend.millitm - timestart.millitm;
   if (msdif > 0)
      milliseconds = msdif;
   else
   {
      seconds--;
      milliseconds = (timeend.millitm + 1000) - timestart.millitm;
   }

   minutes = seconds/60;
   if (minutes == 0)
      COLogFile("elapsed %d.%d seconds.\n", seconds, milliseconds);
   else
   {
      hours = minutes/60;
      seconds = seconds - minutes*60;
      if (hours == 0)
         COLogFile("elapsed %dm%ds%dms\n", minutes, seconds, milliseconds);
      else
      {
         minutes = minutes - hours*60;
         COLogFile("elapsed %dh%dm%ds%dms\n", hours, minutes, seconds, milliseconds);
      }
   }
}
*/

time_t timestart;
time_t timeend;
void timer_start()
{
   timestart = clock();
   COLog("timer start...\n");
}
void timer_stop()
{
   timeend = clock();
   unsigned long seconds, milliseconds;
   seconds = (timeend-timestart)/CLOCKS_PER_SEC;
   milliseconds = ((100*(timeend-timestart))/CLOCKS_PER_SEC) - 100*seconds;
   COLog("timer stop, elapsed %d.%d seconds.\n", seconds, milliseconds);
}

void timer_elapsed(int prnt)
{
   timeend = clock();
   unsigned long hours, minutes, seconds, milliseconds;
   seconds = (timeend-timestart)/CLOCKS_PER_SEC;
   milliseconds = ((100*(timeend-timestart))/CLOCKS_PER_SEC) - 100*seconds;
   minutes = seconds/60;
   if (minutes == 0)
      if (prnt==0)
         COLog("elapsed %d.%d seconds. \n", seconds, milliseconds);
      else
         COLog("elapsed %d.%d seconds. \n", seconds, milliseconds);
   else
   {
      hours = minutes/60;
      seconds = seconds - minutes*60;
      if (hours == 0)
         if (prnt==0)
            COLog("elapsed %dm%ds%dms\n", minutes, seconds, milliseconds);
         else
            COLog("elapsed %dm%ds%dms\n", minutes, seconds, milliseconds);
      else
      {
         minutes = minutes - hours*60;
         if (prnt==0)
            COLog("elapsed %dh%dm%ds%dms\n", hours, minutes, seconds, milliseconds);
         else
            COLog("elapsed %dh%dm%ds%dms\n", hours, minutes, seconds, milliseconds);
      }
   }
}

/************************************/
/*       Filename Manipulation      */
/************************************/

//adds an extension (label) to a filename
void COAddExtension(char **filename, char *label)
{
	//allocate memory for new filename
	char *new_filename	= new char [strlen(*filename) + strlen(label) + 1], ext[5];


	//copy filename
	strcpy(new_filename, *filename);

	//get extension of filename (e.g. '.txt')
	char *pdest = strchr(new_filename, '.');
	strcpy(ext, pdest);

	//place filename label at the end of the filename
	//followed by extension...
	strcpy(pdest, label);
	strcat(new_filename, ext);

	//delete old filename and replace it with new one...
	delete *filename;
	(*filename)	= new_filename;

	//done.
	return;
}