
#ifndef _CO_DEFAULTS_H
#define _CO_DEFAULTS_H

#define PI 3.1415926535
#define GTRESH PI/6
#define GTRESH_2 PI/4
#define RTRESH 7.0
#define RTRESH_2 20
#define DF_SIGMA 1.0
#define DF_HIGH 0.9
#define DF_LOW 0.5
#define DF_MINN 5

#define FC_ELLIPSE 0
#define FC_VERT_LINE 1
#define FC_HORIZ_LINE 2
#define FC_LINE 3
#define FC_SQUARE_BOX 4
#define FC_CUSTOM 5

// Line detection
#define NE_TRESH 5
#define MIN_DIST 30
#define DMAX 9
#define DMAX2 DMAX*DMAX
#define DMIN2 (DMAX-2)*(DMAX-2)
#define ALPHA_MAX 160.0*PI/180.0

//#define MIN_LINE_REZID 1.2
//#define MIN_LINE_INLIER 7
//#define MIN_LINE_INLIER 15

#define M_FRM_NO_POINT -10000000

#define AUGMENT_MIN_CORN 25
#define AUGMENT_LINE_STEP 10

extern double COSign(double);
extern int COSolveCubic(double, double, double, double, double&, double&, double&);
extern inline int CORound(double);
extern inline int COMax(int, int);
extern inline int COMin(int, int);
extern inline int COMin(double*, int);
extern inline int COIMin(double ,double , double );
extern inline double COMin3(double , double , double );
extern inline double COMax(double, double);
extern inline double COMin(double, double);
extern inline int CORoundSign(double);
extern inline double COEuclidianDistance2D(double*, double*);
extern inline double COEuclidianDistance2D2(double*, double*);
extern double COMedian(double*, int, double);
extern inline double COMedianToSigmaGaussian(double);
extern void COSort(double*, int);
extern void COISort(double*, int, int*);
extern double COFactorial(double);

extern void COLog(const char*, ...);
//extern void COLogFile(const char*, ...);
extern int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
    int cols, char *comment, int maxval);
extern int write_ppm_image(char *outfilename, unsigned char *image_red,
    unsigned char *image_grn, unsigned char *image_blu, int rows,
    int cols, char *comment, int maxval);
extern int write_ppm_image2(char *outfilename, unsigned char *image_all,
    int rows, int cols, char *comment, int maxval);
extern int read_pgm_image(char *infilename, unsigned char *image, int *rows,
    int *cols);
extern void timer_start();
extern void timer_stop();
extern void timer_elapsed(int prnt=0);

// threading
extern int gThreadProgress(float);


extern void long_job_update();
extern void long_job_update(unsigned long);
extern void long_job_end();
extern void long_job_start(unsigned long);

//image sampling functions
extern void COZoomIn(unsigned char**, unsigned char*, int, int, int, bool);
extern void COZoomOut(unsigned char**, unsigned char*, int, int, int, bool);

//file extension function
extern void COAddExtension(char**, char*);

extern "C" int lmdif(int m, int n,double x[], double fvec[] , double ftol, double xtol, double gtol, int maxfev,
						 double epsfcn,double diag[],int mode, double factor, int nprint, int *info, int *nfev,
						 double fjac[],int ldfjac, int ipvt[], double qtf[], double wa1[], double wa2[] ,
						 double wa3[], double wa4[]);

#endif