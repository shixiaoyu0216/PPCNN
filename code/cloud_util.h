#ifndef CLOUD_UTIL_H
#define CLOUD_UTIL_H
#include <NTL/ZZX.h>
#include <string>
#include <sstream>
#include <enc.h>
#include <fstream>
#include <param.h>
#include <util.h>
#include <helib/helib.h>



using namespace std;
using namespace NTL;

class Cloud_util{
  public:
  Cloud_util();
  int** plain_weight(int kerN,int kerD,int kerH, int kerW, string orignalweight,string orignalbias);
  int** plain_weightNobias(int kerN,int kerD,int kerH, int kerW, string orignalweight);
  int** dense_weight(int den_num1,int fc1_w,string w_dense1,string b_dense1);
  void matTo2Dpoly(ZZX **M,int ***ww_conv,int b_conv[],int kerN,int kerD,int kerH,int kerW);
  void matTopoly(ZZX **W_DEN1,int **wb_den1,int a,int dennum1,int fc1_w);
  void readweights(int **w_dense,int fc1_w,int den_num1,string filename);
  void readweights2D_4D(int ****w_conv_int,int kerH,int kerW,int kerD,int kerN,string filename);
  void readbias(int b_conv[], int len,string filename);
  int doubleToInt(double f);
  void matPacking_by_poolsize(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
  void slideToMat_SAME(int ***C, int ***A, int inD,int inH, int inW, int kerH,int kerW,int strid_H,int strid_W);
  int** mat_mult(int **W, int **X, int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W);
 void mat_mult_twoPart(int ***WX2D1, int ***WX2D2, int **W, int **X, int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W);
  void matPacking_by_poolsize12(ZZX** polys, int **dat, int datrows, int datcols, int pcol);
  ZZX * GentNoiseforRelu(int kerN, int inH, int inW,  int range, helib::EncryptedArray ea);
  // ZZX GentNoiseforReluFC(int len,int range, helib::EncryptedArray ea);
  int*mat_fc_mult(int **W, int *X, int X_num, int W_col);
  void GentNoise(ZZX &nr_poly, ZZX &nr_invpoly, int slots, int p,long randnum,helib::EncryptedArray ea);
  void GentNoiseforRelu_Conv(ZZX* nr_poly, ZZX *nr_invpoly,int kerN, int inH, int inW, long p,int range,helib::EncryptedArray ea);
  void GentNoiseforRelu_FC(int len, long p, int range, ZZX* nr_poly, ZZX *nr_invpoly,helib::EncryptedArray ea);
  void Encode_Enc_one(int AA, ZZX *u,ZZX v,helib:: EncryptedArray ea);
 void packingbyPoolWin(int kerH, int kerW, ZZX *poly, int *dat, int datcols, int pcol, int pool_stridH, int pool_stridW, int max_pool_num);


};
#endif // CLOUD_UTIL_H