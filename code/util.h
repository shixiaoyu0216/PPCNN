#ifndef UTIL_H
#define UTIL_H
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>
#include <Binomialdistribution.h>
#include <ctime>
#include <enc_param.h>
#include <param.h>

using namespace std;
using namespace NTL;

class util
{
public:
    util();
    int ReverseInt(int i);
    int *rand1DGen(int lenth, int r);
    int ** rand2DGen(int maxcha,int maxrow, int r);
    int *** rand3DGen(int maxcha,int maxrow, int maxcol, int r);
    int ***abstrc3D(int cha,int row,int col, int ***p, int ***g);
    int * abstrc1D(int col, int *p, int *g);
    void rand(ZZX& r,int n,ZZ q);
    void Module(ZZX& modx,int n,ZZ q);
    void Mult(ZZX& pr,ZZ_pX pm,ZZX& p1,ZZX& p2,int n,ZZ q,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4);
    void Mult_para(ZZX& pr,ZZ_pX pm,ZZX& p1,ZZX& p2,int n,ZZ q,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4);
    void compress(ZZX& pc,int n,ZZ q,int d);
    void decompress(ZZX& pdc,int n,ZZ q,int d);
    int doubleToInt(double f);
    int floatToInt(float f);
    ZZX reangecoeff(ZZX m2,int n);
    int gcdExtended(int a, int b,int *x, int *y);
    enc_param randnum_gen();
    enc_param ***EncRandGen(int kerN, int packedRowCount,int packedColCount);



    inline int *** newArray3(int cha,int row, int col){
        int *** p3;
        p3 = new int**[cha];
        for (int i = 0; i < cha; i++) {
            p3[i] = new int*[row];
            for (int j = 0; j < row; j++)
                p3[i][j] = new int[col];
            }
        return p3;
    }

    inline ZZ *** newZZArray3(int cha,int row, int col){
    ZZ *** p3;
    p3 = new ZZ**[cha];
    for (int i = 0; i < cha; i++) {
        p3[i] = new ZZ*[row];
        for (int j = 0; j < row; j++)
            p3[i][j] = new ZZ[col];
         }
    return p3;
    }

    inline enc_param*** new3(int cha, int row, int col){
        enc_param***  usp = new enc_param**[cha];
        for (int i = 0; i < cha; i++) {
            usp[i] = new enc_param*[row];
            for (int j = 0; j < row; j++)
                usp[i][j] = new enc_param[col];
            }
        return usp;
    }





   inline void initArray3(int ***p3, int cha, int row, int col,int value) {
    for (int i = 0; i < cha; i++)
        for (int j = 0; j < row; j++)
            for (int k = 0; k < col; k++)
                p3[i][j][k] = value;
    }

    inline void printArray3(int ***p3, int cha,int row, int col) {
    cout << "[";
    for (int i = 0; i < cha; i++)
    {
        cout << "[";
        for (int j = 0; j < row; j++)
        {
            for (int k = 0; k < col; k++)
                {cout << p3[i][j][k] << ' ';}
            cout<<endl;

        }
        cout << "]"<<endl;
    }
    cout << "]" << endl;
}
 
   

    inline void deleteArray3(int ***p3, int cha,int row, int col) {
    for (int i = 0; i < cha; i++) {
        for (int j = 0; j < row; j++) {
            delete p3[i][j];
        }
        delete  p3[i];
    }
    delete p3;
    }

   inline void deleteArray2(int **p2,int row, int col) {
    for (int i = 0; i < row; i++) {
            delete[] p2[i];
        }
        delete  p2;
    }


inline int ** newArray2(int row, int col)
{
     //动态开辟空间  
    int **p2 = new int*[row]; //开辟行  
    for (int i = 0; i < row; i++)
        p2[i] = new int[col]; //开辟列  
    return p2;
}

inline void initArray2(int **p2, int row, int col,int value) {
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            p2[i][j] = value;
}


inline void printArray2(int **p2, int row, int col) {
    cout << "[";
    for (int i = 0; i < row; i++)
    {
        cout << "[";
        for (int j = 0; j < col; j++)
        {        
            cout << p2[i][j] << ' ';
        }
        cout << "]"<<endl;
    }
    cout << "]" << endl;
}




inline ZZX ** newZZXArray2(int row, int col)
{
     //动态开辟空间  
    ZZX **p2 = new ZZX*[row]; //开辟行  
    for (int i = 0; i < row; i++)
        p2[i] = new ZZX[col]; //开辟列  
    return p2;
}
 inline void deleteZZXArray2(ZZX **p2,int row, int col) {
    for (int i = 0; i < row; i++) {
            delete[] p2[i];
        }
        delete  p2;
    }

    inline ZZX *** newZZXArray3(int cha,int row, int col){
    ZZX *** p3;
    p3 = new ZZX**[cha];
    for (int i = 0; i < cha; i++) {
        p3[i] = new ZZX*[row];
        for (int j = 0; j < row; j++)
            p3[i][j] = new ZZX[col];
         }
    return p3;
    }

    inline void deleteZZXArray3(ZZX ***p3, int cha,int row, int col) {
    for (int i = 0; i < cha; i++) {
        for (int j = 0; j < row; j++) {
            delete[] p3[i][j];
        }
        delete[]  p3[i];
    }
    delete p3;
    }


   inline int **** newArray4(int num, int cha,int row, int col){
    int ****p4;
    p4=new int***[num];
    for(int z=0;z<num;z++){
        p4[z] = new int**[cha];
        for (int i = 0; i < cha; i++) {
            p4[z][i] = new int*[row];
            for (int j = 0; j < row; j++){
                p4[z][i][j] = new int[col];}
            }
          }
    return p4;
   
    }

    inline ZZX **** newZZXArray4(int num, int cha, int row, int col){
    ZZX ****p4;
    p4 = new ZZX***[num];
    for(int z=0;z<num;z++){
        p4[z] = new ZZX**[cha];
        for (int i = 0; i < cha; i++) {
            p4[z][i] = new ZZX*[row];
            for (int j = 0; j < row; j++){
                p4[z][i][j] = new ZZX[col];
                }
            }
        }

        // for(int i=0;i<num;i++){
        //     for(int j=0;j<cha;j++){
        //         for(int z=0;z<row;z++){
        //             for(int u=0;u<col;u++){
        //                 p4[i][j][z][u].SetLength(256);
        //                 // p4[i][j][z][u]=ZZX();
        //             }
        //         }
        //     }
        // }
    return p4;
}

  

    inline void deleteZZXArray4(ZZX ****p4, int num, int cha,int row, int col) {
        for(int z=0;z<num;z++){
            for (int i = 0; i < cha; i++) {
                for (int j = 0; j < row; j++) {
                    delete[] p4[z][i][j];
                }
                delete[]  p4[z][i];
            }
            delete[] p4[z];
        }
            delete p4;
    }

    inline void deleteArray4(int ****p4, int num, int cha,int row, int col) {
    for(int z=0;z<num;z++){
        for (int i = 0; i < cha; i++) {
            for (int j = 0; j < row; j++) {
                delete[] p4[z][i][j];
            }
            delete[]  p4[z][i];
        }
        delete[] p4[z];
        }
        delete p4;
    }


};
#endif // UTIL_H
