#include <iostream>
#include <math.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>

using namespace std;
using namespace NTL;
int main(){
    // int n = 2048;
    // int inH = 32, inW =32;
    // int kerH=5, kerW=5,  pcol=32;  
    int n = 12;
    int inH = 4, inW = 4;
    int kerH= 1, kerW=2,  pcol = 4; 

    int  strid_H = 1,strid_W = 1, pool_stridH=2,pool_stridW=2;

    int outH = ceil((inH*1.0)/strid_H);   
    int outW = ceil((inW * 1.0)/strid_W);      
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n/(real_nrow);
    int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
    int packed_num = max_pool_num * pool_stridH * pool_stridH; 
    int packed_col = ceil(double(1.0*datcols/packed_num));
    
    int part = ceil(double(datcols) / double(pcol)); 
    int prow = part *  real_nrow; 

     int packedNum_row = pool_stridH * real_nrow; 
    int packedNum_col = pool_stridW;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); 
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));   

    cout<<"dat="<<endl;
    int dat[datrows][datcols];
    for(int z = 0; z < datrows; z++){
        for(int i = 0; i < datcols; i++){
            dat[z][i] = rand() % 100;
            cout<<dat[z][i]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    int **tempMat = new int *[real_nrow];
    for (int i = 0; i < real_nrow; i++){
        tempMat[i] = new int[datcols];
    }  
    for (int i = 0; i < real_nrow; i++){
        for (int j = 0; j < datcols; j++){
            tempMat[i][j] = 0;
        }
    }
    for (int j = 0; j < datcols; j++){
        tempMat[0][j] = 1;
    }  
    for (int i = 1; i < real_nrow; i++){
        for (int j = 0; j < datcols; j++){
            tempMat[i][j] =  dat[i - 1][j];            
        }
    }

    // data每列倒序排列
    int value = 0;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++)
    {
        for (int j = 0; j < datcols; j++)
        {
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }
    
    //每col列一截断
    int tempid = 0;
    int Pdat1[prow][pcol];
    for(int i = 0; i < part; i++){
        int rowId = i*real_nrow;
        for(int j = 0; j < pcol; j++){
            for(int k = 0; k < real_nrow; k++){
                Pdat1[rowId+k][j] = tempMat[k][tempid];
            }
            tempid++;           
        }
    }
    cout<<"Pdat1="<<endl;
    for(int i = 0; i < prow; i++){
        for(int j = 0; j < pcol; j++){
            cout << Pdat1[i][j]<<" ";
        }
    cout<<endl;
    }
    cout<<endl;



    //是否需要填0判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int p_dat2[pp_row][pp_col];
    cout<<"pp_row="<<pp_row<<endl;
    cout<<"packedNum_row="<<packedNum_row<<endl;
    cout<<"packedRowCount="<<packedRowCount<<endl;

      for(int i = 0; i < pp_row; i++){
        for(int j = 0; j < pp_col; j++){
            p_dat2[i][j] = 0;
        }
    }

    for(int i = 0; i < prow; i++){
        for(int j = 0; j < pcol; j++){
            p_dat2[i][j] = Pdat1[i][j];
        }
    }

     cout<<"Pdat2="<<endl;
    for(int i = 0; i < pp_row; i++){
        for(int j = 0; j < pp_col; j++){
            cout << p_dat2[i][j]<<" ";
        }
    cout<<endl;
    }
    cout<<endl;


    int p_dat3[packedRowCount*packedColCount][real_nrow*pool_stridW*pool_stridH]={};
    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
        
            for (int j = startColId; j <= endColId; j++){
                for (int i = startRowId; i <startRowId+real_nrow; i++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
          
            for(int j = startColId; j <= endColId; j++){
                for(int i = startRowId + real_nrow; i <= endRowId; i++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
        }
    } 

    
    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;


     cout<<"Pdat3 ="<<endl;
    for(int i = 0; i < packedRowCount*packedColCount; i++){
        for(int j = 0; j < real_nrow*pool_stridW*pool_stridH; j++){
            cout << p_dat3[i][j]<<" ";
        }
    cout<<endl;
    }
    cout<<endl;


    int p_dat4[packed_col*max_pool_num][pool_stridW*pool_stridH*real_nrow]={};
    for(int i = 0; i < packedRowCount*packedColCount; i++){
        for(int j = 0; j < pool_stridW*pool_stridH*real_nrow; j++){
            p_dat4[i][j] = p_dat3[i][j];        
        }

    }

    ZZX poly[packed_col];
    cout<<"poly="<<endl;
    for(int z = 0; z < packed_col; z++){
        poly[z].SetLength(n);
        int index = 0;
        for(int i = 0; i < max_pool_num; i++){      
            for(int j = 0; j < pool_stridW*pool_stridH*real_nrow; j++){
                SetCoeff(poly[z], index, p_dat4[z*max_pool_num+i][j]);
                index++;
            }
            cout<<poly[z]<<endl;
        }
    }
    
    return 0;
}