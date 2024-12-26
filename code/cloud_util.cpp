#include <cloud_util.h>

using namespace std;
using namespace NTL;





Cloud_util::Cloud_util(){

}

int** Cloud_util::plain_weight(int kerN,int kerD,int kerH, int kerW, string orignalweight,string orignalbias)
{      
   util U;
    Cloud_util cdu;
    int ****w_conv2= U.newArray4(kerN,kerD,kerH,kerW);
    cdu.readweights2D_4D(w_conv2,kerH,kerW,kerD,kerN,orignalweight);
    // cout<<"w_conv2="<<endl;
    // for(int i=0;i<kerN;i++){
    //     for(int j=0;j<kerD;j++){
    //         for(int z=0;z<kerH;z++){
    //             for(int v=0;v<kerW;v++)
    //                 cout<<w_conv2[i][j][z][v]<<" ";
    //             cout<<endl;
    //             }
    //         cout<<endl;
    //         }
    // cout<<endl;
    // } 
 
    // cout<<kerD<<endl;
    // cout<<w_conv2[0][2][3][4]<<endl;
    // cout<<w_conv2[2][1][3][3]<<endl;
    // cout<<w_conv2[3][1][0][2]<<endl;    
    //四维变三维
    int ***w3Dconv2=U.newArray3(kerN,kerD,kerH*kerW);         
    for(int u=0;u<kerN;u++) 
        for(int v=0;v<kerD;v++)
            for(int i=0;i<kerH;i++)
                for(int j=0;j<kerW;j++)
                    w3Dconv2[u][v][i*kerW+j]=w_conv2[u][v][i][j];
          

    // cout<<"4D变3D"<<endl;
    // cout<<"w3Dconv2="<<endl;
    // for(int i=0;i<kerN;i++){
    //     for(int j=0;j<kerD;j++){
    //         for(int z=0;z<kerH*kerW;z++)
    //             cout<<w3Dconv2[i][j][z]<<" ";
    //         cout<<endl;
    //     }
    //         cout<<endl;
    // }
     
   int **w2Dconv2=U.newArray2(kerN,kerD*kerH*kerW);       
    for(int u=0;u<kerN;u++) {
        for(int v=0;v<kerD;v++){
            for(int j=0;j<kerH*kerW;j++){
                w2Dconv2[u][v*kerH*kerW+j]=w3Dconv2[u][v][j];
            }
        }
    }
    // cout<<"3D变2D"<<endl;

    // for(int i=0;i<kerN;i++){
    //     for(int j=0;j<kerD*kerH*kerW;j++){
    //         cout<<w2Dconv2[i][j]<<" ";
    //     }
    //     cout<<"%%%%%%%%%%"<<endl;
    //     }
    //读取b
    int b_conv2[kerN];
    cdu.readbias(b_conv2,kerN,orignalbias);
    // cout<<"b="<<endl;
    // for(int i=0;i<kerN;i++)
    //     cout<<b_conv2[i]<<" ";
    // cout<<endl;

    Cloud_util cud;
      // cout<<3456*double(1.0/100)<<endl;
     // cout<<cud.doubleToInt(0.0193*w_precision)<<endl;
      
    int **wb2Dconv2=U.newArray2(kerN,kerD*kerH*kerW+1);       
    for(int u=0;u<kerN;u++)
    {
        wb2Dconv2[u][0]=b_conv2[u];
    }
    for(int u=0;u<kerN;u++) {
        for(int v=1;v<kerD*kerH*kerW+1;v++){
            wb2Dconv2[u][v]=w2Dconv2[u][v-1];
            }
    }
    // for(int u=0;u<kerN;u++) {
    //     for(int v=0;v<kerD*kerH*kerW+1;v++){
    //          cout<<wb2Dconv2[u][v]<<" ";
    //         }
    //         cout<<endl;
    // }
    // cout<<endl;

    U.deleteArray4(w_conv2,kerN,kerD,kerH,kerW);
    U.deleteArray3(w3Dconv2,kerN,kerD,kerH*kerW);
    U.deleteArray2(w2Dconv2,kerN,kerD*kerH*kerW);
    return wb2Dconv2;
}

int** Cloud_util::plain_weightNobias(int kerN,int kerD,int kerH, int kerW, string orignalweight)
{      
    util U;    
    int **w_int = U.newArray2(kerN, kerD * kerH * kerW);     
    double data[kerN][kerD * kerH * kerW];
    ifstream file(orignalweight, ios::in);
    if(!file){
        cout<<"Cannot open weights!\n";
        exit(0);
    }
    if (file.is_open()){
        for(int i = 0; i < kerN; i++){
            for (int j = 0; j < kerD * kerH * kerW; j++){
                file >> data[i][j];
            }
        }
    }  
    file.close();  
  
    for(int i = 0; i < kerN; i++){
        for (int j = 0; j< kerD * kerH * kerW; j++){                                             
            w_int[i][j] = U.doubleToInt(w_precision * data[i][j]);
        }
    }
    return w_int;
}

int**  Cloud_util::dense_weight(int den_num1,int fc1_w,string w_dense1,string b_dense1)
{   
    // util U;
     //int den_num1=120;
    //int fc1_w=784;
    enc_param usp;
    Cloud_util cdu;

    int **w_den1=new int *[den_num1];
    for (int i=0;i<den_num1;i++){
        w_den1[i]=new int[fc1_w];
    }
    // cout<<"den_num1  "<<den_num1<<endl;
    // cout<<"fc1_num1  "<<fc1_w<<endl;
    cdu.readweights(w_den1,den_num1,fc1_w,w_dense1);
    // double **w_den1=new double *[fc1_w];
    // for (int i=0;i<fc1_w;i++){
    //     w_den1[i]=new double[den_num1];
    // }
    // cdu.readweights()
    // cdu.readweights(w_den1,fc1_w,den_num1,w_dense1);

    // cout<<"reading weights in FC1 ============="<<endl;
    // cout<<"2-5"<<w_den1[2][5]<<endl;
    // cout<<"5-2"<<w_den1[5][2]<<endl;
    // cout<<"W3="<<endl;
    // for(int i=0;i<fc1_w;i++){
    //     for(int j=0;j<den_num1;j++){
    //     cout<<w_den1[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
  
    int b_fc1[den_num1];
    cdu.readbias(b_fc1,den_num1,b_dense1); 
    // cout<<"b_fc1="<<endl;
    // for(int j=0;j<den_num1;j++)
    //     cout<<b_fc1[j]<<" ";   
    // cout<<endl;
  

    int **wb_fc1=new int *[den_num1];
    for (int i=0;i<den_num1;i++){
        wb_fc1[i]=new int[fc1_w+1];
    }   

    for(int i=0; i< den_num1; i++){
        wb_fc1[i][0] = b_fc1[i];
        for(int j=1; j<fc1_w+1; j++){
            wb_fc1[i][j] = w_den1[i][j-1];
        }
    }

    // for(int j=0;j<den_num1;j++){
    //     wb_fc1[0][j]=b_fc1[j];
    //     }
    // for(int i=1;i<fc1_w+1;i++){
    //     for(int j=0;j<den_num1;j++){
    //         wb_fc1[i][j]=w_den1[i-1][j];           
    //     }
    // }
    // cout<<"Wb="<<endl;
    // double **wb_fc1=new double *[fc1_w+1];
    // for (int i=0;i<fc1_w+1;i++){
    //     wb_fc1[i]=new double[den_num1];
    // }   
              
    // for(int j=0;j<den_num1;j++){
    //     wb_fc1[0][j]=b_fc1[j];
    //     }
    // for(int i=1;i<fc1_w+1;i++){
    //     for(int j=0;j<den_num1;j++){
    //         wb_fc1[i][j]=w_den1[i-1][j];           
    //     }
    // }
    // cout<<"Wb="<<endl;
    // for(int i=0;i<fc1_w+1;i++){
    //  for(int j=0;j<den_num1;j++){
    //     cout<<wb_fc1[i][j]<<" ";
    //     }
    // cout<<endl;
    // }
    return wb_fc1;
}

 int** Cloud_util::mat_mult(int **W, int **X, int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W)
 {  util U;
    const int outH = ceil(inH/strid_H);
    const int outW = ceil(inW/strid_W);
    int**WX = U.newArray2(kerN, outH*outW);
    for(int i = 0; i < kerN; i++){
        for(int z = 0; z < outH*outW;z++){  
            WX[i][z]=0;       
            for(int j = 0; j < kerD*kerH*kerW; j++){                                 
                    WX[i][z] += W[i][j]*X[j][z];            
            }
        }
    }

    // int **XT = U.newArray2(outH*outW, kerD*kerH*kerW);
    // for(int i = 0; i <kerD*kerH*kerW; i++){
    //     for(int j = 0; j < outH*outW; j++){
    //         XT[j][i] = X[i][j];
    //     }
    // }
    //  for(int i = 0; i < kerN; i++){
    //     for(int z = 0; z < outH*outW;z++){  
    //         WX[i][z]=0;       
    //         for(int j = 0; j < kerD*kerH*kerW; j++){                                 
    //                 WX[i][z] += W[i][j]*XT[z][j];            
    //         }
    //     }
    // }

    return WX;
 }

 
 void Cloud_util::mat_mult_twoPart(int ***WX2D1, int ***WX2D2, int **W, int **X, int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W)
 {
    util U;
    const int outH=ceil(inH/strid_H);
    const int outW=ceil(inW/strid_W);
    int**WX_1=U.newArray2(kerN,outH*outW);
    int**WX_2=U.newArray2(kerN,outH*outW);
    cout<<"i am in 236"<<endl;     
    for(int i=0;i<kerN;i++){
        for(int j=0;j<outH*outW;j++){
            WX_1[i][j]=0;
            WX_2[i][j]=0;
        }
    }
    cout<<"i am in 243"<<endl; 
    for(int i=0;i<kerN;i++){
        for(int z=0;z<outH*outW;z++){
            for(int j=0;j<kerH*kerW+1;j++){                       
                WX_1[i][z]+=W[i][j]*X[j][z];          
            }
            for(int j=kerH*kerW+1;j<kerD*kerH*kerW;j++){                        
                WX_2[i][z]+=W[i][j]*X[j][z];          
            }
        }
    }
   cout<<"i am in 254"<<endl; 
 
  int**WX_precise1=U.newArray2(kerN,outH*outW);
  int**WX_precise2=U.newArray2(kerN,outH*outW);
    Cloud_util cdu;
    for(int i=0;i<kerN;i++){
        for(int j=0;j<outH*outW;j++){
            WX_precise1[i][j]=cdu.doubleToInt(WX_1[i][j]*double(1.0/w_precision));
            WX_precise2[i][j]=cdu.doubleToInt(WX_2[i][j]*double(1.0/w_precision));
        }
    }
 cout<<"i am in 265"<<endl; 
    // cout<<"二维变三维"<<endl;
    // int ***WX2D1=U.newArray3(kerN,outH,outW);      
    // int ***WX2D2=U.newArray3(kerN,outH,outW);    
    for(int z=0;z<kerN;z++){ 
        for(int i=0;i<outH;i++){
            for(int j=0;j<outW;j++){
                    WX2D1[z][i][j]=WX_precise1[z][i*outW+j]; 
                    WX2D2[z][i][j]=WX_precise2[z][i*outW+j];
                    cout<<"i am in 274"<<endl;
            }

        }
    }
     cout<<"i am in 277"<<endl;                  
    U.deleteArray2(WX_1,kerN,outH*outW);   
    U.deleteArray2(WX_2,kerN,outH*outW); 
    // return WX2D1;
 }





 int Cloud_util::doubleToInt(double f){
    int i=0;
    if(f>0)
    i=(f*10+5)/10;
    else if(f<0)
    i=(f*10-5)/10;
    else i=0;
    return i;
}


//W按照方法2装配为多项式
void Cloud_util::matTo2Dpoly(ZZX **W,int ***ww_conv,int b_conv[],int kerN,int kerD,int kerH,int kerW)
{
  
    for(int i = 0; i < kerN; i++){
        W[i][0].SetLength(n);

        SetCoeff(W[i][0], 0, b_conv[i]);

        for(int z = 0; z < kerH * kerW; z++){
            SetCoeff(W[i][0], z + 1, ww_conv[i][0][z]);
        }
    }

     //其余列不用装配偏置,但是第一个位置装配为0
    for (int i = 0; i < kerN; i++){

        for(int j = 1; j < kerD; j++){

            W[i][j].SetLength(n);

            for(int z = 0; z < kerH * kerW; z++){  

                SetCoeff(W[i][j], 0, 0);

                SetCoeff(W[i][j], z + 1, ww_conv[i][j][z]);
            }
        }
    }
}

void Cloud_util::matTopoly(ZZX **W_DEN1,int **wb_den1,int a,int dennum1,int fc1_w)
{   
    util U;
    int **CC=new int*[dennum1];
    for(int i=0;i<dennum1;i++)
        CC[i]=new int[a*n];
    //将权值装配入多项式，不足位补0.
    
    for(int i=0;i<dennum1;i++){
        for(int j=0;j<a*n;j++){
            CC[i][j]=0;      
            }
    }   
    for(int i=0;i<dennum1;i++){
        for(int j=0;j<fc1_w+1;j++){
            CC[i][j]=wb_den1[i][j];
        }
    }
    for(int i=0;i<dennum1;i++){
        for(int j=0;j<a;j++){
            W_DEN1[i][j].SetLength(n);
            for(int z=0;z<n;z++){
                SetCoeff(W_DEN1[i][j],z,CC[i][j*n+z]);
            }
          //  cout<<W_DEN1[i][j]<<endl;
        }        
    }
    U.deleteArray2(CC,dennum1,a*n);
}

void Cloud_util::readweights(int **w_dense,int den_num1,int fc1_w,string filename){    
  
    util U;
    double **w=new double*[den_num1];
    for (int i=0;i<den_num1;i++){
        w[i]=new double[fc1_w];
    } 
   
   
    ifstream file(filename, ios::in);
    if(!file){
        cout<<"Cannot open originalweights.csv in fc!\n";
        exit(0);
    }
   
    if (file.is_open()){
        for(int i=0;i<den_num1;i++){
            for (int j=0;j<fc1_w;j++){           
                file>>w[i][j];             
                w_dense[i][j]=U.doubleToInt(w_precision * w[i][j]);
            }
        }
    }
    file.close();

    for (int i = 0; i < den_num1; i++) {
            delete[] w[i];
        }
        delete  w;
    
}

void Cloud_util::readweights2D_4D(int ****w_conv_int,int kerH,int kerW,int kerD,int kerN,string filename)
{    util U;
    int size2D=kerD*kerH*kerW;
    double data[kerN][size2D];
    double w_conv[kerN][kerD][kerH][kerW];
    ifstream file(filename, ios::in);
    if(!file){
        cout<<"Cannot open weights!\n";
        exit(0);
    }

    if (file.is_open()){
        for(int i=0;i<kerN;i++){
            for (int j=0;j<size2D;j++){
                file>>data[i][j];
            }
        }
    }  
    file.close();
  
    int count=0;
    for(int i=0;i<kerN;i++){
        count=0;
        for (int j=0;j<kerD;j++){  
            for (int k1=0;k1<kerH;k1++){
                for (int k2=0;k2<kerW;k2++) {                
                    w_conv[i][j][k1][k2] = data[i][count];                        
                    w_conv_int[i][j][k1][k2] = U.doubleToInt(w_precision*w_conv[i][j][k1][k2]);
                    count++;
                }
        //cout<<endl;
            }
        //cout<<endl;
        }
    }
}

void Cloud_util::slideToMat_SAME(int ***C, int ***A, int inD,int inH, int inW, int kerH,int kerW,int strid_H,int strid_W)
{
    util U;
    const int outH=ceil((inH*1.0)/strid_H);
    const int outW=ceil((inW*1.0)/strid_W);
    const int pad_W=max((outW-1)*strid_W+kerW-inW,0);
    const int pad_H=max((outH-1)*strid_H+kerH-inH,0);
    const int pad_top=int((pad_H)/2);
    const int pad_left=int((pad_W)/2);
 
   
    int ***B=U.newArray3(inD,inH+pad_H,inW+pad_W);
    
    U.initArray3(B,inD,inH+pad_H,inW+pad_W,0);
    //pading B in A
   
    for(int z=0;z<inD;z++){
        int i=0;
        for (int u=0;u<inH;u++){ 
            int j=0;
            for(int v=0;v<inW;v++){ 
                B[z][i+pad_top][j+pad_left]=A[z][u][v];
                j=j+1;
            }
            i=i+1;
        }
    }  

    // //cout A and B
    // cout<<"A=";
    // for (int i=0;i<inH;i++){
    //     for(int j=0;j<inW;j++){
    //         cout<<A[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"B=";
    // for (int i=0;i<inH+pad_H;i++){
    //     for(int j=0;j<inW+pad_W;j++){
    //         cout<<B[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    
     //after padding,slide window on B and transform to C
    for(int z=0;z<inD;z++){
        int y=0;
        for(int i=0;i<outH;i++){
            int start_u = strid_H * i;
            for(int j=0;j<outW;j++){
                int x=0;
                int start_v = strid_W * j;
                for(int u=start_u;u<start_u+kerH;u++){
                    for(int v=start_v;v<start_v+kerW;v++){
                        C[z][x][y]=B[z][u][v];
                        x=x+1;
                    }
                }
                y=y+1;
                }
            }    
        }
    U.deleteArray3(B,inD,inH+pad_H,inW+pad_W);
}
   
void Cloud_util::matPacking_by_poolsize(ZZX** polys, int **dat, int datrows, int datcols, int pcol){
    
    util U;
    //dat:datrow*datcols; 加上bias 行数多一变为 tempMat:real_nrow*（pcol*part）i.e.26*784; Pmat:prow*pcol

    //  int datrows=25,datcols=784;//输入矩阵dat 的行和列
    
    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows + 1;
    int prow = part * (datrows + 1);                 //将dat28列一截断生成的矩阵Pdat的行和列
    
    int** Pdat = U.newArray2(prow,pcol);
    U.initArray2(Pdat,prow,pcol,0);

    
    //tempMat是第一行全为1，其余值dat矩阵的值大小为26*784
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
    // cout <<value<<endl;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++){
        for (int j = 0; j < datcols; j++){
            value = tempMat[i][j];
            tempMat[i][j] = tempMat[real_nrow - 1 - i][j];
            tempMat[real_nrow - 1 - i][j] = value;
        }
    }


    int eleCeil = datrows + 1;       // 1:bias
    int packedNum_row = 2 * eleCeil; //26*2
    int packedNum_col = 4;
    // int unitNum = 2 * 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); //7
    int packedRowCount = ceil(prow / (packedNum_row * 1.0)); //14
 
    int tempid = 0;
    for(int i=0; i<part; i++){
        int rowId = i*real_nrow;
        for(int j=0; j<pcol; j++){
            for(int k=0; k<real_nrow; k++){
                Pdat[rowId+k][j] = tempMat[k][tempid];
            }
            tempid++;           
        }
    }
  
    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for (int j = startColId; j <= endColId; j++){
                for (int i = startRowId; i <startRowId+eleCeil; i++){
                    SetCoeff(polys[row][col], index1, Pdat[i][j]);
                    index1++;
                }
            }
            for(int j=startColId; j<=endColId; j++){
                for(int i=startRowId+eleCeil; i<=endRowId; i++){
                    SetCoeff(polys[row][col], index1, Pdat[i][j]);
                    index1++;
                }
            }
        }
    }

    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }
    //  cout<<polys[3][6]<<endl;

    for (int i = 0; i < real_nrow; i++)
    {
        delete[] tempMat[i];
    }
    delete[] tempMat;

    U.deleteArray2(Pdat,26*28,28);

    // return polys;
}


void Cloud_util::readbias(int b_conv[], int len,string filename ){  
    util U;
    double b[len];
    ifstream file1(filename, ios::in);
    if (!file1){
        cout<<"Cannot open originalbias!\n";
        exit(0);
    }
   
    if (file1.is_open()){
        for (int z=0;z<len;z++){
            file1>>b[z];  
            b_conv[z]=U.doubleToInt(w_precision*b[z]);
        }
    }
    file1.close();
}


// cloud_enc_param Cloud_util::randNumGen(){
//      cloud_enc_param usp;
//     for(int i=0;i<k;i++){
//         usp.r1[i].SetLength(n);
//         usp.Bin_dis.Binomial_distribution_gen(usp.r1[i],n,ita);
//     }
//     for(int i=0;i<k;i++){
//        usp.e11[i].SetLength(n);
//        usp.Bin_dis.Binomial_distribution_gen(usp.e11[i],n,ita);
//     }
//     usp.Bin_dis.Binomial_distribution_gen(usp.e12,n,ita);
//     return usp;  
// }

void Cloud_util::matPacking_by_poolsize12(ZZX** polys, int **dat, int datrows, int datcols, int pcol){
    util U;
    int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为几个部分
    int real_nrow = datrows;
    int prow = part * datrows; //将dat 7 列一截断生成的矩阵Pdat的行和列
    int** Pdat = U.newArray2(prow,pcol);
    U.initArray2(Pdat,prow,pcol,0);
    // data每列倒序排列
    int value = 0;
    int iters = real_nrow / 2;
    for (int i = 0; i < iters; i++){
        for (int j = 0; j < datcols; j++){
            value = dat[i][j];
            dat[i][j] = dat[real_nrow - 1 - i][j];
            dat[real_nrow - 1 - i][j] = value;
        }
    }
    // for (int i = 0; i < real_nrow; i++){
    //     for (int j = 0; j < datcols; j++){            
    //         cout<<dat[i][j]<<" ";            
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"每part=4个部分一截断"<<endl;    
    int eleCeil = datrows;  
    int packedNum_row = 2 * eleCeil; //37*2
    int packedNum_col = 4;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); //4
    int packedRowCount = ceil(prow / (packedNum_row * 1.0)); //4
    // ZZX **polys = U.newZZXArray2(packedRowCount, packedColCount);//6*4*4

    int tempid = 0;
    for(int i=0; i<part; i++){
        int rowId = i*real_nrow;
        for(int j=0; j<pcol; j++){
            for(int k=0; k<real_nrow; k++){
                Pdat[rowId+k][j] = dat[k][tempid];
            }
            tempid++;           
        }
    }
    //  cout<<"Pdat="<<endl;
    //     for(int i=0;i<prow;i++){
    //         for(int j=0;j<pcol;j++){
    //             cout<<Pdat[i][j]<<" ";
    //         }
    //         cout<<endl;
    //     }
    //     cout<<endl;

    //是否需要填0 判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int** p_dat2 = U.newArray2(pp_row, pp_col);
    U.initArray2(p_dat2,pp_row,pp_col,0);
    for(int i = 0; i<prow; i++){
        for(int j=0; j<pcol; j++){
            p_dat2[i][j] = Pdat[i][j];
        }
    }   

    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for (int j = startColId; j <= endColId; j++){
                for (int i = startRowId; i <startRowId+eleCeil; i++){
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }
          
            for(int j=startColId; j<=endColId; j++){
                for(int i=startRowId+eleCeil; i<=endRowId; i++){
                    SetCoeff(polys[row][col], index1, p_dat2[i][j]);
                    index1++;
                }
            }
        }
    }
 
    // cout<<"packing result***************"<<endl;
    // for(int i=0; i<packedRowCount; i++){
    //     for(int j=0; j<packedColCount; j++){
    //         cout<<i<<"--"<<j<<":"<<polys[i][j]<<endl;
    //     }
    //     cout<<endl;
    // }
    //  cout<<polys[3][6]<<endl;
    U.deleteArray2(Pdat,prow,pcol);
    U.deleteArray2(p_dat2,pp_row,pp_col);
}

   //生成噪声，并计算出相应的乘法模逆
    void Cloud_util::GentNoise(ZZX &nr_poly, ZZX &nr_invpoly, int slots, int p,long randnum,helib::EncryptedArray ea){
        util U;
        std::vector<long > nr(slots,0);
        std::vector<long > nr_inv(slots,0); 
        for(int i=0; i<n; i++){
            nr[i]=RandomBnd(randnum)%1000+1;
            // nr[i]=1; //随机数为1
            int x=0,y=0;
            U.gcdExtended(nr[i],p,&x,&y);
            //nr*x=1 mod p,  p2和p互素，即p2取小于p的大于0的整数，p2的乘法模逆唯一。
            nr_inv[i]=x;          
        }
        ea.encode(nr_poly, nr);
        ea.encode(nr_invpoly, nr_inv);
    }



    void Cloud_util::GentNoiseforRelu_Conv(ZZX* nr_poly, ZZX* nr_invpoly,int kerN, int inH, int inW, long p,int range,helib::EncryptedArray ea){
       
        Cloud_util cdu;
        int slots = ea.size();
        if(inH * inW <= slots){
            for(int i = 0; i < kerN; i++){
                cdu.GentNoise(nr_poly[i], nr_invpoly[i], slots, p, range, ea);
            }
        }        
    }

    void Cloud_util::GentNoiseforRelu_FC(int len, long p, int range, ZZX* nr_poly, ZZX *nr_invpoly,helib::EncryptedArray ea){
       
        Cloud_util cdu;
        int slots=ea.size();
        int num=ceil(1.0*len/slots);
       
        for(int i=0;i<num;i++){
            cdu.GentNoise(nr_poly[i], nr_invpoly[i], slots, p,range,ea);
        }
    
        

    }

    // ZZX Cloud_util::GentNoiseforReluFC(int len,int range, helib::EncryptedArray ea){
    
    // int slots=ea.size();
    // int *d=U.rand1DGen(len,range);
    // std::vector<long> p1(slots,0); 
    // //每一个通道编码为一个多项式 
    // for(int z1=0;z1<len;z1++){                          
    //     p1[z1]= d[z1];
    // }        
    // ZZX encode1_poly;
    // ea.encode(encode1_poly,p1);
    // return encode1_poly;
    // }

void Cloud_util::Encode_Enc_one(int AA, ZZX *u,ZZX v,helib:: EncryptedArray ea){
    util U;
    long slots=ea.size();
    std::vector<long > p(slots,0);
    ZZX encode_poly;
    p[0]=AA;
    ea.encode(encode_poly,p);
    
  
    ZZX t[k];
    ifstream inFile1;
    inFile1.open("./file/key/pk-t.csv",ios::in); 
    if(!inFile1){
        cout<<"Cannot open pk-t in conv1!\n";
        exit(0);
    }   
    for(int i=0;i<k;i++){
        inFile1>>t[i];       
    }    
    inFile1.close();    
    for(int i=0;i<k;i++){
        U.decompress(t[i],n,q,dt);
    }
    //read A in pk-A
    ZZX **A=U.newZZXArray2(k,k);
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv",ios::in);
    if(!inFile2){
        cout<<"Cannon open pk-A in conv1!\n";
        exit(0);
    }
    for(int i=0;i<k;i++){
        for(int j=0;j<k;j++){
            inFile2>>A[i][j];     
        }
    }    
    inFile2.close();        
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
    //Enc
    Enc enc;
    //q/2d
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,dp);//two_d=2^dp
    two_d_rr=to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d=RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)
    
   
    enc_param usp;
    Cloud_util cdu;
    usp= U.randnum_gen();
    enc.encrypt(encode_poly,u,v,n,q,k,q2d,usp.r1, A, usp.e11,usp.e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 
    

    // ofstream outFile1;
    // outFile1.open(cipherfile,ios::ate);
    // for(int z=0;z<inD;z++){
    //     for(int j=0;j<k;j++)
    //         outFile1<<u[z][j]<<endl;
    //     outFile1<<v[z]<<endl;
    // }    
    // outFile1.close(); 
    U.deleteZZXArray2(A,k,k); 
 
}

    int*Cloud_util:: mat_fc_mult(int **W, int *X, int X_num, int W_col){
    // cout<<"W*X="<<endl;
   
    int*WX=new int[X_num];
    for(int i=0;i<X_num;i++)
        WX[i]=0;

    for(int i=0;i<X_num;i++){
        for(int j=0;j<W_col+1;j++){
            WX[i] = WX[i] + W[i][j] * X[j];                         
            // XW[j]=XW[j]+X[i]*W[i][j];            
            }
        }
    // double*XW=new double[W_col];
    // for(int j=0;j<W_col;j++)
    //     XW[j]=0;

    // for(int j=0;j<W_col;j++){
    //     for(int i=0;i<X_num+1;i++){                           
    //         XW[j]=XW[j]+X[i]*W[i][j];            
    //         }
    
    //     }
    // cout<<"The fc1 X*W="<<endl;

    // int*WX_precise=new int[X_num];
    // Cloud_util cdu;
    // for(int i=0;i<X_num;i++){        
    //         WX_precise[i]=cdu.doubleToInt(WX[i]*double(1.0/w_precision));
    // }
   
    // Client_util cu;
    // int *relu_WX=new int[X_num];
    // for(int i=0;i<X_num;i++){
    //     relu_WX[i]=cu.relu(WX_precise[i]);
        // cout<<XW[i]<<" ";
            // }
            // cout<<endl;
    // cout<<"relu="<<endl;
    // for(int i=0;i<W_col;i++)
        // cout<<relu_XW[i]<<" ";
        // cout<<endl;
  return  WX;
 }


 void Cloud_util::packingbyPoolWin(int kerH, int kerW, ZZX *poly, int *dat, int datcols, int pcol, int pool_stridH, int pool_stridW, int max_pool_num)
{
  int part = ceil(double(datcols) / double(pcol)); //part是将dat的列要截断为32个部分
    int prow = part ;   //将dat 32 列一截断生成的矩阵Pdat的行和列   

    //每col列一截断
    int Pdat1[prow][pcol];
    for(int i = 0; i < prow; i++){       
        for(int j = 0; j < pcol; j++){            
            Pdat1[i][j] = dat[i*pcol+j];
        }
    }
    //  cout<<"Pdat1="<<endl;  
    // for(int i = 0; i < prow; i++){
    //     for(int j = 0; j < pcol; j++){
    //         cout << Pdat1[i][j]<<" ";           
    //     }
    //        cout<<endl;
    // }
    //    cout<<endl<<"end Pdat1"<<endl; 

      
    int packedNum_row = pool_stridH; 
    int packedNum_col = pool_stridW;
    int packedColCount = ceil(pcol / (packedNum_col * 1.0)); 
    int packedRowCount = ceil(prow / (packedNum_row * 1.0));   

    //是否需要填0判断
    int pp_row = packedRowCount * packedNum_row;
    int pp_col = packedColCount * packedNum_col;
    int p_dat2[pp_row][pp_col];

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


    int p_dat3[packedRowCount*packedColCount][pool_stridW*pool_stridH];
    for (int row = 0; row < packedRowCount; row++){
        int startRowId = row * packedNum_row;
        int endRowId = (row + 1) * packedNum_row - 1;
        for (int col = 0; col < packedColCount; col++){
            int startColId = col * packedNum_col;
            int endColId = (col + 1) * packedNum_col - 1;
            int index1 = 0;
            for(int i = startRowId; i<=endRowId; i++){
                for(int j = startColId; j<=endColId; j++){
                    p_dat3[row*packedColCount+col][index1] = p_dat2[i][j];
                    index1++;
                }
            }
 
        }
    }

    int poly_num = ceil(1.0 * packedRowCount * packedColCount /max_pool_num);
    int p_dat4[poly_num*max_pool_num][pool_stridW*pool_stridH]={};

   
    for(int i = 0; i < packedRowCount*packedColCount; i++){
        for(int j = 0; j < pool_stridW*pool_stridH; j++){
            p_dat4[i][j] = p_dat3[i][j];        
        }

    }

    for(int z = 0; z < poly_num; z++){
        poly[z].SetLength(n);
        int index = 0;
        for(int i = 0; i < max_pool_num; i++){      
            for(int j = 0; j < pool_stridW*pool_stridH; j++){
                SetCoeff(poly[z], kerH*kerW+(kerH*kerW+1)*index, p_dat4[z*max_pool_num+i][j]);
                index++;
            }
        }
    }

}


