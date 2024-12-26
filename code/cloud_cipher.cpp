#include "cloud_cipher.h"
#include <sys/time.h>
#include <helib/helib.h>
using namespace std;
using namespace NTL;

Cloud_cipher::Cloud_cipher(){

}


void Cloud_cipher::homMultConv_Conv1(int packedRowCount, int packedColCount, int packedRowUnit, int packedColUnit, string outputcipherfile, string weightplainfile, const char *inputfile, int kerN, int kerD, int kerH, int kerW, int inH, int inW, int strid_H, int strid_W, string padding){
  
    util U;
    int outH = ceil(inH / strid_H);
    int outW = ceil(inW / strid_W);
    ZZX W[kerN]; //6个元素是多项式
    for(int i=0; i<kerN; i++){
        W[i] = ZZX();
    }

    ifstream inFile1;
    inFile1.open(weightplainfile, ios::in);
    if (!inFile1){
        cout << "Cannot open weight plain file in homMult!" << endl;
        exit(0);
    }

    for (int z = 0; z < kerN; z++){
        inFile1 >> W[z];
    }
    inFile1.close();

        //read input cipher
    ZZX ***u2 = U.newZZXArray3(packedRowCount, packedColCount, k);
    ZZX **v2 = U.newZZXArray2(packedRowCount, packedColCount);

    fstream fr1(inputfile);
    if (!fr1.is_open()){
        cout << "Cannot open input file in homMult!" << endl;
        exit(1);
    }

    //云读取客户的密文文件并解压缩
    timeval t1,t2;
    gettimeofday(&t1,NULL);
    for (int i = 0; i < packedRowCount; i++){
            for (int j = 0; j < packedColCount; j++)//7
            {
                for(int z=0; z<k; z++){
                    fr1 >> u2[i][j][z];
                }
                fr1 >> v2[i][j];
            }
        }
        fr1.close();

    // cout<<"finish reading the X"<<endl;
    for (int z = 0; z < packedRowCount; z++){
        for (int i = 0; i < packedColCount; i++){
            for (int j = 0; j < k; j++){
                    U.decompress(u2[z][i][j], n, q, du);
            }
            U.decompress(v2[z][i], n, q, dv);
        }
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "conv1 云读取客户的密文文件并解压缩:" << deltaT1/1000<<"ms"<< endl;

    
        
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX pp_temp1, pp_temp2, pp_temp3, pp_temp4;

    timeval t3, t4;
    gettimeofday(&t3, NULL);
    ZZX*** p1 = U.newZZXArray3(kerN, packedRowCount,packedColCount); 
    ZZX*** p2 = U.newZZXArray3(kerN, packedRowCount,packedColCount);
    ZZX*** p3 = U.newZZXArray3(kerN, packedRowCount,packedColCount);  
    for (int i = 0; i < kerN; i++){
        int index = 0;
        for (int j = 0; j < packedRowCount; j++){
            for (int u = 0; u < packedColCount; u++){
                U.Mult(p1[i][j][u], pm, W[i], v2[j][u], n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
                U.Mult(p2[i][j][u], pm, W[i], u2[j][u][0], n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
                U.Mult(p3[i][j][u], pm, W[i], u2[j][u][1], n, q2, pp_temp1, pp_temp2, pp_temp3, pp_temp4);
            }
        }
    }
    gettimeofday(&t4, NULL);
    double deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
    cout << "conv1:Homorphic eval time is:" << deltaT2/1000<<"ms"<< endl;

    ofstream outFile4;
    outFile4.open(outputcipherfile, ios::ate);
    for (int z = 0; z < kerN; z++){
        for(int i=0; i<packedRowCount; i++){
            for (int j = 0; j < packedColCount; j++){
                outFile4 << p1[z][i][j] << endl;
                outFile4 << p2[z][i][j] << endl;
                outFile4 << p3[z][i][j] << endl;
            }
        }
    }
    outFile4.close();

    U.deleteZZXArray3(u2, packedRowCount, packedColCount, k);
    U.deleteZZXArray2(v2, packedRowCount, packedColCount);
    U.deleteZZXArray3(p1,kerN, packedRowCount,packedColCount);
    U.deleteZZXArray3(p2, kerN, packedRowCount,packedColCount);
    U.deleteZZXArray3(p3, kerN, packedRowCount,packedColCount);
    
}

void Cloud_cipher::homMultConv_Conv2(string outputcipherfile, string weightplainfile, const char *inputfile, int kerN, int kerD, int kerH, int kerW, int inH, int inW, int strid_H, int strid_W, int pool_strid_W, int pool_strid_H, string padding,ZZ_pX pm,  ZZ_pX p_temp1,  ZZ_pX p_temp2,  ZZ_pX p_temp3,  ZZ_pX p_temp4)
{   
    util U;
    timeval t1, t2;
  
    int outH = ceil((inH*1.0)/ strid_H);
    int outW = ceil((inW*1.0)/ strid_W);
      
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n/(real_nrow);
    int max_pool_num = max_packed_num / (pool_strid_H*pool_strid_W);
    int packed_num = max_pool_num * pool_strid_H*pool_strid_H; 
    int packed_col = ceil(double(1.0*datcols/packed_num));
    
    
    
    ZZX W[kerN][kerD];   
    ifstream inFile1;
    inFile1.open(weightplainfile, ios::in);
    if (!inFile1){
        cout << "Cannot open weight plain file in homMult!" << endl;
        exit(0);
    }

    for (int z1 = 0; z1 < kerN; z1++){
        for(int z2=0; z2 < kerD; z2++){
            inFile1 >> W[z1][z2];
        }
    }
    inFile1.close();
 

    //read input cipher
    ZZX ***u2 = U.newZZXArray3(kerD, packed_col, k);
    ZZX **v2 = U.newZZXArray2(kerD, packed_col);
    fstream fr1(inputfile);
    if (!fr1.is_open()){
        cout << "Cannot open input file in homMult!" << endl;
        exit(0);
    }
    for(int i = 0; i < kerD; i++){
        for (int j = 0; j < packed_col; j++){            
            for(int z = 0; z < k; z++){
                fr1 >> u2[i][j][z];
                U.decompress(u2[i][j][z], n, q, du);
            }
            fr1 >> v2[i][j];
            U.decompress(v2[i][j], n, q, dv);
        }
    }    
    fr1.close();
  
       
  

    gettimeofday(&t1,NULL); 

    ZZX*** p1 = U.newZZXArray3(kerN, kerD, packed_col); 
    ZZX*** p2 = U.newZZXArray3(kerN, kerD, packed_col);
    ZZX*** p3 = U.newZZXArray3(kerN, kerD, packed_col);
   
    // ZZ_p::init(q2);
    // ZZ_pX ppm;
    // SetCoeff(ppm,0,1);
    // SetCoeff(ppm,n,1);
    // ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
    for(int i = 0; i < kerN; i++){ 
        for(int j = 0; j < kerD; j++){           
            for(int z1 = 0; z1 < packed_col; z1++){                         
                U.Mult(p1[i][j][z1],pm,W[i][j],u2[j][z1][0],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
                U.Mult(p2[i][j][z1],pm,W[i][j],u2[j][z1][1],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
                U.Mult(p3[i][j][z1],pm,W[i][j],v2[j][z1],n,q,p_temp1,p_temp2,p_temp3,p_temp4);
            }
        }     
    }
    //对于通道1到通道kerD,每个通道的多项式的对应位置的系数相加
    ZZX**CsumMat1 = U.newZZXArray2(kerN, packed_col); 
    ZZX**CsumMat2 = U.newZZXArray2(kerN, packed_col); 
    ZZX**CsumMat3 = U.newZZXArray2(kerN, packed_col); 
  
    for(int i = 0; i < kerN; i++){
        for(int j1 = 0; j1 < packed_col; j1++){        
            ZZX value1 = ZZX(0);
            ZZX value2 = ZZX(0);
            ZZX value3 = ZZX(0);
            for(int z=0; z<kerD; z++){
                value1 = value1 + p1[i][z][j1];
                value2 = value2 + p2[i][z][j1];
                value3 = value3 + p3[i][z][j1];
            }
            U.Module(value1,n,q);
            U.Module(value2,n,q);
            U.Module(value3,n,q);      
            CsumMat1[i][j1] = value1;
            CsumMat2[i][j1] = value2;
            CsumMat3[i][j1] = value3;    
        }//end for j1
    }//end for i


    
    //     ZZX s[k];
    //     ifstream inFile22;
    //     inFile22.open("./file/key/sk.csv",ios::in); 
    //     if(!inFile22){
    //         cout<<"Cannot open sk in dec one ploy!"<<endl;
    //         exit(0);
    //     }   
    //         for(int i=0;i<k;i++){
    //         inFile22>>s[i];            
    //     }    
    //     // cout<<endl;
    //     inFile22.close();
    //     ZZX h5[kerN][packed_col];
    //     ZZX s0u10w1[kerN][packed_col],
    //     s1u11w1[kerN][packed_col];
    //     for(int i=0; i<kerN; i++){
    //     for(int j = 0; j < packed_col; j++){
          
    //         U.Mult(s0u10w1[i][j],ppm,s[0],CsumMat1[i][j],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    //         U.Mult(s1u11w1[i][j],ppm,s[1],CsumMat2[i][j],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    //         h5[i][j]=0;
    //         h5[i][j]=CsumMat3[i][j]-s0u10w1[i][j]-s1u11w1[i][j];  
    //         }
    //     } 
    // ZZ two_d;
    // RR two_d_rr;
    // power2(two_d,dp);
    // two_d_rr=to_RR(two_d);
    // RR q_rr;
    // q_rr = to_RR(q);
    // ZZ q2d;
    // q2d=RoundToZZ(q_rr/two_d_rr);

    // RR two_dp_rr;
    // power2(two_dp_rr,dp);
    // RR q2_rr;
    // ZZ two_d2;
    // ZZ two_d3;
    // two_d2=two_d/2;
    // two_d3=two_d2*(-1);
    // q2_rr = q_rr;


    // RR coeff_rr;
    // ZZ coeff_zz;
    // for(int i=0; i<kerN; i++){
    //     for(int j=0; j<packed_col; j++){            
    //         for(int j3=0; j3<n; j3++){
    //             GetCoeff(coeff_zz, h5[i][j], j3);
    //             coeff_rr=to_RR(coeff_zz);
    //             coeff_zz=RoundToZZ(coeff_rr*two_dp_rr/q2_rr)%two_d;
    //             if(coeff_zz>=two_d2){
    //                 coeff_zz= coeff_zz-two_d;
    //             }
    //             else if(coeff_zz<=two_d3){
    //                 coeff_zz=coeff_zz+two_d;
    //             }
    //             SetCoeff(h5[i][j],j3,coeff_zz);
    //         }
    //     }
    // }
    // cout<<"jieguowei"<<endl;
    // for(int i=0; i<kerN; i++){
    //     for(int j1=0; j1<packed_col; j1++){
    //             for(int z1=0;z1<packed_num;z1++){
    //                 cout<<coeff(h5[i][j1],25+26*z1)<<" ";
    //             }
    //         }
    //         cout<<endl;
    //     }
    //     cout<<"&&&&&&&&&&"<<endl;
 
    // cout<<"end@@@@@@@@@@@@@@@"<<endl;



    gettimeofday(&t2,NULL);
    double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout<<"conv2:Homorphic eval time is:"<<deltaT/1000<<"ms"<<endl<<endl;
    
    clock_t startTime11;
    clock_t endTime11;
    startTime11=clock();

    ofstream outFile3;
    outFile3.open(outputcipherfile, ios::ate);   

    for(int i = 0; i < kerN; i++){
        for (int j1 = 0; j1 < packed_col; j1++){           
            outFile3<<CsumMat1[i][j1]<<endl;
            outFile3<<CsumMat2[i][j1]<<endl;
            outFile3<<CsumMat3[i][j1]<<endl;                 
        }
    }
    outFile3.close();      
    
    // cout<<"Homomorphic res store time is "<<(double)(endTime11-startTime11)/CLOCKS_PER_SEC<<"s"<<endl;
    U.deleteZZXArray3(u2, kerD, packed_col, k);
    U.deleteZZXArray2(v2, kerD, packed_col);
    U.deleteZZXArray3(p1,kerN,kerD, packed_col);
    U.deleteZZXArray3(p2, kerN,kerD, packed_col);
    U.deleteZZXArray3(p3, kerN,kerD, packed_col);  
    U.deleteZZXArray2(CsumMat1,kerN, packed_col); 
    U.deleteZZXArray2(CsumMat2,kerN, packed_col);
    U.deleteZZXArray2(CsumMat3,kerN, packed_col);  
        
}



//云服务器加密四维权重
void Cloud_cipher::weight_plainGen(int kerN,int kerD,int kerH, int kerW, string orignalweight,string orignalbias,string plainfile) 
{      
    util U;
    Cloud_util cdu;
    //读取权重然后round
    int ****w_conv2 = U.newArray4(kerN,kerD,kerH,kerW);
    cdu.readweights2D_4D(w_conv2,kerH,kerW,kerD,kerN,orignalweight);
    //读取偏置b,然后round
    int b_conv2[kerN];
    cdu.readbias(b_conv2,kerN,orignalbias); 

    //将W四维变三维
    int ***ww_conv2=U.newArray3(kerN,kerD,kerH*kerW);         
    for(int u=0;u<kerN;u++) 
        for(int v=0;v<kerD;v++)
            for(int i=0;i<kerH;i++)
                for(int j=0;j<kerW;j++)
                    ww_conv2[u][v][i*kerW+j] = w_conv2[u][v][i][j];                 
    // 三维数组ww_conv2变二维多项式W,同时装配上b
    ZZX **W = U.newZZXArray2(kerN,kerD);
    cdu.matTo2Dpoly(W,ww_conv2,b_conv2,kerN,kerD,kerH,kerW);  
    //  多项式进行存储
    ofstream outFile2;
    outFile2.open(plainfile, ios::ate);
    for (int z = 0; z < kerN; z++){	
        for(int j = 0; j < kerD; j++){         	
		    outFile2 << W[z][j] << endl;
        }
	}
    outFile2.close();
    //释放内存
    U.deleteArray3(ww_conv2,kerN,kerD,kerH*kerW);
    U.deleteZZXArray2(W,kerN,kerD);
    U.deleteArray4(w_conv2,kerN,kerD,kerH,kerW); 
}

/*void Cloud_cipher::weight_plainGen2(int kerN,int kerD,int kerH, int kerW, string orignalweight,string orignalbias,string plainfile) 
{      
  
    Cloud_util cdu;
    //读取权重然后round
    int ****w_conv2= U.newArray4(kerN,kerD,kerH,kerW);
    cdu.readweights2D_4D(w_conv2,kerH,kerW,kerD,kerN,orignalweight);
    //读取偏置b,然后round
    int b_conv2[kerN];
    cdu.readbias(b_conv2,kerN,orignalbias); 
    //将W四维变三维
    int ***ww_conv2=U.newArray3(kerN,kerD,kerH*kerW);         
    for(int u=0;u<kerN;u++) 
        for(int v=0;v<kerD;v++)
            for(int i=0;i<kerH;i++)
                for(int j=0;j<kerW;j++)
                    ww_conv2[u][v][i*kerW+j]=w_conv2[u][v][i][j];                 
    // 三维数组ww_conv2变二维多项式W,同时装配上b
    // ZZX **W=U.newZZXArray2(kerN,kerD);

    // cdu.matTo2Dpoly(W,ww_conv2,b_conv2,kerN,kerD,kerH,kerW);  
    
    
     int **w2Dconv2=U.newArray2(kerN,kerD*kerH*kerW);       
    for(int u=0;u<kerN;u++) {
        for(int v=0;v<kerD;v++){
            for(int j=0;j<kerH*kerW;j++){
                //cout<<w3Dconv2[u][v][j]<<" ";
                w2Dconv2[u][v*kerH*kerW+j]=ww_conv2[u][v][j];
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

    cout<<"装配上b的w"<<endl;
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
    
    // // 二维矩阵进行存储
    ofstream outFile2;
    outFile2.open(plainfile,ios::ate);
    for (int z=0;z<kerN;z++){	
        for(int j=0;j<kerD*kerH*kerW+1;j++){	
		    outFile2<<wb2Dconv2[z][j]<<endl;
           }
	}
    outFile2.close();
    //释放内存
    U.deleteArray3(w_conv2,kerN,kerD,kerH*kerW);
    U.deleteArray3(ww_conv2,kerN,kerD,kerH*kerW);
    U.deleteArray2(w2Dconv2,kerN,kerD*kerH*kerW);
}*/



void Cloud_cipher::homMultHaveAdd(string outputcipherfile,string weightplainfile, const char* inputfile,int kerN,int kerD,int kerH,int kerW,int inH,int inW,int strid_H,int strid_W,string padding)
{   util U;
    timeval t1, t2;
    if (padding=="SAME"){
        int outH=ceil(inH/strid_H);
        int outW=ceil(inW/strid_W);
        int packed_num1 = n/(kerH*kerW+1);
        int packed_col1 =ceil(float(outH*outW)/packed_num1);

        int packed_num2 = n/(kerH*kerW);
        int packed_col2 =ceil(float(outH*outW)/packed_num2);
        // cout<<"packed_num1="<<packed_num1<<endl;
        // cout<<"packed_col1="<<packed_col1<<endl;
        // cout<<"packed_num2="<<packed_num2<<endl;
        // cout<<"packed_col2="<<packed_col2<<endl;

        //read weight plaintext
        ZZX W[kerN][kerD]; 
        ifstream inFile1;
        inFile1.open(weightplainfile,ios::in);
        if(!inFile1){
            cout<<"Cannot open weight plain file in homMult!"<<endl;
            exit(0);
        }

        for(int z=0;z<kerN;z++){
            for(int i=0;i<kerD;i++){          
                    inFile1>>W[z][i];   
                }          
        } 
        inFile1.close(); 

        
        //read input cipher
        ZZX ***u2=U.newZZXArray3(kerD,packed_col1,k);
        ZZX **v2=U.newZZXArray2(kerD,packed_col1); 

        // long cu2[kerD][packed_col1][k][n];
        // long cv2[kerD][packed_col1][n];

        
        // FILE *pFile2=fopen(inputfile,"rb");
        // fread(cu2,sizeof(cu2[0][0][0][0]),kerD*packed_col1*k*n,pFile2);
        // fread(cv2,sizeof(cv2[0][0][0]),kerD*packed_col1*n,pFile2);
        // fclose(pFile2);

        // cout<<"read test======================="<<endl;
        // cout<<cu2[0][0][0][0]<<"   "<<cu2[0][1][1][4]<<endl;
        // cout<<cv2[0][0][0]<<"   "<<cv2[0][1][3]<<endl;
        // cout<<"****************************************"<<endl;
        // for(int z=0;z<kerD;z++){
        //     for(int i=0;i<packed_col1;i++){
        //         for (int j=0;j<k;j++){
        //             for(int u=0;u<n;u++){
        //                 SetCoeff(u2[z][i][j],u,conv<ZZ>(cu2[z][i][j][u]));        
        //             }     
        //         } 
        //         for(int v=0;v<n;v++){
        //             SetCoeff(v2[z][i],v,conv<ZZ>(cv2[z][i][v]));   
        //         }
        //     }      
        // } 

        // for(int z=0;z<kerD;z++){
        //     for(int i=0;i<packed_col1;i++){
        //         for (int j=0;j<k;j++){                        
        //             U.decompress(u2[z][i][j],n,q,du);        
        //         }    
        //         U.decompress(v2[z][i],n,q,dv);        
        //     }
        // } 

        // //read input cipher
        // ZZX ***u2=U.newZZXArray3(kerD,packed_col1,k);
        // ZZX **v2=U.newZZXArray2(kerD,packed_col1);  
    
        // ifstream inFile2;
        // inFile2.open(inputfile,ios::in);  
        // if(!inFile2){
        //     cout<<"Cannot open input cipher file in homMult!"<<endl;
        //     exit(0);
        // }
        // for(int z=0;z<kerD;z++)
        //     for(int i=0;i<packed_col1;i++){
        //         for (int j=0;j<k;j++){
        //             inFile2>>u2[z][i][j];       
        //             U.decompress(u2[z][i][j],n,q,du);        
        //         }     
        //         inFile2>>v2[z][i]; 
        //         U.decompress(v2[z][i],n,q,dv);        
        // } 
        // inFile2.close();

        // cout<<"recover=============================="<<endl;
        // cout<<coeff(u2[0][0][0],0)<<"   "<<coeff(u2[0][1][1],4)<<endl;
        // cout<<coeff(v2[0][0],0)<<"  "<<coeff(v2[0][1],3)<<endl;
        // cout<<"**************************************"<<endl;
        ZZX p1[kerD][kerN][packed_col1],p2[kerD][kerN][packed_col1],p3[kerD][kerN][packed_col1];    
       
        ZZ_p::init(q2);
        ZZ_pX ppm;
        SetCoeff(ppm,0,1);
        SetCoeff(ppm,n,1);
        ZZ_pXModulus pm(ppm);
        ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
        gettimeofday(&t1,NULL);
        for(int j=0; j<kerD; j++){
            for(int i=0; i<kerN;i++){      
                if(j==0){
                    for(int u=0; u<packed_col1; u++){
                        // U.Mult(v1v2[j][i][u],pm,v1[i][0],v2[0][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        // U.Mult(u11u21[j][i][u],pm,u1[i][0][1],u2[0][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p1[j][i][u],pm,W[i][0],v2[0][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p2[j][i][u],pm,W[i][0],u2[0][u][0],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p3[j][i][u],pm,W[i][0],u2[0][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    }//end for u
                }
                else{
                        for(int u=0; u<packed_col2; u++){
                            // U.Mult(v1v2[j][i][u],pm,v1[i][j],v2[j][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                            // U.Mult(u11u21[j][i][u],pm,u1[i][j][1],u2[j][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                            U.Mult(p1[j][i][u],pm,W[i][j],v2[j][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                            U.Mult(p2[j][i][u],pm,W[i][j],u2[j][u][0],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                            U.Mult(p3[j][i][u],pm,W[i][j],u2[j][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        }
                }
            }// end for i
        }//end for j
        //对于通道2到通道kerD,每个通道的多项式的对应位置的系数相加
        ZZX CsumMat1[kerN][packed_col2];
        ZZX CsumMat2[kerN][packed_col2];
        ZZX CsumMat3[kerN][packed_col2];
        // ZZX CsumMat9[kerN][packed_col2];

        for(int i=0; i<kerN; i++){
            for(int j=0; j<packed_col2; j++){
                ZZX value1 = ZZX(0);
                ZZX value2 = ZZX(0);
                ZZX value3 = ZZX(0);

                for(int z=1; z<kerD; z++){
                    value1 = value1 + p1[z][i][j];
                    value2 = value2 + p2[z][i][j];
                    value3 = value3 + p3[z][i][j];

                }
                CsumMat1[i][j] = value1;
                CsumMat2[i][j] = value2;
                CsumMat3[i][j] = value3;
        
            }
        }
        gettimeofday(&t2,NULL);
        double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        cout<<"1:Homorphic eval time is:"<<deltaT/1000<<"ms"<<endl<<endl;
        
        clock_t startTime11;
        clock_t endTime11;
        startTime11=clock();

        if(kerD==1){
            ofstream outFile4;
            outFile4.open(outputcipherfile,ios::ate);   
            for(int i=0;i<kerN;i++){
                for (int j=0;j<packed_col1;j++){
                    outFile4<<p1[0][i][j]<<endl;    
                    outFile4<<p2[0][i][j]<<endl;    
                    outFile4<<p3[0][i][j]<<endl;                  
                }
            }
            outFile4.close();
        }

        else{
            ofstream outFile3;
            outFile3.open(outputcipherfile,ios::ate);   
            for(int i=0;i<kerN;i++){
                for (int j=0;j<packed_col1;j++){
                    outFile3<<p1[0][i][j]<<endl;    
                    outFile3<<p2[0][i][j]<<endl;    
                    outFile3<<p3[0][i][j]<<endl; 
                        
                }
            }
            for(int i=0; i<kerN; i++){
                for(int j=0; j<packed_col2; j++){
                    outFile3<<CsumMat1[i][j]<<endl;
                    outFile3<<CsumMat2[i][j]<<endl;
                    outFile3<<CsumMat3[i][j]<<endl;
        
                }
            }
            outFile3.close(); 
        } 
        endTime11=clock();
        cout<<"Homomorphic res store time is "<<(double)(endTime11-startTime11)/CLOCKS_PER_SEC<<"s"<<endl;
        U.deleteZZXArray3(u2,kerD,packed_col1,k);
        U.deleteZZXArray2(v2,kerD,packed_col1);
    }

    
    else if(padding=="VALID"){
        int outH=ceil(float(inH-kerH+1)/float(strid_H));
        int outW=ceil(float(inW-kerW+1)/float(strid_W));
        cout<<"outH="<<outH<<endl;
        int packed_num1 = n/(kerH*kerW+1);
        int packed_col1 =ceil(double(outH*outW)/packed_num1);
        int packed_num2 = n/(kerH*kerW);
        int packed_col2 =ceil(double(outH*outW)/packed_num2);
        cout<<"********************************"<<endl;
        cout<<"packed_col1="<<packed_col1<<endl;
        cout<<"packed_col2="<<packed_col2<<endl;
        cout<<"********************************"<<endl;

    //read weight plaintext
    ZZX W[kerN][kerD]; 
    ifstream inFile1;
    inFile1.open(weightplainfile,ios::in);
    if(!inFile1){
        cout<<"Cannot open weight plain file in homMult!"<<endl;
        exit(0);
    }

    for(int z=0;z<kerN;z++){
        for(int i=0;i<kerD;i++){          
                inFile1>>W[z][i];   
            }          
    } 
    inFile1.close(); 
 cout<<"read weight plaintext!"<<endl;

    
    //read input cipher
    ZZX ***u2=U.newZZXArray3(kerD,packed_col1,k);
    ZZX **v2=U.newZZXArray2(kerD,packed_col1); 

    long cu2[kerD][packed_col1][k][n];
    long cv2[kerD][packed_col1][n];

    
    FILE *pFile2=fopen(inputfile,"rb");
    fread(cu2,sizeof(cu2[0][0][0][0]),kerD*packed_col1*k*n,pFile2);
    fread(cv2,sizeof(cv2[0][0][0]),kerD*packed_col1*n,pFile2);
    fclose(pFile2);

    for(int z=0;z<kerD;z++){
        for(int i=0;i<packed_col1;i++){
            for (int j=0;j<k;j++){
                for(int u=0;u<n;u++){
                    SetCoeff(u2[z][i][j],u,conv<ZZ>(cu2[z][i][j][u]));   
                } 
                 
            } 
            for(int v=0;v<n;v++){
                SetCoeff(v2[z][i],v,conv<ZZ>(cv2[z][i][v]));                    
            }
        }      
	} 

    for(int z=0;z<kerD;z++){
            for(int i=0;i<packed_col1;i++){
                for (int j=0;j<k;j++){                        
                    U.decompress(u2[z][i][j],n,q,du);        
                }    
                U.decompress(v2[z][i],n,q,dv);        
            }
    }  
  
    ZZX p1[kerD][kerN][packed_col1],p2[kerD][kerN][packed_col1],p3[kerD][kerN][packed_col1];    
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
   
    
    gettimeofday(&t1,NULL);
    for(int j=0; j<kerD; j++){
        for(int i=0; i<kerN;i++){      
            if(j==0){
                for(int u=0; u<packed_col1; u++){
                    // U.Mult(v1v2[j][i][u],pm,v1[i][0],v2[0][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    U.Mult(p1[j][i][u],pm,W[i][0],v2[0][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    U.Mult(p2[j][i][u],pm,W[i][0],u2[0][u][0],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    U.Mult(p3[j][i][u],pm,W[i][0],u2[0][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    
                }//end for u
            }
            else{
                    for(int u=0; u<packed_col2; u++){
                        // U.Mult(u11u21[j][i][u],pm,u1[i][j][1],u2[j][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p1[j][i][u],pm,W[i][j],v2[j][u],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p2[j][i][u],pm,W[i][j],u2[j][u][0],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                        U.Mult(p3[j][i][u],pm,W[i][j],u2[j][u][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
                    }
            }
        }// end for i
    }//end for j


    //对于通道2到通道kerD,每个通道的多项式的对应位置的系数相加
    ZZX CsumMat1[kerN][packed_col2];
    ZZX CsumMat2[kerN][packed_col2];
    ZZX CsumMat3[kerN][packed_col2];

    for(int i=0; i<kerN; i++){
        for(int j=0; j<packed_col2; j++){
            ZZX value1 = ZZX(0);
            ZZX value2 = ZZX(0);
            ZZX value3 = ZZX(0);

            for(int z=1; z<kerD; z++){
                value1 = value1 + p1[z][i][j];
                value2 = value2 + p2[z][i][j];
                value3 = value3 + p3[z][i][j];

            }
            CsumMat1[i][j] = value1;
            CsumMat2[i][j] = value2;
            CsumMat3[i][j] = value3;
        }
    }
   
   
    gettimeofday(&t2,NULL);
    double deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout<<"1:Homorphic eval time is:"<<deltaT/1000<<"ms"<<endl<<endl;
    
    clock_t startTime11;
    clock_t endTime11;
    startTime11=clock();
    
    if(kerD==1){
    ofstream outFile4;
    outFile4.open(outputcipherfile,ios::ate);   
    for(int i=0;i<kerN;i++){
        for (int j=0;j<packed_col1;j++){
            outFile4<<p1[0][i][j]<<endl;    
            outFile4<<p2[0][i][j]<<endl;    
            outFile4<<p3[0][i][j]<<endl;                  
        }
    }
    outFile4.close();
   }
   else{
        ofstream outFile3;
        outFile3.open(outputcipherfile,ios::ate);   
        for(int i=0;i<kerN;i++){
            for (int j=0;j<packed_col1;j++){
                outFile3<<p1[0][i][j]<<endl;    
                outFile3<<p2[0][i][j]<<endl;    
                outFile3<<p3[0][i][j]<<endl; 
                    
            }
        }
        for(int i=0; i<kerN; i++){
            for(int j=0; j<packed_col2; j++){
                outFile3<<CsumMat1[i][j]<<endl;
                outFile3<<CsumMat2[i][j]<<endl;
                outFile3<<CsumMat3[i][j]<<endl;
    
            }
        }
        outFile3.close();  
   }
    endTime11=clock();
    cout<<"Homomorphic res store time is "<<(double)(endTime11-startTime11)/CLOCKS_PER_SEC<<"s"<<endl;
    U.deleteZZXArray3(u2,kerD,packed_col1,k);
    U.deleteZZXArray2(v2,kerD,packed_col1);
 }

}




void Cloud_cipher::dense_weight_plainGen(int den_num,int fc_w,string w_dense,string b_dense,string fc_weightpliantex)
{    util U;  
    //input is sizeof(fci_w),output size is sizeof(den_num1)
    // int packednum=ceil(float(fc1_w+1)/n); 
    int packednum=ceil(1.0*fc_w/(n-1));  
    enc_param usp;
    Cloud_util cdu;
    
    int **w_den=new int *[den_num];
    for (int i=0;i<den_num;i++){
        w_den[i]=new int[n*packednum];
    }  
    U.initArray2(w_den,den_num,n*packednum,0);
   
    cdu.readweights(w_den,den_num,fc_w,w_dense);  
    int b[den_num];
    cdu.readbias(b,den_num,b_dense); 
    
    ZZX **W_DEN=new ZZX *[den_num]();
    for (int i=0;i<den_num;i++){
        W_DEN[i]=new ZZX[packednum]();
    }  
    for(int i=0;i<den_num;i++){
        W_DEN[i][0].SetLength(n);
        SetCoeff(W_DEN[i][0],0,b[i]);
        for(int z=0;z<n-1;z++)
        SetCoeff(W_DEN[i][0],z+1,w_den[i][z]);
    }
    for(int i=0;i<den_num;i++){
        for(int j=1;j<packednum;j++){
            W_DEN[i][j].SetLength(n);
            SetCoeff(W_DEN[i][j],0,0);
            for(int z=0;z<n-1;z++)
            SetCoeff(W_DEN[i][j],z+1,w_den[i][j*(n-1)+z]);
        }
    }
    ofstream outFile2;
    outFile2.open(fc_weightpliantex,ios::ate);
    for (int z=0;z<den_num;z++){	
        for(int j=0;j<packednum;j++){	
		    outFile2<<W_DEN[z][j]<<endl;
           }
	}
    outFile2.close();
    U.deleteArray2(w_den,den_num,fc_w);
}

// void Cloud_cipher::dense_weight_plainGen(int den_num1,int fc1_w,string w_dense1,string b_dense1,string fc_weightpliantex,helib:: EncryptedArray ea)
// {      
//    //input is sizeof(fc1_w),output size is sizeof(den_num1)
//     // cout<<"fc1_w="<<fc1_w<<endl;
//     int packednum=ceil(float(fc1_w+1)/n); 
//     cloud_enc_param usp;
//     Cloud_util cdu;
//     util U;
//     int **w_den1=new int *[den_num1];
//     for (int i=0;i<den_num1;i++){
//         w_den1[i]=new int[fc1_w];
//     }  
   
//     cdu.readweights(w_den1,den_num1,fc1_w,w_dense1);  
//     int b_fc1[den_num1];
//     cdu.readbias(b_fc1,den_num1,b_dense1); 
    
//     int slots=ea.size();
 
//     int **wb_fc1=new int *[den_num1];
//     for (int i=0;i<den_num1;i++){
//         wb_fc1[i]=new int[fc1_w+1];
//     }   
     
//     for(int i=0;i<den_num1;i++){
//         wb_fc1[i][0]=b_fc1[i];
//         for(int j=1;j<fc1_w+1;j++){   
//             wb_fc1[i][j]=w_den1[i][j-1];           
//         }
//     } 
//     vector <vector <long>> CC;
//     // int ** CC=new int *[den_num1];

//     // for(int i=0;i<den_num1;i++)
//     //     CC[i]=new int[packednum*n];
//     //将权值装配入多项式，不足位补0.
    
//     for(int i=0;i<den_num1;i++){
//         for(int j=0;j<packednum*n;j++){
//             CC[i][j]=0;      
//             }
//     }   
//     for(int i=0;i<den_num1;i++){
//         for(int j=0;j<fc1_w+1;j++){
//             CC[i][j]=wb_fc1[i][j];
//         }
//     }

//     ZZX enc_poly[den_num1][packednum];
//     for(int i=0;i<den_num1;i++){
//         for(int j=0;j<packednum;j++){
//             for(int z=0;z<n;z++){
//                 ea.encode(enc_poly[i][j],CC[i][j*n+z]);
//             }
//         }
//     }
    
//     ofstream outFile2;
//     outFile2.open(fc_weightpliantex,ios::ate);
//     for (int z=0;z<den_num1;z++){	
//         for(int j=0;j<packednum;j++){	
// 		    outFile2<<W_DEN[z][j]<<endl;
//            }
// 	}
//     outFile2.close();
//     U.deleteArray2(w_den,den_num,fc_w);
    

// }



void Cloud_cipher::homMultHaveAdd(string outputfile,string weightfile, string inputfile,int dennum,int fc1_w)
 {  util U;
    int a=ceil(1.0*fc1_w/(n-1));
    // cout<<"In homomorphic con packednum="<<a<<endl;
    //read weightplaintex.csv
    ZZX W_DEN[dennum][a];
    ifstream inFile1;
    inFile1.open(weightfile,ios::in);
    if(!inFile1){
        cout<<"Cannot open cipher weightfile in homomorphic com!"<<endl;
        exit(0);
    }
    for(int z=0;z<dennum;z++){
        for(int i=0;i<a;i++){
            inFile1>>W_DEN[z][i];
        }     
    }     
    inFile1.close(); 

    //read inputcipher.csv
   	ZZX** u2= new ZZX*[a]();//input M has 13(outH*outW) cipher,and each cipher is vector and each element is a poly with n degree
    for (int i = 0; i < a; i++) {
        u2[i] = new ZZX[k]();
    }   
    ZZX v2[a];   
    ifstream inFile2;
    inFile2.open(inputfile,ios::in);  
   if(!inFile2){
        cout<<"Cannot open cipher inputfile in homomorphic com!"<<endl;
        exit(0);
    } 
    // timeval s1,s2;
    // gettimeofday(&s1,NULL);     
    for(int i=0;i<a;i++){
        for (int j=0;j<k;j++){
            inFile2>>u2[i][j];       
            U.decompress(u2[i][j],n,q,du);        
        }     
        inFile2>>v2[i]; 
        U.decompress(v2[i],n,q,dv);        
	} 
    inFile2.close();
   
    ZZ_p::init(q2);  
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,n,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4; 
    //dennum=120,a=4
    ZZX p1[a*dennum], p2[a*dennum], p3[a*dennum];
    int z=0;            
    for(int u=0;u<dennum;u++){ 
        for(int j=0;j<a;j++){
            U.Mult(p1[z],ppm,W_DEN[u][j],v2[j],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
            U.Mult(p2[z],ppm,W_DEN[u][j],u2[j][0],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
            U.Mult(p3[z],ppm,W_DEN[u][j],u2[j][1],n,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
            z=z+1;            
        }    
    }    

    ZZX *sum1=new ZZX[dennum];
    ZZX *sum2=new ZZX[dennum];
    ZZX *sum3=new ZZX[dennum];
    int zz=0;   
    for(int i=0;i<dennum;i++)
    {   sum1[i]=ZZX(0);
        sum2[i]=ZZX(0);
        sum3[i]=ZZX(0);
        while(zz<a*(i+1)){
            sum1[i]=sum1[i]+p1[zz];
            sum2[i]=sum2[i]+p2[zz];
            sum3[i]=sum3[i]+p3[zz];
            zz++;
        }
      
    }

    // gettimeofday(&s2,NULL);
    // double interval = (s2.tv_sec - s1.tv_sec) * 1000000 + s2.tv_usec - s1.tv_usec;;
    // cout<<"FC Homomorphic eval time is "<<interval/1000<<"ms"<<endl; 

    //store cloud compute result in outputfile    
    ofstream outFile3;
    outFile3.open(outputfile,ios::ate);  
    for(int i=0;i<dennum;i++){
            outFile3<<sum1[i]<<endl;    
            outFile3<<sum2[i]<<endl;    
            outFile3<<sum3[i]<<endl;                   
        }
    outFile3.close();
    delete[] sum1,sum2,sum3; 
    U.deleteZZXArray2(u2,a,k);
}



void Cloud_cipher::cipherAdd_Noise(int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, string outcipherfile1, string cipherfile1,string cipherfile2)
{
    util U;
    Cloud_util cdu; 
    int outH = ceil((inH*1.0)/strid_H);
    int outW = ceil((inW*1.0)/strid_W);      
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n/(real_nrow);
    int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
    int packed_num = max_pool_num * pool_stridH * pool_stridH; 
    int packed_col = ceil(double(1.0*datcols/packed_num));
    cout<<"packed_col="<<packed_col<<endl;
    
    //求同态加并加上噪声：w[g]+[R]+[w(x-g)]
    //打开密文文件 w[g],无需解压    
    ZZX ***u1 = U.newZZXArray3(kerN, packed_col, k);
    ZZX **v1 = U.newZZXArray2(kerN, packed_col);   
    ifstream inFile1;
    inFile1.open(cipherfile1, ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   
    for(int z1 = 0; z1 < kerN; z1++){
        for (int i = 0; i < packed_col; i++){
            for (int j = 0; j < k; j++){                                     
                inFile1 >> u1[z1][i][j];                 
            }
                inFile1 >> v1[z1][i];                     
        }              
    }    
    inFile1.close();

    //打开第二个密文文件[w(x-g)]
    ZZX ***u2 = U.newZZXArray3(kerN, packed_col, k);
    ZZX **v2 = U.newZZXArray2(kerN, packed_col);     
    ifstream inFile2;
    inFile2.open(cipherfile2, ios::in); 
    if(!inFile2){
        cout<<"Cannot open cipher2 file!\n";
        exit(0);
    } 
    for(int z1 = 0; z1 < kerN; z1++){  
        for (int i = 0; i < packed_col; i++){         
            for (int z = 0; z < k; z++){               
                inFile2 >> u2[z1][i][z];  
                U.decompress(u2[z1][i][z], n, q, du) ;                  
            }   
            inFile2 >> v2[z1][i]; 
            U.decompress(v2[z1][i], n, q, dv); 
        }
    }
    inFile2.close();
   
    ofstream outFile1;
    outFile1.open(outcipherfile1, ios::ate); 
    if(!outFile1){
        cout << "Cannot open Gencipher file!\n";
        exit(0);
    } 
    for(int z1 = 0; z1 < kerN; z1++){  //直接加上用于maxpool的噪声
        for (int i = 0; i < packed_col; i++){ 
            u2[z1][i][0] = u1[z1][i][0] + u2[z1][i][0] ;
            u2[z1][i][1] = u1[z1][i][1] + u2[z1][i][1] ;
            v2[z1][i] = v1[z1][i] + v2[z1][i] ;             
            outFile1 << u2[z1][i][0] << endl;   
            outFile1 << u2[z1][i][1] << endl;  
            outFile1 << v2[z1][i] << endl;             
        }
    }
    outFile1.close(); 
    U.deleteZZXArray3(u1,kerN, packed_col,k); 
    U.deleteZZXArray2(v1,kerN, packed_col);
    U.deleteZZXArray3(u2,kerN, packed_col,k); 
    U.deleteZZXArray2(v2,kerN, packed_col);   
      
}



void Cloud_cipher::plainEnc_Add(ZZX**R, int**C,int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, const char *outcipherfile1, const char *cipherfile1,string padding,ZZX *t, ZZX **A,ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp)
{  
    util U;
    if(padding=="SAME"){ 
        int outH = ceil((inH * 1.0) / strid_H);
        int outW = ceil((inW * 1.0) / strid_W);      
        int datrows = kerH * kerW;
        int datcols = outH * outW; 
        int real_nrow = datrows + 1;            
        int max_packed_num = n / (real_nrow);
        int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
        int packed_num = max_pool_num * pool_stridH * pool_stridH; 
        int packed_col = ceil(double(1.0*datcols/packed_num));
        
        timeval p1, p2;
        gettimeofday(&p1, NULL);        
        //所有通道按照方法2装配为多项式 
        //outH*ouW列每个列都添加了1，不用倒排,再装配为多项式      
    timeval test1,test2;
    gettimeofday(&test1,NULL);
        int pcol = outW;      
        Cloud_util cdu;
        ZZX **polys = U.newZZXArray2(kerN, packed_col);  
        for (int i = 0; i < kerN; i++){
            cdu.packingbyPoolWin(kerH, kerW, polys[i], C[i], datcols, pcol, pool_stridH, pool_stridW, max_pool_num);          
        }
    gettimeofday(&test2,NULL);
    double deltaTest1 = (test2.tv_sec - test1.tv_sec) * 1000000 
    + test2.tv_usec - test1.tv_usec;
    cout << "PoolWin计算总时间为:" << deltaTest1/1000<<"ms"<< endl;

  
    // for(int i = 0; i < kerN; i++){
    //      for(int j = 0; j < packed_col; j++){                    
    //         polys[i][j] += R[i][j];   
    //     }
    // }
    timeval test3,test4;
    gettimeofday(&test3,NULL);
        ZZX ***u2 = U.newZZXArray3(kerN, packed_col, k);
        ZZX **v2 = U.newZZXArray2(kerN, packed_col);
        Enc enc;  
        for(int z = 0; z < kerN; z++){
            for(int i = 0; i < packed_col; i++){  
                polys[z][i] += R[z][i];              
                enc.encrypt(polys[z][i],u2[z][i],v2[z][i],n,q,k,q2d,usp[0][z][i].r1, A, usp[0][z][i].e11,usp[0][z][i].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 
            }
        }
         gettimeofday(&test4,NULL);
    double deltaTest2 = (test4.tv_sec - test3.tv_sec) * 1000000 
    + test4.tv_usec - test3.tv_usec;
    cout << "encrypt计算总时间为:" << deltaTest2/1000<<"ms"<< endl;

    //求同态加并加上噪声：w[g]+[R]+[w(x-g)]
    //打开密文文件 w[g],无需解压    
    ZZX ***u1 = U.newZZXArray3(kerN, packed_col, k);
    ZZX **v1 = U.newZZXArray2(kerN, packed_col);   
    ifstream inFile1;
    inFile1.open(cipherfile1, ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   
    for(int z1 = 0; z1 < kerN; z1++){
        for (int i = 0; i < packed_col; i++){
            for (int j = 0; j < k; j++){                                     
                inFile1 >> u1[z1][i][j];                 
            }
                inFile1 >> v1[z1][i];                     
        }              
    }    
    inFile1.close();

  timeval test5,test6;
    gettimeofday(&test5, NULL);
    for(int z1 = 0; z1 < kerN; z1++){  
        for (int i = 0; i < packed_col; i++){
            U.decompress(u2[z1][i][0], n, q, du); 
            U.decompress(u2[z1][i][1], n, q, du);
            U.decompress(v2[z1][i], n, q, dv);     
            u2[z1][i][0] = u1[z1][i][0] + u2[z1][i][0] ;
            u2[z1][i][1] = u1[z1][i][1] + u2[z1][i][1] ;
            v2[z1][i] = v1[z1][i] + v2[z1][i] ;                          
        }
    }
    gettimeofday(&test6, NULL);
    double deltaTest3 = (test6.tv_sec - test5.tv_sec) * 1000000 
    + test6.tv_usec - test5.tv_usec;
    cout << "同态加计算时间为:" << deltaTest3/1000<<"ms"<< endl; 

    ofstream outFile1;
    outFile1.open(outcipherfile1, ios::ate); 
    if(!outFile1){
        cout << "Cannot open Gencipher file!\n";
        exit(0);
    } 
    for(int z1 = 0; z1 < kerN; z1++){  
        for (int i = 0; i < packed_col; i++){             
            outFile1 << u2[z1][i][0] << endl;   
            outFile1 << u2[z1][i][1] << endl;  
            outFile1 << v2[z1][i] << endl;             
        }
    }
    outFile1.close(); 
    U.deleteZZXArray3(u1,kerN, packed_col,k); 
    U.deleteZZXArray2(v1,kerN, packed_col);
    U.deleteZZXArray3(u2,kerN, packed_col,k); 
    U.deleteZZXArray2(v2,kerN, packed_col); 
  }
}


//云端逐个整数加密用于FC层
// void Cloud_cipher::dense_input_cipherGen(int* C, int lenth_C,string fc_inputcipherfile)
// {   
         

//     ZZX M[lenth_C];//for f1:13 ploy and each poly's degree n is 256
// 	for(int i=0;i<lenth_C;i++){
//         M[i].SetLength(n);
//         SetCoeff(M[i],0,C[i]);
       
//     }
//     //  gettimeofday(&p2, NULL);
//     // double deltaT = (p2.tv_sec - p1.tv_sec) * 1000000 + p2.tv_usec - p1.tv_usec;
//     // cout << "paking time in FC is " << (double)deltaT / 1000.0 << "ms" << endl;

//     // ofstream outFile2;
//     // outFile2.open(fc_inputplaintex,ios::ate);  
//     // for(int i=0;i<packednum;i++){
//     //     outFile2<<M[i]<<endl;
//     // }
//     // outFile2.close();

//     //加密M
//     //read t in pk.csv 
// 	ZZX t[k];
//     ifstream inFilepk;
//     inFilepk.open("./file/key/pk-t.csv",ios::in);  
//     if(!inFilepk){
//         cout<<"Cannot open pk-t!"<<endl;
//         exit(0);
//          }  
//         for(int i=0;i<k;i++){
//         inFilepk>>t[i];       
//     }    
//     inFilepk.close();	

	
// 	for(int i=0;i<k;i++){
//         U.decompress(t[i],n,q,dt);
//     }

// 	//read A in pk-A
//     ZZX** A = new ZZX*[k]();
//     for (int i = 0; i < k; i++) {
//         A[i] = new ZZX[k]();
//     }
   
// 	//ZZX A[k][k];
//     ifstream inFile2;
//     inFile2.open("./file/key/pk-A.csv",ios::in);
//      if(!inFile2){
//         cout<<"Cannot open pk-A!"<<endl;
//         exit(0);
//          } 
//     for(int i=0;i<k;i++){
//         for(int j=0;j<k;j++){
//             inFile2>>A[i][j];     
//             //cout<<A[i][j]<<" ";  
//         }
// 	}    
//     inFile2.close(); 
	
//     ZZ_p::init(q);
//     ZZ_pX ppm;
//     SetCoeff(ppm,0,1);
//     SetCoeff(ppm,n,1);
//     ZZ_pXModulus pm(ppm);
//     ZZ_pX p_temp1,p_temp2,p_temp3,p_temp4;
//     //Enc
//     Enc enc;
//     //q/2d
//     ZZ two_d;
//     RR two_d_rr;
//     power2(two_d,dp);//two_d=2^dp
//     two_d_rr=to_RR(two_d);
//     RR q_rr;
//     q_rr = to_RR(q);
//     ZZ q2d;
//     q2d=RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)

//      ZZX MM[packednum];
//     Client_util cu;
//     client_enc_param usp[packednum];
//     for(int i=0;i<packednum;i++)
//         usp[i]=cu.randnum_gen();
	
     
//     //u,vand c, k=2;
// 	ZZX** u1= new ZZX*[packednum]();//input M has 13(outH*outW) cipher,and each cipher is vector and each element is a poly with n degree
//     for (int i = 0; i < packednum; i++) {
//         u1[i] = new ZZX[k]();
//     }
//     ZZX v1[packednum];

  

   
    
//     timeval t1,t2;
//     gettimeofday(&t1,NULL);

//     for(int i=0;i<packednum;i++){        
//         MM[i]=U.reangecoeff(M[i],n);  
// 	    enc.encrypt(MM[i],u1[i],v1[i],n,q,k,q2d,usp[i].r1,A,usp[i].e11,usp[i].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv);
//     }
//     gettimeofday(&t2,NULL);
//     double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
//     cout << "FC层编码和加密时间为:" << deltaT1/1000<<"ms"<< endl;
 
//          //store inputcipher file
//     ofstream outFile1;
//     outFile1.open(fc_inputcipherfile,ios::ate);
//     for(int i=0;i<packednum;i++){
//         for (int j=0;j<k;j++){
//             outFile1<<u1[i][j]<<endl;
//             }     
//         outFile1<<v1[i]<<endl; 
// 	}  
//         outFile1.close();  
       
//         U.deleteZZXArray2(A,k,k);
//         U.deleteZZXArray2(u1,packednum,k);
// }

void Cloud_cipher::input_cipherGen_Noslide_twoPart(int***C_1,int***C_2,int inH,int inW,int kerN,int kerH,int kerW, int packedNum_row, int packedNum_col,const char *cipherfile,
string padding){  
    util U;
    if(padding=="SAME")
    {  
        // int packedNum_row = 2 ;
        // int packedNum_col = 4;
        int packedRowCount = ceil(inH / (packedNum_row * 1.0));//7
        int packedColCount = ceil(inW / (packedNum_col * 1.0));//4
        // cout<<C[0][1][4]<<","<<C[0][1][5]<<","<<C[0][1][6]<<","<<C[0][1][7]<<","
        // <<C[0][2][4]<<","<<C[0][2][5]<<","<<C[0][2][6]<<","<<C[0][2][7]<<endl;
        //将C转化为多项式poly，并在第25,51,77,103,129,155,181,207项装配上C_1的2*4的8个值
        ZZX polys1[kerN][packedRowCount][packedColCount], polys2[kerN][packedRowCount][packedColCount];
        for(int i=0;i<kerN;i++){//6
            for(int j=0;j<packedRowCount;j++){//14
                for(int z=0;z<packedColCount;z++){//7
                    polys1[i][j][z].SetLength(n); 
                     polys2[i][j][z].SetLength(n); 
                    int z1=0;                    
                    for(int z2=0;z2<packedNum_row;z2++){//2                           
                        for(int z3=0;z3<packedNum_col;z3++){  //4                             
                            SetCoeff(polys1[i][j][z],kerH*kerW+(kerH*kerW+1)*z1,C_1[i][j*packedNum_row+z2][z*packedNum_col+z3]);
                            SetCoeff(polys2[i][j][z],kerH*kerW-1+kerH*kerW*z1,C_1[i][j*packedNum_row+z2][z*packedNum_col+z3]);
                            z1++;                               
                        }        
                    }                                   
                }               
            }
        } 
            // cout<<polys[0][1][1]<<endl;            
        enc_param usp[kerN][packedRowCount][packedColCount];
        Cloud_util cdu;
        for(int z=0;z<kerN;z++){
            for(int i=0;i<packedRowCount;i++){
                for(int j=0;j<packedColCount;j++){
                    usp[z][i][j] = U.randnum_gen(); 
                }
            }
        }

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
       
        //u,vand c, k=2;
        ZZX ****u2=U.newZZXArray4(kerN,packedRowCount,packedColCount,k);
        ZZX ***v2=U.newZZXArray3(kerN,packedRowCount,packedColCount);
         ZZX ****u3=U.newZZXArray4(kerN,packedRowCount,packedColCount,k);
        ZZX ***v3=U.newZZXArray3(kerN,packedRowCount,packedColCount);

        timeval p1, p2;
        gettimeofday(&p1, NULL); 
        // 加密
        for(int z=0;z<kerN;z++){//16
            for(int i=0;i<packedRowCount;i++){//7
                for(int j=0;j<packedColCount;j++){ //4              
                    enc.encrypt(polys1[z][i][j],u2[z][i][j],v2[z][i][j],n,q,k,q2d,usp[z][i][j].r1, A, usp[z][i][j].e11,usp[z][i][j].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 
                    enc.encrypt(polys2[z][i][j],u3[z][i][j],v3[z][i][j],n,q,k,q2d,usp[z][i][j].r1, A, usp[z][i][j].e11,usp[z][i][j].e12,p_temp1,p_temp2,p_temp3,p_temp4,pm,t,dt,du,dv); 

                }
            }
        }
        gettimeofday(&p2, NULL);
        double deltaT = (p2.tv_sec - p1.tv_sec) * 1000000 + p2.tv_usec - p1.tv_usec;
        cout << "conv2编码和加密时间是： " << (double)deltaT / 1000.0 << "ms" << endl;
        //  store inputcipher file
        ofstream outFile1;
        outFile1.open(cipherfile,ios::ate);
        for(int z=0;z<kerN;z++){
            for(int i=0;i<packedRowCount;i++){
                for (int j=0;j<packedColCount;j++){
                    for(int u=0;u<k;u++){
                        outFile1<<u2[z][i][j][u]<<endl;
                        outFile1<<u3[z][i][j][u]<<endl;
                    }
                outFile1<<v2[z][i][j]<<endl; 
                outFile1<<v3[z][i][j]<<endl;
                }
            }  
        }
        outFile1.close(); 
        U.deleteZZXArray4(u2,kerN,packedRowCount,packedColCount,k);
        U.deleteZZXArray4(u3,kerN,packedRowCount,packedColCount,k);        
        U.deleteZZXArray2(A,k,k); 
        U.deleteZZXArray3(v2,kerN,packedRowCount,packedColCount);  
        U.deleteZZXArray3(v3,kerN,packedRowCount,packedColCount); 
      }
}

void Cloud_cipher::AbsNoise_MulNoise(int**r, ZZX *d, int inD,int inH,int inW,int kerN,int kerH,int kerW,int strid_H,int strid_W, int pool_stridH, int pool_stridW, string outcipherfile,string cipherfile,helib:: EncryptedArray ea,ZZX *t, ZZX **A, ZZ_pX pm, ZZ_pX p_temp1, ZZ_pX p_temp2, ZZ_pX p_temp3, ZZ_pX p_temp4, ZZ q2d, enc_param ***usp)
{
    int outH = ceil((inH*1.0)/strid_H);
    int outW = ceil((inW*1.0)/strid_W);      
    int datrows = kerH * kerW;
    int datcols = outH * outW; 
    int real_nrow = datrows + 1;            
    int max_packed_num = n/(real_nrow);
    int max_pool_num = max_packed_num / (pool_stridH * pool_stridW);
    int packed_num = max_pool_num * pool_stridH * pool_stridH; 
    int packed_col = ceil(double(1.0*datcols/packed_num));

    util U;
    ZZX u2[kerN][k];
    ZZX v2[kerN];     
    ifstream inFile1;
    inFile1.open(cipherfile, ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   


   for(int z = 0; z < kerN; z++){ 
        for(int j = 0; j < k; j++)                
            inFile1 >> u2[z][j];  
        inFile1 >> v2[z]; 
    }        
    inFile1.close();  
    long slots=ea.size();  
   
   timeval t1,t2;
    gettimeofday(&t1,NULL);
    int p1[kerN][max_pool_num*packed_col];//14*19=266
    //每一个通道编码为一个多项式     
    for(int z1 = 0; z1 < kerN; z1++){ 
        int z2 = 0;              
        for(int i = 0; i < packed_col; i++){//14
            for(int j = 0; j < max_pool_num;j++){//19 
                p1[z1][z2]= r[z1][i];                             
                z2++;
           }          
        }                       
    }       

    int p2[kerN][inH/pool_stridH*inW/pool_stridW];//6*256
    std::vector <vector < long >> p3(kerN, vector<long>(slots,0));

    for(int i = 0; i<kerN; i++){
        for(int j = 0; j < inH/pool_stridH*inW/pool_stridW; j++){
            p2[i][j] = p1[i][j];
            p3[i][j] = p2[i][j];
        }
    }

  
    ZZX encode1_poly[kerN];
    for(int i = 0; i < kerN; i++){
        ea.encode(encode1_poly[i], p3[i]);
    }
    ZZX u3[kerN][k];
    ZZX v3[kerN];
    Enc enc;  
    for(int z=0; z<kerN; z++){  
        enc.encrypt(encode1_poly[z], u3[z], v3[z], n, q, k, q2d, usp[0][0][z].r1, A, usp[0][0][z].e11, 
        usp[0][0][z].e12, p_temp1, p_temp2, p_temp3, p_temp4, pm, t, dt, du, dv); 
    }  
        gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "用于同态加的噪声的编码和加密的时间为:" << deltaT1/1000<<"ms"<< endl;

 timeval t3,t4;
    gettimeofday(&t3,NULL);
    for(int i=0; i<kerN; i++){
        u2[i][0] = u2[i][0] - u3[i][0];
        u2[i][1]= u2[i][1] - u3[i][1];
        v2[i] = v2[i] - v3[i];
        U.Module(u2[i][0], slots, q);
        U.Module(u2[i][1], slots, q);
        U.Module(v2[i], slots, q);
    }
    
   
    for(int i=0; i<kerN; i++){
        for(int j=0; j<k; j++){
            U.decompress(u2[i][j],n,q,du);
        }
        U.decompress(v2[i],n,q,dv);
    }
    //乘上用于求relu的噪声
    ZZX mu0[kerN],mu1[kerN],mv[kerN];
    ZZ_p::init(q2);
     ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4;    
   
    for(int i=0;i<kerN;i++){
        U.Mult(mu0[i], ppm, d[i],u2[i][0],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mu1[i], ppm, d[i],u2[i][1],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mv[i], ppm, d[i], v2[i],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    }
    gettimeofday(&t4,NULL);
    double deltaT2 = (t4.tv_sec - t3.tv_sec) * 1000000 + t4.tv_usec - t3.tv_usec;
    cout << "同态减和同态乘法的计算总时间为:" << deltaT2/1000<<"ms"<< endl;


    ofstream outFile1;
    outFile1.open(outcipherfile, ios::ate); 
    for(int i=0;i<kerN;i++){
        outFile1<<mu0[i]<<endl;
        outFile1<<mu1[i]<<endl;
        outFile1<<mv[i]<<endl;
    }
    outFile1.close();
 
    /*
    //乘上d后解密测试
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,dp);
    two_d_rr=to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);

    ZZX s0u10w1[kerN],s1u11w1[kerN];
    RR two_dp_rr;
    power2(two_dp_rr,dp);
    RR q2_rr;
    ZZ two_d2;
    ZZ two_d3;
    two_d2=two_d/2;
    two_d3=two_d2*(-1);
    q2_rr = q_rr;

    ZZX s[k];
    ifstream inFile22;
    inFile22.open("./file/key/sk.csv", ios::in); 
    if(!inFile22){
        cout<<"Cannot open sk in dec one ploy!"<<endl;
        exit(0);
    }   
    for(int i=0;i<k;i++){
        inFile22>>s[i];            
    } 
    std::vector<vector<long> >pp(kerN,vector<long>(slots,0));
    long p=7681;
    ZZX h5[kerN];
    for(int i=0; i<kerN; i++){
        U.Mult(s0u10w1[i],ppm,s[0],mu0[i],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(s1u11w1[i],ppm,s[1],mu1[i],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        h5[i]=mv[i]-s0u10w1[i]-s1u11w1[i]; 
        for(int j=0;j<n;j++){
            RR coeff_rr;
            ZZ coeff_zz;
            GetCoeff(coeff_zz,h5[i],j);
            coeff_rr=to_RR(coeff_zz);
            coeff_zz=RoundToZZ(coeff_rr*two_dp_rr/q2_rr)%two_d;
            if(coeff_zz>=two_d2){
                coeff_zz= coeff_zz-two_d;
            }
            else if(coeff_zz<=two_d3){
                coeff_zz=coeff_zz+two_d;
            }
            SetCoeff(h5[i],j,coeff_zz);           
        }  
        U.Module(h5[i], n, conv<ZZ>(p)); 
    
        //解码    
        ea.decode(pp[i],h5[i]);      
        for(int j=0; j<slots; j++){
            pp[i][j]=pp[i][j] % p;
            if(pp[i][j]>=p/2){
                pp[i][j]= pp[i][j] - p;
            }
            else if(pp[i][j] <= -p/2){
                pp[i][j]=pp[i][j] + p;
            }
            cout<<pp[i][j]<<" ";
        }
    }
    cout<<"end WX"<<endl;
    cout<<"d[0]_ploy="<<d[0]<<endl;
    cout<<"d[0]="<<endl;
    std::vector<vector<long> >dd(kerN,vector<long>(slots,0));
    ea.decode(dd[0],d[0]);
    for(int i=0; i<slots; i++){
        cout<<dd[0][i]<<" ";
    }    */
}



void Cloud_cipher::AbsNoise_MulNoiseFC(ZZX *d, ZZX **u1, ZZX *v1, int outputlen,string outcipherfile,string cipherfile,helib:: EncryptedArray ea){
     util U;
    ZZX u2[outputlen][k];
    ZZX v2[outputlen];     
    ifstream inFile1;
    inFile1.open(cipherfile,ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   
    for(int z1=0;z1<outputlen;z1++){ 
        for(int j=0;j<k;j++)                
            inFile1>>u2[z1][j];  
        inFile1>>v2[z1];                   
    }        
    inFile1.close();

    long slots=ea.size();
    
   

    for(int i=0;i<outputlen;i++){
        U.decompress(u1[i][0],n,q,du);
        U.decompress(u1[i][0],n,q,du);
        U.decompress(v1[i],n,q,dv);
        u2[i][0]=u2[i][0]+u1[i][0];
        u2[i][1]=u2[i][1]+u1[i][0];
        v2[i]=v2[i]+v1[i];
        U.Module(u2[i][0],n,q2);
        U.Module(u2[i][1],n,q2);
        U.Module(v2[i],n,q2);
    }
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4;  
   
    ZZX mu0[outputlen],mu1[outputlen],mv[outputlen]; 
    ofstream outFile1;
    outFile1.open(outcipherfile,ios::ate); 
    for(int i=0;i<outputlen;i++){
        U.Mult(mu0[i],ppm,d[i],u2[i][0],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mu1[i],ppm,d[i],u2[i][1],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mv[i],ppm,d[i],v2[i],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        outFile1<<mu0[i]<<endl;
        outFile1<<mu1[i]<<endl;
        outFile1<<mv[i]<<endl;
    }
    outFile1.close();
}

void Cloud_cipher::divNoise(ZZX *d, int kerN, string outputcipherfile, string inputcipher, helib:: EncryptedArray ea){
    util U;
    int slots=ea.size();
    ZZX u[kerN][k];
    ZZX v[kerN];     
    ifstream inFile1;
    inFile1.open(inputcipher,ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   
    for(int z1 = 0; z1 < kerN; z1++){ 
        for(int j = 0; j < k; j++){                
            inFile1 >> u[z1][j];            
        }
        inFile1 >> v[z1];               
    }        
    inFile1.close();
    
    timeval t1,t2;
    gettimeofday(&t1,NULL);
    ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4;  
    ZZX mu0[kerN],mu1[kerN],mv[kerN]; 
    for(int i = 0; i < kerN; i++){
        U.decompress(u[i][0],n,q,du);
        U.decompress(u[i][1],n,q,du);
        U.decompress(v[i],n,q,dv); 
        U.Mult(mu0[i],ppm,d[i],u[i][0],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mu1[i],ppm,d[i],u[i][1],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
        U.Mult(mv[i],ppm,d[i],v[i],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    }
    gettimeofday(&t2,NULL);
    double deltaT1 = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    cout << "同态除法的计算总时间为:" << deltaT1/1000<<"ms"<< endl;



    ofstream outFile1;
    outFile1.open(outputcipherfile, ios::ate); 
    for(int i = 0; i < kerN; i++){
        outFile1 << mu0[i] << endl;
        outFile1 << mu1[i] << endl;
        outFile1 << mv[i] << endl;
    }
    outFile1.close();
}

void Cloud_cipher::divNoise(ZZX *d,  string outputcipherfile, string inputcipher, helib:: EncryptedArray ea){
    util U;
    int slots=ea.size();
    ZZX u[k];
    ZZX v;     
    ifstream inFile1;
    inFile1.open(inputcipher,ios::in); 
    if(!inFile1){
        cout<<"Cannot open cipher1 file!\n";
        exit(0);
    }   
   
    for(int j=0;j<k;j++)                
        inFile1>>u[j];  
    inFile1>>v;                
          
    inFile1.close();
     ZZ_p::init(q2);
    ZZ_pX ppm;
    SetCoeff(ppm,0,1);
    SetCoeff(ppm,slots,1);
    ZZ_pX pp_temp1,pp_temp2,pp_temp3,pp_temp4;  
 
    ZZX mu0,mu1,mv; 
    ofstream outFile1;
    outFile1.open(outputcipherfile,ios::ate);    
    U.Mult(mu0,ppm,d[0],u[0],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(mu1,ppm,d[0],u[1],slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    U.Mult(mv,ppm,d[0],v,slots,q2,pp_temp1,pp_temp2,pp_temp3,pp_temp4);
    outFile1<<mu0<<endl;
    outFile1<<mu1<<endl;
    outFile1<<mv<<endl;
    
    outFile1.close();
}
