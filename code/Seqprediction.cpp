#include "client_cipher.h"
#include <iostream>
#include <keyparaGen.h>
#include <sys/time.h>
#include <cloud_cipher.h>
using namespace std;
using namespace NTL;
util U;
Client_cipher cc;
Client_util cu; 
Cloud_cipher cdc;
Cloud_util cdu;
   
int main()
{   
    timeval off_t1, off_t2;
    gettimeofday(&off_t1,NULL);
    long m = 4096, p=12289, r=1;
    ZZX G;
    G.SetLength(2);
    SetCoeff(G,1,1);
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                            .m(m)
                            .p(p)
                            .r(r)
                            .build();
    // Print the context
    // context.printout(); 
    helib:: EncryptedArray ea(context,G);  
    KeyparaGen keyparaGen; 
    keyparaGen.keygen(); 


    ZZX t[k];
    ifstream inFile1;
    inFile1.open("./file/key/pk-t.csv",ios::in); 
    if(!inFile1){
        cout<<"Cannot open pk-t in conv1!\n";
        exit(0);
    }   
    for(int i=0; i<k; i++){
        inFile1 >> t[i];       
    }    
    inFile1.close();    
    for(int i = 0; i < k; i++){
        U.decompress(t[i], n, q, dt);
    }
  
    ZZX **A = U.newZZXArray2(k, k);
    ifstream inFile2;
    inFile2.open("./file/key/pk-A.csv",ios::in);
    if(!inFile2){
        cout<<"Cannon open pk-A in conv1!\n";
        exit(0);
    }
    for(int i = 0; i < k; i++){
        for(int j = 0; j < k; j++){
            inFile2 >> A[i][j];   
        }
    }    
    inFile2.close();  

 
//     生成加密随机数
//     enc_param *** barusp=U.EncRandGen(16, 16, 16);/
//     store the random number
//     ofstream file1;
//     file1.open("./file/random_num.csv", ios::ate); 
//     for(int i = 0; i<16; i++){
//         for(int j=0; j<16; j++){
//             for(int z1=0; z1<16; z1++){
//                 for(int z2=0; z2<k; z2++){
//                     file1 << barusp[i][j][z1].r1[z2] <<endl;
//                     file1 << barusp[i][j][z1].e11[z2] <<endl;
//                 }
//                 file1 << barusp[i][j][z1].e12 <<endl;
//             }
//         }
//     }
//     file1.close();


    enc_param *** barusp = U.new3(16, 16, 16);
    ifstream file2;
    file2.open("./file/random_num.csv", ios::in);   
    if (!file2)
    {
        cout << "Cannot open random_num in the file!" << endl;
        exit(0);
    }
    for(int i = 0; i<16; i++){
        for(int j=0; j<16; j++){
            for(int z1=0; z1<16; z1++){
                for(int z2=0; z2<k; z2++){
                    file2 >> barusp[i][j][z1].r1[z2];
                    file2 >> barusp[i][j][z1].e11[z2];
                }
                file2 >> barusp[i][j][z1].e12;
            }
        }
    }
    file2.close();

    enc_param *** usp1 = U.new3(6,16,8);
    enc_param *** usp11 = U.new3(1,1,6);
  
    enc_param *** usp2 = U.new3(16,8,4);
    enc_param *** usp12 = U.new3(1,1,16);

    enc_param *** usp13 = U.new3(1,6,14);
    enc_param *** usp14 = U.new3(1,16,4);

     enc_param *** usp01 = U.new3(1,3,14);
    enc_param *** usp02 = U.new3(1,6,4);

    for(int i=0; i<6; i++)
    for(int j=0;j<16; j++)
    for(int z=0; z<8; z++)
     usp1[i][j][z]=barusp[i][j][z];
  
    for(int i = 0; i<6; i++)
        usp11[0][0][i]=barusp[0][0][i];

    for(int i=0; i<16; i++)
    for(int j=0;j<8; j++)
    for(int z=0; z<4; z++)
        usp2[i][j][z]=barusp[i][j][z];

    for(int i = 0; i<16; i++)
        usp12[0][0][i]=barusp[0][0][i];

    for(int i=0; i<6; i++)
    for(int j=0;j<14; j++)
        usp13[0][i][j]=barusp[0][i][j];
    
    for(int i=0; i<16; i++)
    for(int j=0;j<4; j++)
        usp14[0][i][j]=barusp[0][i][j];

    for(int i=0; i<3; i++)
    for(int j=0;j<14; j++)
        usp01[0][i][j]=barusp[0][i][j];

    for(int i=0; i<6; i++)
    for(int j=0;j<4; j++)
        usp02[0][i][j]=barusp[0][i][j];

   

    enc_param *** usp3 = U.new3(1,1,5);
    enc_param *** usp4 = U.new3(1,1,5);
    enc_param *** usp5 = U.new3(1,1,5);
 
    for(int i=0; i<5; i++){
        usp3[0][0][i] = barusp[1][1][i];
        usp4[0][0][i] = usp3[0][0][i];
        usp5[0][0][i] = usp3[0][0][i];
       }





    ZZ two_d;
    RR two_d_rr;
    power2(two_d, dp);//two_d=2^dp
    two_d_rr = to_RR(two_d);
    RR q_rr;
    q_rr = to_RR(q);
    ZZ q2d;
    q2d = RoundToZZ(q_rr/two_d_rr);//q2d=round(q/2^dp)
    ZZ_p::init(q);
    ZZ_pX ppm;
    SetCoeff(ppm, 0, 1);
    SetCoeff(ppm, n, 1);
    ZZ_pXModulus pm(ppm);
    ZZ_pX p_temp1, p_temp2, p_temp3, p_temp4;



    //云处理按照方法1将权重整理为整数矩阵，权重值扩大w_precise倍。不装配偏执 W:kerN × (kerD*kerH*kerW)。 不装配为多项式.
    int** W1 = cdu.plain_weightNobias(6,3,5,5,"./cell/conv1_weight_arr.csv");   
    int** W2= cdu.plain_weightNobias(16,6,5,5,"./cell/conv2_weight_arr.csv");   

    //云处理权重W为多项式矩阵 size:kerN×kerD,按照方法2装配。权重和偏执值扩大w_precise倍三个文件分别是原来的权重，原来的偏执，生成新的权重多项式明文
    cdc.weight_plainGen(6,3,5,5,"./cell/conv1_weight_arr.csv","./cell/conv1_bias_arr.csv","./cell/plain/plain_conv1_weight.csv");
    cdc.weight_plainGen(16,6,5,5,"./cell/conv2_weight_arr.csv","./cell/conv2_bias_arr.csv","./cell/plain/plain_conv2_weight.csv");
 
    // //right
    // cdc.dense_weight_plainGen(120,1024,"./cell/fc1_weight_arr.csv","./cell/fc1_bias_arr.csv","./cell/plain/plain_fc1_weight.csv");
    // cdc.dense_weight_plainGen(84,120,"./cell/fc2_weight_arr.csv","./cell/fc2_bias_arr.csv","./cell/plain/plain_fc2_weight.csv");
    // cdc.dense_weight_plainGen(4,84,"./cell/fc3_weight_arr.csv","./cell/fc3_bias_arr.csv","./cell/plain/plain_fc3_weight.csv");

    //离线阶段客户端生成随机数g
    int ***g1 = U.rand3DGen(3, 32, 32, 10);//10为随机数的范围
    int ***g2 = U.rand3DGen(6, 16, 16, 10);
 
    //客户端按照方法2扩展g并加密得到[g],并存储到"./file/cipher/cipher_conv1_input.csv"
    cc.input_cipherGen_Conv2(g1,3,32,32,6,5,5,1,1,2,2,"./file/cipher/cipher_conv1_input.csv","SAME",t, A,  pm, p_temp1,  p_temp2, p_temp3, p_temp4, q2d, usp01);
    cc.input_cipherGen_Conv2(g2,6,16,16,16,5,5,1,1,2,2,"./file/cipher/cipher_conv2_input.csv","SAME", t, A,  pm, p_temp1,  p_temp2, p_temp3, p_temp4, q2d, usp02);
   
    //云端离线计算w[g]+b,right test on 21.11.25
    cdc.homMultConv_Conv2("./file/homocipher/cipher_conv1_res.csv","./cell/plain/plain_conv1_weight.csv",
    "./file/cipher/cipher_conv1_input.csv",6,3,5,5,32,32,1,1,2,2,"SAME", pm, p_temp1,  p_temp2, p_temp3, p_temp4);
    cdc.homMultConv_Conv2("./file/homocipher/cipher_conv2_res.csv","./cell/plain/plain_conv2_weight.csv",
    "./file/cipher/cipher_conv2_input.csv",16,6,5,5,16,16,1,1,2,2,"SAME",pm, p_temp1,  p_temp2, p_temp3, p_temp4);
 

//生成用于maxpool的噪声
    // cout<<"r="<<endl;
    int **r1 = U.rand2DGen(6, 14, 10);//产生1-10（包括1和10）为随机数的范围 
    ZZX **R1 = U.newZZXArray2(6, 14);
    for(int i = 0; i < 6; i++){
         for(int j = 0; j < 14; j++){           
            R1[i][j].SetLength(n);          
            for(int z1 = 0; z1 < n; z1++){ 
                SetCoeff(R1[i][j], z1, r1[i][j]);    
            }   
         }
    }

    //生成用于maxpool的噪声
    // cout<<"r="<<endl;
    int **r2 = U.rand2DGen(16, 4, 10);//产生1-10（包括1和10）为随机数的范围 
    ZZX **R2 = U.newZZXArray2(16, 4);
    for(int i = 0; i < 16; i++){
         for(int j = 0; j < 4; j++){           
            R2[i][j].SetLength(n);          
            for(int z1 = 0; z1 < n; z1++){ 
                SetCoeff(R2[i][j], z1, r2[i][j]);    
            }   
        }
    }


    //云端离线阶段生成用于扰动relu操作的噪声 
    ZZX encode_d1[6],encode_invd1[6];
    ZZX encode_d2[16],encode_invd2[16];
    // ZZX encode_d3[1], encode_invd3[1];
    // ZZX encode_d4[1], encode_invd4[1];
    // ZZX encode_d5[1], encode_invd5[1];
    
    cdu.GentNoiseforRelu_Conv(encode_d1, encode_invd1, 6, 16, 16, p, 2, ea); 
    cdu.GentNoiseforRelu_Conv(encode_d2, encode_invd2, 16, 8, 8, p, 2, ea);
                        // cdu.GentNoiseforRelu_FC(120, p, 1, encode_d3, encode_invd3, ea);
                        // cdu.GentNoiseforRelu_FC(84, p, 1, encode_d4, encode_invd4, ea);
                        // cdu.GentNoiseforRelu_FC(4, p, 1, encode_d5, encode_invd5, ea);
                        
 
    clock_t startTimeSum,endTimeSum;
    startTimeSum=clock();

    gettimeofday(&off_t2,NULL);
    double deltaT1 = (off_t2.tv_sec - off_t1.tv_sec) * 1000000 + 
    off_t2.tv_usec - off_t1.tv_usec;
    cout << "离线的计算总时间为:" << deltaT1/1000<<"ms"<< endl;

    //在线阶段
    timeval online_t1, online_t2;
    gettimeofday(&online_t1, NULL);
    vector<vector<float>>images;      
    cu.read_data_Images("./file/img_0.csv", images);
    vector<vector<float>>labels;
	cu.read_data_Images("./file/labels.csv", labels);       
    int right_count=0;
    int wrong_count=0;
    int startID=10;
    int endID=11; 

    int ***p1=U.newArray3(3, 32, 32);    
    double predictionID=0;
    int z;
    cout<<"reading image........."<<endl; 
    for(int i1=0; i1<3; i1++) {
        z=0;
        for(int i=0; i<32; i++){
            for(int j=0; j<32; j++){
                p1[i1][i][j]=U.floatToInt(images[i1][z]);
                z++;
            }
        }
    }
    cout<<"The conv1 ....."<<endl;
    cout<<"==========Input size:3×32×32; kern size:6×5×5,strid=1 =========="<<endl;
    //客户端减去扰动x-g发给云端
    timeval m_1, m_2;
    gettimeofday(&m_1, NULL);
    p1 = U.abstrc3D(3, 32, 32, p1, g1); 
    gettimeofday(&m_2, NULL);
    double deltaM1 = (m_2.tv_sec - m_1.tv_sec) * 1000000 + 
    m_2.tv_usec - m_1.tv_usec;
    cout<<"Module A time is "<<deltaM1/1000<<"ms"<<endl;

    timeval m_3, m_4;
    gettimeofday(&m_3, NULL);
    //云端将x-g按照方法1展开//W里不装配偏执b，X处不装备1 
    int** X1= cu.plain_input(p1,3,32,32,5,5,1,1);  //X1:75*1024,W1:6*75    
    //云端明文计算w(x-g)         
    int** WX1=cdu.mat_mult(W1,X1,6,3,5,5,32,32,1,1);   //Wx1:6*1024  
    //云端加密  [w(x-g)+R]    
    cdc.plainEnc_Add(R1,WX1,3,32,32,6,5,5,1,1,2,2,"./file/cipher/cipher_conv1add_input.csv","./file/homocipher/cipher_conv1_res.csv","SAME",
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp13); 
    // 云端同态加 得到w[x]
    // cdc.cipherAdd_Noise(3,32,32,6,5,5,1,1,2,2,"./file/cipher/cipher_conv1add_input.csv",
    //  "./file/homocipher/cipher_conv1_res.csv","./file/cipher/cipher_aconv1_input.csv");
   gettimeofday(&m_4, NULL);
    double deltaM2 = (m_4.tv_sec - m_3.tv_sec) * 1000000 + 
    m_4.tv_usec - m_3.tv_usec;
    cout<<"Module B time is "<<deltaM2/1000<<"ms"<<endl;

    timeval m_5, m_6;
    gettimeofday(&m_5, NULL);
    //客户端解密并求maxpool
    // cout<<"result_maxpool1="<<endl;
    int ***result_maxpool1 = cc.DecCNN2D_Conv1(3,32,32,6,5,5,1,1,2,2,"./file/cipher/cipher_conv1add_input.csv","./file/key/sk.csv",
    "SAME");
    //客户端对maxpool的结果用BGV编码后,再加密成一条密文            
    cu.Encode_Enc(result_maxpool1, 6, 16, 16, "./file/cipher/cipher_maxpool_res1.csv", ea,
    t, A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp11); 
    gettimeofday(&m_6, NULL);
    double deltaM3 = (m_6.tv_sec - m_5.tv_sec) * 1000000 + 
    m_6.tv_usec - m_5.tv_usec;
    cout<<"Module C time is "<<deltaM3/1000<<"ms"<<endl; 


    timeval m_7, m_8;
    gettimeofday(&m_7, NULL);
    //云端对maxpool加密后的结果减去扰动，再乘上一个新扰动  
    cdc.AbsNoise_MulNoise(r1,encode_d1, 3,32,32,6,5,5,1,1,2,2,"./file/cipher/absnoise_mulnoise_res1.csv",
    "./file/cipher/cipher_maxpool_res1.csv",  ea,
     t, A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp11);
      gettimeofday(&m_8, NULL);
    double deltaM4 = (m_8.tv_sec - m_7.tv_sec) * 1000000 + 
    m_8.tv_usec - m_7.tv_usec;
    cout<<"Module D time is "<<deltaM4/1000<<"ms"<<endl; 


    timeval m_9, m_10;
    gettimeofday(&m_9, NULL);
    // 客户端解密
     ZZX h1[6];
    cu.Genl_Dec(h1,6,"./file/cipher/absnoise_mulnoise_res1.csv","./file/key/sk.csv",n,p);
    // 客户端解码，带着正倍数求relu.  
    int*** relu_result1=cu.Decode_Relu(h1, 6, 16, 16, ea, p);
    //客户端加密带倍数的结果
    cu.Encode_Enc(relu_result1, 6, 16, 16, "./file/cipher/cipher_times_relu1.csv", ea,
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp11); 
    gettimeofday(&m_10, NULL);
    double deltaM5 = (m_10.tv_sec - m_9.tv_sec) * 1000000 + 
    m_10.tv_usec - m_9.tv_usec;
    cout<<"Module E time is "<<deltaM5/1000<<"ms"<<endl; 

    timeval m_11, m_12;
    gettimeofday(&m_11, NULL);
    //云同态去掉倍数
    cdc.divNoise(encode_invd1, 6, "./file/cipher/cipher_relu_res1.csv", 
    "./file/cipher/cipher_times_relu1.csv", ea);
    gettimeofday(&m_12, NULL);
    double deltaM6 = (m_12.tv_sec - m_11.tv_sec) * 1000000 + 
    m_12.tv_usec - m_11.tv_usec;
    cout<<"Module F time is "<<deltaM6/1000<<"ms"<<endl; 
    
    timeval m_13, m_14;
    gettimeofday(&m_13, NULL);
    //客户端解密    
    ZZX y1[6];
    cu.Genl_Dec(y1, 6, "./file/cipher/cipher_relu_res1.csv","./file/key/sk.csv", n, p);
    int *** f1=cu.Decode(y1, 6, 16, 16, ea, p);
      gettimeofday(&m_14, NULL);
    double deltaM7 = (m_14.tv_sec - m_13.tv_sec) * 1000000 + 
    m_14.tv_usec - m_13.tv_usec;
    cout<<"Module G time is "<<deltaM7/1000<<"ms"<<endl; 
 
    cout<<"===================================="<<endl;     
    cout<<"The conv2 ...  input size:6×16×16;kern size:16×6×5×5,strid=1"<<endl;
    f1=U.abstrc3D(6,16,16,f1,g2); 
    int** X2 = cu.plain_input(f1,6,16,16,5,5,1,1);
    int** WX2 = cdu.mat_mult(W2,X2,16,6,5,5,16,16,1,1);    
    cdc.plainEnc_Add(R2,WX2,6,16,16,16,5,5,1,1,2,2,"./file/cipher/cipher_conv2add_input.csv","./file/homocipher/cipher_conv2_res.csv","SAME",
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp14);
    int ***result_maxpool2=cc.DecCNN2D_Conv1(6,16,16,16,5,5,1,1,2,2,"./file/cipher/cipher_conv2add_input.csv","./file/key/sk.csv",
    "SAME");

   
    cu.Encode_Enc(result_maxpool2, 16, 8, 8, "./file/cipher/cipher_maxpool_res2.csv", ea,
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp12); 
    cdc.AbsNoise_MulNoise(r2,encode_d2, 6,16,16,16,5,5,1,1,2,2,"./file/cipher/absnoise_mulnoise_res2.csv",
    "./file/cipher/cipher_maxpool_res2.csv", ea,
    t, A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp12);
    ZZX h2[16];
    cu.Genl_Dec(h2, 16, "./file/cipher/absnoise_mulnoise_res2.csv", "./file/key/sk.csv", n, p);
    int*** relu_result2=cu.Decode_Relu(h2, 16, 8, 8, ea, p);
    //客户端加密带倍数的结果
    cu.Encode_Enc( relu_result2, 16, 8, 8, "./file/cipher/cipher_relu_res2.csv", ea,
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp12); 
    cdc.divNoise(encode_invd2, 16, "./file/cipher/cipher_relu_res2.csv", "./file/cipher/cipher_relu_res2.csv", ea);
    ZZX y2[16];
    cu.Genl_Dec(y2, 16, "./file/cipher/cipher_relu_res2.csv", "./file/key/sk.csv", n, p);
    int *** f2=cu.Decode(y2, 16, 8, 8, ea, p);
    int *F1=cu.flatten(f2, 16, 8, 8);   

    cout<<"================================================="<<endl;
    cout<<"The dense 1.."<<endl;
    cout<<"==========kern size: 120×1024; Input size: (16*8*8)×1=1024×1=========="<<endl;
    cc.dense_input_cipherGen(F1, 1024, "./file/cipher/cipher_fc1_input.csv",
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp3[0][0]);
    cdc.homMultHaveAdd("./file/homocipher/cipher_fc1_res.csv","./cell/plain/plain_fc1_weight.csv","./file/cipher/cipher_fc1_input.csv",120,1024);       
    int *relu_fc1=cc.DecDenseHaveAdd_relu(120, 1024, "./file/homocipher/cipher_fc1_res.csv");

    cout<<endl<<endl;
    cout<<"The dense 2.."<<endl;
    cout<<"==========kern size: 84×120; Input size: 120×1=========="<<endl;
    cc.dense_input_cipherGen(relu_fc1,120,"./file/cipher/cipher_fc2_input.csv",
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp4[0][0]);
    cdc.homMultHaveAdd("./file/homocipher/cipher_fc2_res.csv","./cell/plain/plain_fc2_weight.csv","./file/cipher/cipher_fc2_input.csv",84,120); 
    int *relu_fc2=cc.DecDenseHaveAdd_relu(84,120,"./file/homocipher/cipher_fc2_res.csv");
    
    cout<<endl<<endl;
    cout<<"The dense 3.."<<endl;
    cout<<"==========kern size: 4×84;Input size: 84×1=========="<<endl;
    cc.dense_input_cipherGen(relu_fc2, 84, "./file/cipher/cipher_fc3_input.csv",
    t,A, pm, p_temp1, p_temp2, p_temp3, p_temp4, q2d, usp5[0][0]);
    cdc.homMultHaveAdd("./file/homocipher/cipher_fc3_res.csv","./cell/plain/plain_fc3_weight.csv","./file/cipher/cipher_fc3_input.csv",4,84);
    predictionID=cc.DecDenseHaveAdd_softmax(4, 84, "./file/homocipher/cipher_fc3_res.csv");
    gettimeofday(&online_t2, NULL);
    double deltaT2 = (online_t2.tv_sec - online_t1.tv_sec) * 1000000 + 
    online_t2.tv_usec - online_t1.tv_usec;
    cout<<"To predict One image time is "<<deltaT2/1000<<"ms"<<endl;
    cout<<"label value is :"<<labels[0][0]<<endl;
    cout<<"***********************************"<<endl;
    cout<<"***********************************"<<endl<<endl;
    
    if(predictionID==labels[0][0])
        right_count++;
    else {
        wrong_count++;
        cout<<"预测的编号为： "<<predictionID <<endl;
        }   

    //释放内存        
    U.deleteArray2(W1, 6, 75); 
    U.deleteArray2(W2, 16, 150);
    U.deleteArray2(X1, 75, 1024);
    U.deleteArray2(X2, 150, 256);
    U.deleteArray2(WX1, 6, 1024);
    U.deleteArray2(WX2, 16, 256);
    U.deleteArray3(g1, 3, 32, 32);
    U.deleteArray3(g2, 6, 16, 16);
    U.deleteArray3(p1, 3, 32, 32);
    U.deleteArray3(result_maxpool1, 6, 16, 16);
    U.deleteArray3(relu_result1, 6, 16, 16);
    U.deleteArray3(result_maxpool2, 16, 8, 8);
    U.deleteArray3(relu_result2, 16, 8, 8);
    U.deleteArray3(f1, 6, 16, 16);
    U.deleteArray3(f2, 16, 8, 8);  
    
   
    delete[] F1;
    delete[] relu_fc1;
    delete[] relu_fc2;   
    
     
    cout<<"The right count = "<<right_count<<endl;
    cout<<"The wrong count = "<<wrong_count<<endl;
    cout<<"test result================================="<<endl;
    double right_rate=double(right_count)/(endID-startID);
    double wrong_rate=double(wrong_count)/(endID-startID);
    cout<<"right_rate = "<< right_rate << endl;
    cout<<"wrong_rate = " << wrong_rate << endl;  

    endTimeSum=clock();
    cout<<"To predict "<< endID-startID << " images total time is "<<(double)(endTimeSum-startTimeSum)/CLOCKS_PER_SEC<<"s"<<endl;

return 0;
}