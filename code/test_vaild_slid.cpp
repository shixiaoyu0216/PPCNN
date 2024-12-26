
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;
int main(){ 
    // int inD=1;
    // int inH=4, inW=4, kerH=2,kerW=2,strid_H=1,strid_W=1;
  
    // int A[inD][inH][inW];
    // int u=0;
    // for(int i=0;i<inD;i++){
    //     for(int j=0;j<inH;j++){
    //         for(int z=0;z<inW;z++){                
    //             A[i][j][z]=u+1;
    //             u++;
    //             cout<<A[i][j][z]<<" ";
                
    //             }
    //             cout<<endl;
    //     }
    // }
    // cout<<endl;

    // int outH=ceil(float(inH-kerH+1)/float(strid_H));
    // int outW=ceil(float(inW-kerW+1)/float(strid_W));
    //   int C[inD][kerH*kerW][outH*outW];
   
    
    // //after padding,slide window on A and transform to C
    // for(int z=0;z<inD;z++){
    //     int y=0;
    //     for(int i=0;i<outH;i++){
    //         for(int j=0;j<outW;j++){
    //             int x=0;
    //             for(int u=0;u<kerH;u++)
    //                 for(int v=0;v<kerW;v++){
    //                     C[z][x][y]=A[z][u+i][v+j];
    //                     cout<<C[z][x][y]<<" ";
    //                     x=x+1;
    //                     }
    //                 y=y+1;
    //                 cout<<endl;
    //             }
    //         }    
    //     }
    //     cout<<endl;

   int kerN=1, outH=4, outW=4,pool_sizeH=3,pool_sizeW=3,
   pool_stridH=2, pool_stridW=2;

    //padding when padding=SAME
    int in2_H=ceil(outH/pool_stridH);
    int in2_W=ceil(outW/pool_stridW);
    int pad2_H=max((in2_H-1)*pool_stridH+pool_sizeH-outH,0);
    int pad2_W=max((in2_W-1)*pool_stridW+pool_sizeW-outW,0);
   
    int pad2_top=int((pad2_H)/2);
    int pad2_left=int((pad2_W)/2);
    
    int input[kerN][outH][outW];
    int u=0;
    cout<<"*********input********"<<endl;
    for(int i=0;i<kerN;i++){
        for(int j=0;j<outH;j++){
            for(int z=0;z<outW;z++){                
                input[i][j][z]=u+1;
                u++;
                cout<<input[i][j][z]<<" ";
                
                }
                cout<<endl;
        }
    }
    cout<<endl;
    cout<<"*******************"<<endl;
  
    //B is padding 0 in P
    int B[kerN][outH+pad2_H][outW+pad2_W];
    for(int z=0;z<kerN;z++){
        for (int i=0;i<outH+pad2_H;i++){
            for (int j=0;j<outW+pad2_W;j++){
                B[z][i][j]=0;
            }
        }
    }  
    for(int z=0;z<kerN;z++){     
        int i=0;
        for (int u=0;u<outH;u++) { 
            int j=0;
            for(int v=0;v<outW;v++){ 
                B[z][i+pad2_top][j+pad2_left]=input[z][u][v];
                j=j+1;
            }
            i=i+1;
        }
    }    
         
    cout<<"Test padding..."<<endl;
    cout<<"B[u][i][j]="<<endl;
    for(int u=0;u<kerN;u++){
        for (int i=0;i<outH+pad2_H;i++){
            for (int j=0;j<outW+pad2_W;j++){
                cout<<B[u][i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl; 
    }  
    cout<<"outH+pad2_H="<<outH+pad2_H<<endl;
      cout<<"outW+pad2_W="<<outW+pad2_W<<endl;
    //maxpool in every fearture map    
    int P_res[kerN][in2_H][in2_W];    
    int max;
    for(int z=0;z<kerN;z++){
        int u1=0,u2=0;
        for(int k1=0;k1<=outH+pad2_H-pool_sizeH;k1=k1+2){
            for(int k2=0;k2<=outW+pad2_W-pool_sizeW;k2=k2+2){
                max = B[z][k1][k2];
                cout<<B[z][k1][k2]<<" "<<endl;
                for(int i=k1;i<k1+pool_sizeH;i++){
                    for(int j=k2;j<k2+pool_sizeW;j++){
                        cout<<B[z][i][j]<<" ";
                        if(B[z][i][j]>max){
                            max = B[z][i][j];
                        }
                    }
                    cout<<endl;
                }
                cout<<endl;
                P_res[z][u1][u2]=max;
                u2++;
                if(u2%in2_W==0){
                    u1++;
                    u2=0;
                }

            }
            
        }        
    }
       
    cout<<"test maxpool.."<<endl;
    cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    cout<<"P_res[0][i][j]="<<endl;
    for (int i=0;i<in2_H;i++){
        for (int j=0;j<in2_W;j++){
            cout<<P_res[0][i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
   cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
}
