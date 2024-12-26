#include <fstream>
#include <iostream>
#include <stdlib.h> 


using namespace std;

int doubleToInt(double f){
    int i=0;
    if(f>0)
    i=(f*10+5)/10;
    else if(f<0)
    i=(f*10-5)/10;
    else i=0;
    return i;
}

int main(){
    int kerN=6;
    int outH=16;
    int outW=16;
    string weightfile="./file/maxpool.csv";
    // int w_precision=50;
   
    // int w_int[kerN][outN];     
    int data[kerN][outH][outW];

    ifstream file(weightfile, ios::in);
    if(!file){
        cout<<"Cannot open weights!\n";
        exit(0);
    }
    if (file.is_open()){
        for(int i=0;i<kerN;i++){
            for (int j=0; j<outH; j++){
                for(int z=0;z<outW;z++)
                file>>data[i][j][z];
            }
        }
    }  
    file.close();  
  
    // for(int i=0;i<kerN;i++){
    //     for (int j=0; j<outN; j++){                                             
    //         w_int[i][j]= abs(doubleToInt(w_precision * data[i][j]));
    //     }
    // }
    // int max=w_int[0][0];
    int max=data[0][0][0];
    for(int i=0; i<kerN; i++){
        for(int j=0; j<outH; j++){
            for(int z=0;z<outW;z++)
            if(abs(data[i][j][z])>max){
                max=abs(data[i][j][z]);
            }
        }
    }
    cout<<"max="<<max<<endl;
    return 0;
}