#include <iostream>
#include <fstream>
#include <math.h>
// #include <NTL/ZZ.h>
// #include <NTL/ZZ_p.h>
// #include <NTL/ZZ_pX.h>
// #include <NTL/ZZX.h>
// #include <NTL/RR.h>
// #include "Binomialdistribution.h"
// #include "util.h"
#include <ctime>


using namespace std;
// using namespace NTL;
int main(){
    // util U;
    // ofstream outFile1;
    // outFile1.open("./file/1.csv", ios::ate); 
    // enc_param* usp3=new enc_param [5];
    // enc_param* usp4=new enc_param [5];
    // enc_param* usp5=new enc_param [5];
    // for(int i=0; i<5; i++){
    //     usp3[i]=U.randnum_gen();
    //     usp4[i]=U.randnum_gen();
    //     usp5[i]=U.randnum_gen();
    //     outFile1 << usp3[i].r1 << endl;
    //     outFile1 << usp3[i].e11 << endl;
    //     outFile1 << usp3[i].e12 << endl;
    // }
       ofstream outFile1;
    outFile1.open("./file/1.csv", ios::ate); 
    int *a = new int [10]();
    for(int i=0; i<10; i++){
        a[i] = rand()%100;
        cout<<a[i]<<" ";
        outFile1 << a[i] << " ";
    }
    outFile1.close();

cout<<endl;
    int *b = new int[5]();
    ifstream inFile1;
    inFile1.open("./file/1.csv",ios::in); 
    if(!inFile1){
        cout<<"Cannot open file1!\n";
        exit(0);
    }   
    for(int i = 0; i < 10; i++){
        inFile1 >> b[i]; 
        cout<<b[i]<<" ";      
    }    
    inFile1.close(); 

return 0;

}