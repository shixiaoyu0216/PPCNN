#include "util.h"
#include <stdlib.h>


util::util()
{

}

enc_param util::randnum_gen(){
    enc_param usp;
	for(int i=0;i<k;i++){
        usp.r1[i].SetLength(n);
        usp.Bin_dis.Binomial_distribution_gen(usp.r1[i],n,ita);
        }
    for(int i=0;i<k;i++){
       usp.e11[i].SetLength(n);
       usp.Bin_dis.Binomial_distribution_gen(usp.e11[i],n,ita);
    }
    usp.Bin_dis.Binomial_distribution_gen(usp.e12,n,ita);
    return usp;
}


enc_param ***util::EncRandGen(int cha, int row,int col){
    enc_param ***usp;
    usp = new enc_param**[cha];
        for (int i = 0; i < cha; i++) {
            usp[i] = new enc_param*[row];
            for (int j = 0; j < row; j++)
                usp[i][j] = new enc_param[col];
            }
       util U;
        for(int z=0; z<cha; z++){
            for(int i=0; i<row; i++){
                for(int j=0; j<col; j++){
                    usp[z][i][j] = U.randnum_gen(); 
                }
            }
        }
    return usp;
}


int ***util::abstrc3D(int cha,int row,int col, int ***p, int ***g){
    for(int i=0; i<cha; i++){
        for(int j=0; j<row; j++){
            for(int u=0; u<col; u++){
                p[i][j][u] -= g[i][j][u];
            }
        }
    }
    return p;
}

int *util::abstrc1D(int col, int *p, int *g){

    for(int u = 0; u < col; u++){
        p[u] -= g[u];
    }
    return p;
}

int ** util::rand2DGen(int maxcha,int maxrow, int r){
    util U;
    int **g=U.newArray2(maxcha, maxrow);
    srand(time(NULL));
    for(int z=0;z<maxcha;z++){
        for(int i=0;i<maxrow;i++){            
            g[z][i]=std::rand()%(r+1);  
            // g[z][i]=0;         
            }
        }
   
    return g;
}

int *** util::rand3DGen(int maxcha,int maxrow, int maxcol, int r){
    util U;
    int ***g = U.newArray3(maxcha, maxrow, maxcol);
    srand(time(NULL));
    for(int z = 0; z < maxcha; z++){
        for(int i = 0; i < maxrow; i++){
            for(int j = 0; j < maxcol; j++){            
                g[z][i][j] = std::rand()%(r+1);        
            }
        }
    }
    return g;
}




int util::gcdExtended(int a, int b,int *x, int *y){
    if(a==0){
        *x=0;
        *y=1;
        return b;
    }
    int x1,y1;
    int gcd=util::gcdExtended(b%a, a, &x1, &y1);
    *x=y1-(b/a)*x1;
    *y=x1;
    return gcd;
}








int *util::rand1DGen(int lenth, int r){
    util U;
    int *g;
    g=new int [lenth];
    srand(time(NULL));
    for(int z=0;z<lenth;z++){
        g[z]=std::rand()%r+1;           
    }
    return g;
}




int util::ReverseInt(int i){
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 255;
	ch2 = (i >> 8) & 255;
	ch3 = (i >> 16) & 255;
	ch4 = (i >> 24) & 255;
	return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}
void util::rand(ZZX& r, int n, ZZ q){
    for (int i = 0; i < n; i++){
        SetCoeff(r,i,-q / 2+RandomBnd(q));
    }
}

int util::doubleToInt(double f){
    int i=0;
    if(f>0)
    i=(f*10+5)/10;
    else if(f<0)
    i=(f*10-5)/10;
    else i=0;
    return i;
}
int util::floatToInt(float f){
    int i=0;
    if(f>0)
    i=(f*10+5)/10;
    else if(f<0)
    i=(f*10-5)/10;
    else i=0;
    return i;
}

void util::Module(ZZX &modx, int n, ZZ q){
    ZZ x;
    for(int i = 0;i<n;i++){
        GetCoeff(x,modx,i);
        x=x%q;
        if ((x <= (q-1) / 2) && (x >= -(q-1) / 2)) {
        } else if (x > 0) {
            x = x - q;
        } else if (x < 0) {
            x = x + q;
        }
        SetCoeff(modx,i,x);
    }
}

ZZX util::reangecoeff(ZZX m2,int n){
    ZZ tempm2; 
    ZZX m3; 
    m3.SetLength(n);
    int j=1;
    GetCoeff(tempm2,m2,0);
    SetCoeff(m3,0,tempm2);
    for(int i=n-1;i>0;i--){
        GetCoeff(tempm2,m2,i);
        SetCoeff(m3,j,-1*tempm2);
        j++;
    }
    return m3;
 } 

//   void util::Mult(ZZX &pr, ZZ_pX F,ZZX& p1, ZZX& p2,int n,ZZ q,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4){
//     p_temp1=to_ZZ_pX(p1);
//     p_temp2=to_ZZ_pX(p2);    
//     MulMod(p_temp3,p_temp1,p_temp2,F);
//     pr=to_ZZX(p_temp3);
//     Module(pr,n,q); 
//    }

void util::Mult(ZZX &pr, ZZ_pX pm,ZZX& p1, ZZX& p2,int n,ZZ q,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4){
     p_temp1=to_ZZ_pX(p1);
     p_temp2=to_ZZ_pX(p2);
     p_temp3=p_temp1*p_temp2;
    // p_temp3 = to_ZZ_pX(p1) * to_ZZ_pX(p2);
    // cout<<"p_tmp3 finish==============="<<endl;
    rem(p_temp4,p_temp3,pm);//p_temp4=p_temp3%pm
    pr=to_ZZX(p_temp4);
    // cout<<"I am out to_ZZX"<<endl;
    // this->Module(pr,n,q);
    // cout<<"module out............."<<endl;
  
}

void util::Mult_para(ZZX &pr, ZZ_pX pm,ZZX& p1, ZZX& p2,int n,ZZ q,ZZ_pX p_temp1,ZZ_pX p_temp2,ZZ_pX p_temp3,ZZ_pX p_temp4){
    
    //  p_temp1=to_ZZ_pX(p1);
    //  p_temp2=to_ZZ_pX(p2);
    //  p_temp3=p_temp1*p_temp2;
    p_temp3 = to_ZZ_pX(p1) * to_ZZ_pX(p2);
    cout<<"p_tmp3 finish==============="<<endl;
    rem(p_temp4,p_temp3,pm);//p_temp4=p_temp3%pm
    pr=to_ZZX(p_temp4);
    cout<<"I am out to_ZZX"<<endl;
    this->Module(pr,n,q);
    cout<<"module out............."<<endl;
  
}

void util::compress(ZZX &pc, int n, ZZ q, int d){
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,d);
    two_d_rr=to_RR(two_d);
    ZZ temp_compress;
    RR temp_compress_rr;
    RR q_rr;
    q_rr = to_RR(q);
    ZZ two_d2;
    ZZ two_d3;
    two_d2=two_d/2;
    two_d3=two_d2*(-1);
    for (int i = 0; i < n; i++) {
        GetCoeff(temp_compress,pc,i);
        temp_compress_rr=to_RR(temp_compress);
        temp_compress=RoundToZZ(two_d_rr/q_rr*temp_compress_rr);
        temp_compress=temp_compress%two_d;      

         if(temp_compress>=two_d2){
             temp_compress= temp_compress-two_d;
        }else if(temp_compress<=two_d3){
            temp_compress=temp_compress+two_d;
        }

          SetCoeff(pc,i,temp_compress);
    }
    //Module(pc,n,q);
}

void util::decompress(ZZX &pdc, int n, ZZ q, int d){
    ZZ two_d;
    RR two_d_rr;
    power2(two_d,d);
    two_d_rr=to_RR(two_d);
    ZZ temp_decompress;
    RR temp_decompress_rr;
    RR q_rr;
    q_rr = to_RR(q);
    for(int i=0;i<n;i++){
        GetCoeff(temp_decompress,pdc,i);
        temp_decompress_rr=to_RR(temp_decompress);
        temp_decompress=RoundToZZ(q_rr/two_d_rr*temp_decompress_rr);
        SetCoeff(pdc,i,temp_decompress);
    }
}


