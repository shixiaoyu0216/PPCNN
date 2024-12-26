#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
vector<vector<float> >a;  //二维数组存储读入变量
vector<float>b;
int c[6][32][32];

inline void file_to_string(vector<string> &record, const string& line, char delimiter);
inline float string_to_float(string str);

inline void file_to_string(vector<string> &record, const string& line, char delimiter)
{
    int linepos=0;
    char c;
    int linemax=line.length();
    string curstring;
    record.clear();
    while(linepos<linemax)
    {
        c = line[linepos];
        // cout<<c;       //输出文件每行内容，分析后进行提取
        if(isdigit(c)||c=='.'){
            curstring+=c;
        }
        else if(c==delimiter&&curstring.size()){
            record.push_back(curstring);
            curstring="";
        }
        ++linepos;
    }
    cout<<'\n';
    if(curstring.size())
        record.push_back(curstring);
    return;
}
 
inline float string_to_float(string str){
    int i=0,len=str.length();
    float sum=0;
    while(i<len){
        if(str[i]=='.') break;
        sum=sum*10+str[i]-'0';
        ++i;
    }
    ++i;
    float t=1,d=1;
    while(i<len){
        d*=0.1;
        t=str[i]-'0';
        sum+=t*d;
        ++i;
    }
    return sum;
}

int floatToInt(float f){
    int i=0;
    if(f>0)
    i=(f*10+5)/10;
    else if(f<0)
    i=(f*10-5)/10;
    else i=0;
    return i;
}


int main()
{  
    vector<string> row;
    string line;
    string filename;
    ifstream in("./file/f0.csv");  //读取工程文件下的.csv文件
    if (in.fail())  
    { cout << "File not found" <<endl; return 0; } //读取失败，返回为假
    while(getline(in, line)  && in.good() )
    {
        file_to_string(row, line, ',');  //把line里的单元格数字字符提取出来，“,”为单元格分隔符
        for(int i=0, leng=row.size(); i<leng; i++){
            // cout << row[i] << " ";
            b.push_back(string_to_float(row[i]));
        }
        a.push_back(b);
        b.clear();
    }
    in.close();
    // cout<<a.size()<<endl;
    // cout<<a[0].size()<<endl;
    
    ofstream out("conv_original.csv");  //写入到文件conv
    for(int i=0; i<a.size(); ++i){
        int u=0;
        for(int j=0; j<16; ++j){
            for(int z=0; z<16; ++z){
                c[i][j][z]=a[i][u];
                out<<floatToInt(c[i][j][z]*10)<<" ";
                u++;
            }
            out<<'\n';
        }
        out<<'\n';
    }
    out<<'\n';
    out<<'\n';
    out.close();
    return 0;  
}