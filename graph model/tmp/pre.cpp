#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

int main(int argc,char* argv[])
{
    ifstream infile("email");
    ofstream outfile("pre_mail",ios::out);
    //int times=0;
    //int input,index=0;
    //while(!infile.eof())
    //{
        //infile>>input;
        //if(times%2==0)index=input-1;
        //else
        //{
            //model[index].push_back(input-1);
            //model[input-1].push_back(index);
        //}
        //times++;
    //}

    //for(int i=0; i<model.size();i++)
    //{
        //sort(model[i].begin(),model[i].end());
        //outfile<<i;
        //for(int j=0;j<model[i].size();j++)
        //{
            //outfile<<" "<<model[i][j];
        //}
        //outfile<<endl;
    //}
    char c;
    while(!infile.eof())
    {
        infile.get(c);
        if(c=='\n')outfile<<" ;";
        outfile<<c;
    }

    return 0;
}
