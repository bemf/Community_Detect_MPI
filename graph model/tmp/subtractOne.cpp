#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

using namespace std;

int main(int argc,char* argv[])
{
    ifstream infile("pre_mail");
    ofstream outfile("result",ios::out);
    vector<vector<int> >model;
    model.resize(atoi(argv[1]));
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
    int line=0;
    string input;
    while(!infile.eof())
    {
        infile>>input;
        if(input==";")
        {
            line++;

        }
        else
            model[line].push_back(atoi(input.c_str())-1);
    }

    for(int i=0; i<model.size();i++)
    {
        for(int j=0;j<model[i].size();j++)
        {
            outfile<<model[i][j];
            if(j!=model[i].size()-1)outfile<<" ";
        }
        outfile<<endl;
    }
    return 0;
}
