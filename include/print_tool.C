#include "print_tool.h"

void print_vector_info(string message, vector<int> vec)
{
    for(std::size_t i=0; i!=vec.size(); ++i){
        if(i==0) printf("%s = %d ", message.c_str(), vec[i]);
        else     printf(" %d", vec[i]);
    }
}

void print_vector_info(string message, vector<double> vec)
{
    for(std::size_t i=0; i!=vec.size(); ++i){
        if(i==0) printf("%s = %7.4f ", message.c_str(), vec[i]);
        else     printf(" %7.4f", vec[i]);
    }
}

