#include <iostream>
#include <gperftools/profiler.h>
using namespace std;

int main(){
    ProfilerStart("profile/test.prof");
    cout<<"Hello World"<<endl;
    ProfilerEnable();
    return 0;
}