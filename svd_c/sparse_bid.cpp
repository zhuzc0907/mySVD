#include "mySparse.h"

int main(int argc, char** argv){

    SpMat mat = load_mtx("data/worms20_10NN/worms20_10NN.mtx");
    // SpMat mat = load_mtx("data/test/test.mtx");

    cout<<mat.nonZeros()<<endl;

    golub_kahan(mat);
    return 0;
}