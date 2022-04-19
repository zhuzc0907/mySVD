#ifndef __MYSPARSE_H__
#define __MYSPARSE_H__

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <chrono>
#include <gperftools/profiler.h>

using namespace Eigen;
using std::cout;
using std::endl;
using std::string;

typedef SparseMatrix<double> SpMat; //declares a column-major sparse matrix type of double
typedef Triplet<double> T;
typedef Matrix<double, Dynamic, 1> VectorXd;


SpMat load_mtx(string file_name){
    std::ifstream fin(file_name);
    // M: number of rows
    // N: number of columns
    // L: number of lines
    int M, N, L;
    
    // Ignore headers and comments
    while (fin.peek() == '%'){
        fin.ignore(2048, '\n');
    }
    fin >> M >> N >> L;

    std::vector<T> tripletList;
    tripletList.reserve(L); // Input L for estimation of entries
    // Fill the triplet list
    for(int i=0; i<L; i++){
        int x, y;
        double value;
        fin >> x >> y >> value;

        // Check if the original data matrix is 0-index-based or not!!
        tripletList.push_back(T(x-1, y-1, value));
    }
    fin.close();
    
    // Convert the triplet list to sparse matrix
    SpMat mat(M, N);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

VectorXd build_equal_vector(int size, double value){
    VectorXd v(size);
    for(int i=0; i<size; i++){
        v(i) = value;
    }
    return v;
}

void save_vector_to_csv(std::vector<double> data, string filename){
    std::ofstream rtn(filename);
    for(int i=0; i<data.size(); i++){
        rtn << data[i] <<endl;
    }
}

SpMat build_bid_matrix(std::vector<double> alpha_list, std::vector<double> beta_list){
    int a_size = alpha_list.size();
    int b_size = beta_list.size();

    std::vector<T> tripletList;
    tripletList.reserve(a_size+b_size); 
    // Fill the triplet list
    for(int i=0; i<a_size; i++){
        tripletList.push_back(T(i, i, alpha_list[i]));
    }
    for(int i=0; i<b_size; i++){
        tripletList.push_back(T(i, i+1, beta_list[i]));
    }
    SpMat mat(a_size, a_size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

// suppopse to be very slow..
double bid_error(SpMat B, std::vector<VectorXd> U_temp, std::vector<VectorXd> V_temp, SpMat A){

}

double error_delt(double alpha_2, double beta_1, VectorXd v_2, VectorXd u_2, VectorXd u_1, SpMat A, double last_delt){
    double delt = A*v_2 - u_1*beta_1 - u_2*alpha_2;
    return (last_delt + delt);

}

void golub_kahan(SpMat mat){
    const int M = mat.rows();
    const int N = mat.cols();


    VectorXd v = build_equal_vector(N, 1/sqrt(N));
    VectorXd u = build_equal_vector(M, 0);
    double alpha = 0;
    double beta = 0;

    // U, V ...
    std::vector<VectorXd> V_temp;
    std::vector<VectorXd> U_temp;

    std::vector<double> alpha_list;
    std::vector<double> beta_list;

    std::vector<double> u_orth_error;
    std::vector<double> v_orth_error;

    double error_delt = 0.0;
    std::vector<double> delt_list;
    beta_list.push_back(0.0);
    // std::vector<double> iteration_time;
    // ProfilerStart("profile/profile_loop.prof");
    for(int i=0; i<N; i++){
        // std::chrono::steady_clock::time_point end;
        // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        V_temp.push_back(v);        
        u = mat*v - beta*u;
        alpha = u.norm();
        u = u/alpha;
        U_temp.push_back(u);
        v = mat.transpose()*u - alpha*v;
        beta = v.norm();
        v = v/beta;

        alpha_list.push_back(alpha);
        beta_list.push_back(beta);

        // cout<<"Iteration "<<i<<": "<<alpha<<" "<<beta<<endl;
        // cout<<U_temp[0].transpose()*u<<" "<<V_temp[0].transpose()*v<<endl;
        // u_orth_error.push_back(U_temp[0].transpose()*u);
        // v_orth_error.push_back(V_temp[0].transpose()*v);

        // WHEN TO STOP?

        // iteration_time.push_back(std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count());

        error_delt = error_delt(alpha, beta_list.end()[-2], v, u, U_temp.end()[-2], A, error_delt);
        delt_list.push_back(error_delt);
    }
    // ProfilerStop();
    
   
    cout<<"OK"<<endl;
    // cout<<v(1)<<endl;

    // save_vector_to_csv(u_orth_error, "u_orth_error.csv");
    // save_vector_to_csv(v_orth_error, "v_orth_error.csv");
    // save_vector_to_csv(iteration_time, "iteration_time.csv");
    save_vector_to_csv(error_delt, "results/error_delt.csv")

    
}



#endif

