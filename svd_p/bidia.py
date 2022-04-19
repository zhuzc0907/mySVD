import numpy as np
import scipy.sparse as scsp

def GKL_bidiagonalization(A):
    A = np.matrix(A)
    m,n = A.shape
    
    p = max(m, n)
    alpha = np.zeros(p)
    beta = np.zeros(p-1)
    U = np.zeros((p, m))  # in fact it's U.H
    V = np.zeros((p, n))  # in fact it's V.H
    # Transposed matrices are used for easier slicing
    U[0, 0] = 1
    
    for i in range(p):
        V[i] = A.H @ U[i] - beta[i-1]*V[i-1]
        alpha[i] = np.linalg.norm(V[i])
        V[i] /= alpha[i]

        if i > p - 2: continue
        U[i+1] = A @ V[i] - alpha[i]*U[i]
        beta[i] = np.linalg.norm(U[i+1])
        U[i+1] /= beta[i]
    U,V = map(np.matrix, (U, V))
        
    return U.H, (alpha, beta), V.H 

if __name__ == '__main__':
    # A = [[1, 2], [2, 3], [3, 5]]
    # A = np.array([1,2,3,4,5,6,7,8,9,10,11,12]).reshape(6,-1)
    A = np.array([1,2,2,3,3,5,5,7]).reshape(-1,2)
    U, B, V = GKL_bidiagonalization(A)
    B = scsp.diags(B, [0, -1]).toarray()
    print('original\n', np.array(A))
    print('reconstructed\n', U @ B @ V.H)