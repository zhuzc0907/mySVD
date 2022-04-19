import numpy as np

def golub_kahan(a):
  n = a.shape[1]
  v = np.ones(n, dtype="float64") / np.sqrt(n)
  u = np.zeros(a.shape[0], dtype="float64")
  beta = 0
  U, V = np.zeros_like(a, dtype="float64"), np.zeros((n,n), dtype="float64")

  alpha_list = []
  beta_list = []

  for i in range(n):
    V[:, i] = v
    u = a @ v - beta * u
    alpha = np.linalg.norm(u)
    u /= alpha
    U[:, i] = u
    v = a.T @ u - alpha * v
    beta = np.linalg.norm(v)
    v /= beta

    alpha_list.append(alpha)
    beta_list.append(beta)

  return U, V, alpha_list, beta_list

def build_bid_mat(alpha_list, beta_list):
  n = len(alpha_list)

  B = np.zeros(n*n).reshape(n, n)

  for i in range(n):
    B[i,i] = alpha_list[i]
    if(i!=n-1):
      B[i,i+1] = beta_list[i]

  return B


def bid_validation(a, U, B, V): 

  A_re = np.matmul(np.matmul(U,B), V.T)

  return np.linalg.norm(A_re-a)
  # return A_re


# A = np.ones(100).reshape(-1,5)
A = np.array(range(12)).reshape(-1,3)
U, V, alpha_list, beta_list = golub_kahan(A)
B = build_bid_mat(alpha_list,beta_list)

# compute B_re
B_re = U.T @ A @ V

error = bid_validation(A, U, B, V)


# Check SVD of A and B
U_A, S_A, Vh_A = np.linalg.svd(A)
U_B, S_B, Vh_B = np.linalg.svd(B)
U_B_re, S_B_re, Vh_B_re = np.linalg.svd(B_re)


error_svd = np.linalg.norm(S_A - S_B)


print('Input A:')
print(A)
print('---------')
print('B constructed by alpha and beta:')
print(B)
print('---------')
print('B reconstructed by U and V')
print(B_re)
print('---------')
print('Errors between A and B: Bid_error = ' +str(error))
print('---------')
print('Singluar Values of A')
print(S_A)
print('---------')
print('Singluar Values of B')
print(S_B)
print('---------')
print('Singluar Values of B_re')
print(S_B_re)
print('---------')
print('Errors between Singluar Values of A and B: Svd_error = ' + str(error_svd))
