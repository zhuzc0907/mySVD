import numpy as np

a1 = np.array([1,0,1,1,0,1])
a2 = np.array([1,0,0,1,0,1])
a3 = np.array([1,0,0,0,0,1])

a1 = a1.reshape((3,-1))
a2 = a2.reshape((3,-1))
a3 = a3.reshape((3,-1))

aa = a1 + a2 +a3

# apply svd

u1, s1, vh1 = np.linalg.svd(a1, full_matrices=True)
u2, s2, vh2 = np.linalg.svd(a2, full_matrices=True)
u3, s3, vh3 = np.linalg.svd(a3, full_matrices=True)

u, s, vh = np.linalg.svd(aa, full_matrices=True)

print(u1)
print(s1)
print(vh1)

print('--------')

print(u2)
print(s2)
print(vh2)

print('--------')

print(u3)
print(s3)
print(vh3)

print('--------')

print(u)
print(s)
print(vh)

print('--------')
print(s1+s2+s3)