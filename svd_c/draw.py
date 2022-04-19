import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams["figure.autolayout"] = True

df = pd.read_csv("u_orth_error.csv", nrows=1000,names=["U error"])
# df = pd.read_csv("v_orth_error.csv", nrows=1000, names=["V error"])

df.plot()

plt.show()