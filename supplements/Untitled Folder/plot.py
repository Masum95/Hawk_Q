from matplotlib  import pyplot as plt 
import pandas as pd 
import numpy as np 
import math
import random 


#df = pd.read_csv('kmerdata303982630', header=None, sep=',')
#plt.scatter(x=np.log(np.array(df[0])), y=np.array(df[1]))
# print( np.log(np.array(df[0])) )
#plt.show()



xin =  pd.read_csv('pvals_top.txt', header=None, sep=',')
yin =  pd.read_csv('pvals_top2.txt', header=None, sep=',')

random_size = 500
print( min(random_size, len(xin)) )
random_sample = random.sample(list(zip(xin[0],yin[0])), min(random_size, len(xin) ))
xx,yy = list(zip(*random_sample))

plt.scatter(x=np.negative(np.log(np.array(xx))), y=np.negative(np.log(np.array(yy))))

plt.show()

