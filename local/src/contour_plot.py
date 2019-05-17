import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if len(sys.argv) < 2:
    print("Missing Arguments")
    exit(1)

if os.path.isfile(sys.argv[1]) == False:
    print("FileNotFoundException")
    exit(1)

segcopy = pd.read_csv(sys.argv[1], sep="\t")
segcopy_matrix = segcopy.iloc[:, 4:].copy().values
xlist = np.arange(0, len(segcopy_matrix[0]))
ylist = np.arange(0, len(segcopy_matrix))

X,Y = np.meshgrid(xlist,ylist)
plt.figure()

cp = plt.contourf(X, Y, segcopy_matrix)
plt.colorbar(cp)
plt.title('CNVs Contours Plot')
plt.xlabel('cells')
plt.ylabel('positions')
plt.savefig(sys.argv[2] + '-sc/contoutCN.png')

