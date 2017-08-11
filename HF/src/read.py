def ReadHF():
    import numpy as np 
    data=[]
    with open("MULBERSHG.DAT") as f:
        for line in f:
            strData=line.split()
            for elem in strData:
                data.append(float(elem))

    S=np.zeros((6,6))
    H=np.zeros((6,6))
    G2=np.zeros((6,6,6,6))

    K=6
    iData=0
    for i in range(K):
        for j in range(K):
            S[i][j]=data[iData]
            iData=iData+1

    for i in range(K):
        for j in range(K):
            H[i][j]=data[iData]
            iData=iData+1

    for i in range(K):
        for j in range(K):
            for k in range(K):
                for l in range(K):
                    G2[i][j][k][l]=data[iData]
                    iData=iData+1

    return S,H,G2
