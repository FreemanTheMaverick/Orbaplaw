import numpy as np


def ReadMwfnMat(file): # Very ugly codes to parse matrices in very ugly format. Do not try to understand.
    lines=None
    with open(file,'r') as f:
        lines=f.readlines()

    matrix_heads=[]
    for iline,line in enumerate(lines):
        if "**********" in line:
            matrix_heads.append(iline+1)
    matrix_texts=[]
    for imatrix,head in enumerate(matrix_heads):
        tail=len(lines) if head==matrix_heads[-1] else matrix_heads[imatrix+1]-2
        matrix_texts.append(lines[head:tail])
    matrices=[]
    for imatrix,matrix_text in enumerate(matrix_texts):
        nbasis=0
        block_heads=[]
        for iline,line in enumerate(matrix_text):
            if '.' not in line:
                block_heads.append(iline)
                nbasis=int(line.split()[-1])
        blocks=[]
        for jblock,head in enumerate(block_heads):
            tail=len(matrix_text) if head==block_heads[-1] else block_heads[jblock+1]
            blocks.append(matrix_text[head:tail])
        matrix=np.zeros([nbasis,nbasis])
        for jblock,block_text in enumerate(blocks):
            cols=None
            for kline,line in enumerate(block_text):
                if kline==0:
                    cols=np.array(block_text[0].split(),dtype="int")-1
                else:
                    row=int(line.split()[0])-1
                    nelements=len(line.split())-1
                    matrix[row,cols[:nelements]]=np.array(line.split()[1:])
        matrices.append(matrix+matrix.T-np.diag(np.diag(matrix)))
    return matrices
