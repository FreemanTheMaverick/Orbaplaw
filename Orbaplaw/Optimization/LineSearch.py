import numpy as np
import scipy.interpolate as si


def MeshCubic(x1,x2,v1,v2,g1,g2,func,func_para,num):
    displacement=x2-x1
    length=np.linalg.norm(displacement)
    direction=displacement/length
    ls=np.linspace(0,length,num=num)
    xs=np.linspace(x1,x2,num=num)
    vs=np.array([v1]+[func(x,*func_para) for x in xs[1:-1]]+[v2])
    gl1=np.sum(g1*direction)
    gl2=np.sum(g2*direction)
    poly=si.CubicSpline(ls,vs,bc_type=((1,gl1),(1,gl2)))
    roots=poly.derivative().solve()
    possible_ls=[root for root in roots if root>0 and root<length]+[length]
    possible_xs=[x1+direction*l for l in possible_ls]
    possible_vs=[func(x,*func_para) for x in possible_xs[:-1]]+[v2]
    best_index=np.argmax(possible_vs)
    best_x=possible_xs[best_index]
    return best_x

