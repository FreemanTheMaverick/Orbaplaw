import numpy as np


def ConjugateGradient(lasth,g,lastg,cgtype):
    gamma=0
    if lasth is None or lastg is None:
        return g
    if cgtype=="FR":
        gamma=np.linalg.norm(g)/np.linalg.norm(lastg)
    elif cgtype=="PR":
        gamma=np.sum(g*(g-lastg))/np.linalg.norm(lastg)
    h=g+gamma*lasth
    if np.sum(h*g)<0:
        h=g
    return h
