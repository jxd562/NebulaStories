import numpy as np

def silly_sorry(xvec):
    """Hardcoded answer for the parse challenge"""

    Granite = xvec[0]
    Preston = xvec[1]
    Blue = xvec[2]
    Thomas = xvec[3]
    Checkers = xvec[4]
    Fluffy = xvec[5]
    Colors = xvec[6]
    Chancy = xvec[7]
    Snicks = xvec[8]
    Daphnie = xvec[9]
    Wiskers = xvec[10]
    Socks = xvec[11]

    t = np.array([[Granite,   Preston,   Blue,      Thomas, 0],
                  [Checkers,  Thomas,    Preston,   Fluffy, 0],
                  [Checkers,  Checkers,  Socks,     Thomas, 0],
                  [Fluffy,    Preston,   Checkers,  Thomas, 0],
                  [Thomas,    Colors,    Checkers,  Chancy, 0],
                  [Chancy,    Snicks,    Checkers,  Preston,   Preston],
                  [Wiskers,   Fluffy,    Granite,   Daphnie, 0],
                  [Fluffy,    Preston,   Chancy, 0, 0]])
    return np.sum(t, axis=-1)

if __name__ == '__main__':
    result = silly_sorry(np.arange(12))
    print('result is \n')
    print(result)
