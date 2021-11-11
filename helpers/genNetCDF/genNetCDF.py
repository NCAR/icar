import numpy as np

def fixType(val):
    val_t = type(val)
    if val_t == int:
        return np.int32(val)
    elif val_t == float:
        return np.single(val)
    else:
        return val
