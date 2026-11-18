import numpy as np
from math import *

# ===========================================
def math_get_sort_order (a):
# ===========================================
    sort_indices = sorted(list(range(len(a))), key=lambda x: a[x])
    return sort_indices

# ===========================================
def math_concatenate_arrays (a,b,tol=1.0e-8):
# ===========================================
    c = np.concatenate((a, b), axis=None)
    info = np.zeros (c.size,dtype='int')
    for i in range (c.size):
        info [i] = i
    order     = math_get_sort_order (c)

    done      = False
    first_val = c [order[0]]
    k         = 0
    info [order [0]] = 0
    count     = 1

    #info [index_in_c] ==> index in new condensed array

    while k < c.size-1:
        k =  k + 1
        next_val = c [order[k]]
        if abs(next_val-first_val) < tol:
            info[order[k]] = count-1
        else:
            first_val = next_val
            info[order[k]] = count
            count = count + 1


    d = np.zeros (count)
    d_info = []
    for i in range (count):
        d_info.append ([])

    for i in range (c.size):
        #pos in new list
        k = info [i]
        d_info [k].append (i)
        d [k] = c [i]

    #be careful because d_info operates on concatneted list a+b
    #hence each element in d_info let say   c (integer index)
    # if c < a.size ===> a[c]
    #else                b[c-a.size]
    return d, d_info


#=====================================================
def main():
#=====================================================
    b = [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0]
    a = [-3.0, -4.0, -5.0]

    print("a",a)
    print("b",b)
    c, c_info = math_concatenate_arrays (b,a,1.0e-6)
    print(c)
    print(c_info)


if __name__ == '__main__':
    main()