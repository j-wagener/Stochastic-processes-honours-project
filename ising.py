from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from random import randint , random

@njit
def exp(x):
    return np.float( np.exp(x) )

@njit
def sum_rows(matrix):
    l = len(matrix)
    K = [] 
    t = len(matrix[0])
    for i in range (t):
        r=0
        for j in range (l):
            r += matrix[j][i]
        K.append(r)
    return K

@njit
def sum_(vector):
    l = len(vector)
    s = 0
    for i in range (l):
        s+= vector[i]
    return s

@njit
def sum_matrix(matrix):
    rows = sum_rows(matrix)
    s = sum_(rows)
    return s

# print( sum_matrix( np.array([ [1,1] , [2,2]     ] )  )  )

@njit
def one_initial():
    value = 2*randint(0, 1)-1
    return np.int(value)

@njit
def rand_row(L):
    row=[0]*L
    for i in range (len(row)):
        row[i] = one_initial()
    return row

@njit
def lalala():
    value = random()
    return np.float(value)
print(lalala())

@njit
def ran_el(L):
    x = randint(0,L-1)
    return int(x)

@njit
def lattice_initial(L):
    lattice = []
    # lattice = np.zeros(L)
    for i in range (L):
        lattice.append(rand_row(L))
        # lattice[i] = rand_row()
    return lattice


@njit
def update(N,lattice,x,y):
    matrix = []

    for r in range(0, N):
        matrix.append([lattice[r][c] for c in range(0, N)])

    matrix[x][y] = - lattice[x][y]

    return matrix


@njit
def one_flip(L,Temp):
    lattice_ = lattice_initial(L)
    U=0
    Time = int(1e6)
    for i in range (Time):
        x,y = ran_el(L) , ran_el(L)
        l,r,u,d = 0,0,0,0
        if x == len(lattice_)-1:
            r=0
        else:
            r= x+1

        if y == len(lattice_)-1:
            u=0
        else:
            u= y+1


        l = x-1 #python indexing does these bcs for us
        d = y-1
        # sbar = (lattice_[l][y] + lattice_[r][y] + lattice_[x][u] + lattice_[x][d] )/4
        U = 2*1*lattice_[x][y] *(lattice_[l][y] + lattice_[r][y] + lattice_[x][u] + lattice_[x][d]) 
                        #+ lattice_[l][u] + lattice_[r][u] + lattice_[l][d] + lattice_[r][d] )  #for diagnals too
        U = float(U)
        kk = lattice_
        if U<0:
            lattice_ = update(L,kk, x,y)
        
        else:
            if exp(-U/Temp) > lalala():
                lattice_ = update(L,kk, x,y)

    return lattice_


def plotter(lattice_):   
    L = len(lattice_)
    for i in range (L):
        for j in range (L):
            if lattice_[i][j] ==1:
                plt.arrow(i,j,0,0.2 , head_width=0.2, color='red')
            else:
                plt.arrow(i,j,0,-0.2 , head_width=0.2 , color='black')
    # plt.title('Plot showing a random initial \n state of a 30x30 lattice')
    plt.title("Plot showing the 30x30 lattice after $10^6$ iterations \n T=1000")
    plt.show()

    return

plotter(np.array(one_flip(30,1e1000 )))
@njit
def mag_plot():
    M = []
    T = []
    for temp in range (1,48):
        print(temp)
        L = 10
        temp = temp/12
        m_i = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_j = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_k = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_l = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_m = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_n = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_q = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_r = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        m_s = abs( sum_matrix( np.array( one_flip(L, temp))) ) / (L**2)
        MM = [m_i,m_j,m_k,m_l,m_m,m_n,m_q,m_r,m_s]
        if temp<1:
            M.append(max(MM))
        else:
            M.append((m_i+m_j+m_k+m_l+m_m+m_n+m_q+m_r+m_s)/9)
        T.append(temp)

    return T,M
# mag_plt = mag_plot()
# plt.plot(mag_plt[0], mag_plt[1] , 'r.' )
# plt.xlabel("temperature")
# plt.ylabel("magentisation proportion")
# plt.title("Plot of how increasing the temperature \n the net magnetisation decreases \n eight nearest neighbours")

# plt.show()