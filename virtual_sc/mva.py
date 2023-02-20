import numpy as np



# implemented minimum variance analysis for finding the normal direction to a discontinuity (such as a shock) 
# Reference: "Minimum and Maximum Variance Analysis" BENGT U.  ̈O. SONNERUP AND MAUREEN SCHEIBLE (1998, 2000)
# for propagating the points, need to calculate the normal velocity see Andreeova et al. (2011), eq. 2

def mva(B_comp_list):
    # B_comp_list is a list of B-vectors, representing individual measurements
    # each B-vector in the list is a numpy array, e.g. [3.0 6.2] in 2D or [-2.1 3.9 0.1] in 3D
    npts = len(B_comp_list)
    ndim = B_comp_list[0].size
    M = np.zeros([ndim, ndim])

    # construct matrix M_mu_nu, Sonnerup and Scheible eq. 8.8


    # solve Sonnerup and Scheible eq. 8.7
        #l, v = np.linalg.eig(  [matrix]) # l=eigenvalues, v = eigenvectors (2D array)
    return normal_vector


def propagate_points( point_list, rho, v_comp_list, B_comp_list):





