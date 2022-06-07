import numpy as np

def read_vsc_data(filename):  
    
    # Reads data from a file created by vlsvintpol
    # which extracted data from a Vlasiator run at different time steps
    # at one or several given points in the simulation
    
    dict_AMR_version = {'proton/vg_rho':'proton/rho','proton/vg_temperature':'proton/temperature',
                        'proton/vg_V.x':'proton/V.x','proton/vg_V.y':'proton/V.y','proton/vg_V.z':'proton/V.z',
                        'vg_b_vol.x':'B.x','vg_b_vol.y':'B.y','vg_b_vol.z':'B.z'}
    
    with open(filename) as f:
        
        # Read the header to get the variable names
        header = f.readline()[1:]
                
        # Read all data into a large array
        all_data = np.loadtxt(f, skiprows=0)

        # Get the number of virtual spacecraft grouped in the file
        n_sc = 1
        while (all_data[n_sc,0] - all_data[n_sc-1,0]) < 0.2:
            n_sc = n_sc + 1
            
        
        dict_data = {}
        nlines = all_data.shape[0]
                
        for i in range(0,len(header.split())):    
            if header.split()[i] in dict_AMR_version:
                print(dict_AMR_version[header.split()[i]])
                dict_data[dict_AMR_version[header.split()[i]]] = np.reshape(all_data[:,i], (int(nlines/n_sc),n_sc))
            else:
                dict_data[header.split()[i]] = np.reshape(all_data[:,i], (int(nlines/n_sc),n_sc))
            
        return dict_data,n_sc
