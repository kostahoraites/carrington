import numpy as np

def read_vsc_data(filename, n_sc = None, tvar='t'):  
    
    # Reads data from a file created by vlsvintpol
    # which extracted data from a Vlasiator run at different time steps
    # at one or several given points in the simulation
    
    dict_AMR_version = {'proton/vg_rho':'proton/rho','proton/vg_temperature':'proton/temperature',
                        'proton/vg_V.x':'proton/V.x','proton/vg_V.y':'proton/V.y','proton/vg_V.z':'proton/V.z',
                        'vg_b_vol.x':'B.x','vg_b_vol.y':'B.y','vg_b_vol.z':'B.z'}
    
    with open(filename) as f:

        if filename[-4:] == '.csv':
            delimiter = ','
        else:
            delimiter = None
        
        # Read the header to get the variable names
        header = f.readline()[1:]
        varnames = header.split(delimiter)
        if varnames[-1][-1] == "\n":
           varnames[-1] = varnames[-1][:-1]     # get rid of new line character in last field name

        converters = {}
        for c in range(len(varnames)):
           converters[c] = lambda s : float(s.strip() or np.nan)  #handle empty data as nan

        # Read all data into a large array
        all_data = np.loadtxt(f, skiprows=0, delimiter=delimiter, converters=converters)

        if n_sc is None:
        # Get the number of virtual spacecraft grouped in the file
            t_ind = np.where(np.array(varnames) == tvar)
            n_sc = 1
            while (all_data[n_sc,t_ind] - all_data[n_sc-1,t_ind]) < 0.2:
                n_sc = n_sc + 1
            #n_sc = 1
            #while (all_data[n_sc,0] - all_data[n_sc-1,0]) < 0.2:
            #    n_sc = n_sc + 1

        print(varnames)
        print(tvar)
        print(n_sc)

        dict_data = {}
        nlines = all_data.shape[0]

        for i in range(0,len(varnames)):    
            if varnames[i] in dict_AMR_version:
                print(dict_AMR_version[varnames[i]])
                dict_data[dict_AMR_version[varnames[i]]] = np.reshape(all_data[:,i], (int(nlines/n_sc),n_sc))
            else:
                dict_data[varnames[i]] = np.reshape(all_data[:,i], (int(nlines/n_sc),n_sc))
            
        return dict_data,n_sc
