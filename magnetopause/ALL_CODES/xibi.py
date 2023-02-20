import pytools as pt
import numpy as np
import pandas as pd
from myutils import get_vlsvfile_fullpath, save, sidecar
import matplotlib.pyplot as plt
#matplotlib.use('module://drawilleplot')
import argparse


@sidecar
def get_b_xaxis(run, fileIndex, save = False):
    vlsvfile = get_vlsvfile_fullpath(run, fileIndex) 
    f = pt.vlsvfile.VlsvReader(vlsvfile)
    bvec = f.read_variable('fg_b')
    RE = 6371000.
    x = np.linspace(f.read_parameter('xmin'),f.read_parameter('xmax'), bvec.shape[0] ) / RE
    i = np.where(x > 0)
    xi = x[i]
    B = (bvec[:,:,:,0]**2 + bvec[:,:,:,1]**2 + bvec[:,:,:,2]**2 )**0.5 
    Bi = B[i,int(B.shape[1]/2),int(B.shape[2]/2)]        # isolate the x-axis where x>0
    Bi = Bi.flatten()
    #save x [meters], B [Tesla] along the +x axis
    df = pd.DataFrame(data={'x':xi, 'B':Bi})
    if save:
        df.to_csv('xibi_{}_{}.csv'.format(run, fileIndex), index=False)
    return df


    ##READ DATA LATER:
    #import pandas as pd
    #df = pd.read_csv('xibi.csv')
    
    
    #plot the data
    
    #plt.plot([1,10],[1e-4,1e-7])
    #plt.scatter(df['xi'], df['Bi'])
    #plt.xscale('log',base=10) 
    #plt.yscale('log',base=10) 
    
    #plt.savefig('B_r.png')


if __name__ == '__main__':

    run = 'EGL'
    fileIndex = 1760
    
    # Input parameters
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-nproc', default=1, help="number of processors to use " )
    global ARGS
    ARGS = parser.parse_args()
    
    
    dummy = get_b_xaxis(run, 1760, save = True)  #save a .csv file
    
    b_list = []
    
    times = list(range(1200,1761))
    #times = list(range(641,1761))   #test
    
    def get_b_pool(fileIndex):
        run = 'EGL'
        return get_b_xaxis(run, fileIndex, save=False)
    
    ## Parallel processing
    
    
    from multiprocessing import Pool
    pool = Pool(int(ARGS.nproc))
    b_list = pool.map(get_b_pool, times)
    pool.close()
    pool.join()
    
    
    
    #for fileIndex in range(641,times):
    #    b_list.append(get_b_axis(run, fileIndex))
    
    B_ave = b_list[0]['B'] * 0
    
    nb = len(b_list)
    for i in range(nb):
        B_ave += b_list[i]['B']
    
    B_ave = B_ave / nb
    
    df = pd.DataFrame(data={'x':np.array(b_list[0]['x']), 'B':np.array(B_ave)})
    df.to_csv('xibi_ave.csv', index=False)
    save('xibi.pck' , b_list=b_list, times = times  )
    
    
    
    
    
    

