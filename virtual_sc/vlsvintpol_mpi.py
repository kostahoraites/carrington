
import pytools as pt
import numpy as np
import sys
import argparse
from mpi4py import MPI

#outfile="output.txt"
#outfile="output_egi.txt"

def chunk(a, n):
    k, m = divmod(len(a), n)
    tmp=list(a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))
    return np.array(tmp).transpose().tolist()

def extract_file(filename):
    out = []
    try:
        cellids=[]
        values=[]
        f=pt.vlsvfile.VlsvReader(filename)
        f.optimize_open_file()
        t=f.read_parameter("time")
        if t == None:
	        t=f.read_parameter("t")
        if t == None:	    
                print("Unknown time format in file " + filename)
        
        for coord in coords:
            if(inRE):
                cellids.append(f.get_cellid(coord * 6371000))
            else:
                cellids.append(f.get_cellid(coord))

        for i,var in enumerate(variables):
            values.append(f.read_variable(variables[i],operator=operators[i],cellids=cellids))
        
        for i,id in enumerate(cellids):
            out_line = str(t) + " " +  ' '.join(map(str, coords[i])) + " " + str(id)
            for j,varval in enumerate(values):
                out_line = out_line +  " " + str(varval[i])
            out.append([filename, out_line])
        f.optimize_close_file()
    except Exception as e:
        print(filename,e,file=sys.stderr)
        out.append([filename, "#Could not read " + filename])
        pass
    
    return out

def parse():

    parser = argparse.ArgumentParser()
    parser.add_argument('-var', nargs='*', help="a list of variable.operator's to output, e.g., v.magnitude rho B.x " )
    parser.add_argument('-i', nargs='*', help="a list of vlsv files")
    parser.add_argument('-c', help="A file with coordinates (can also be give from stdin)")
    parser.add_argument('-re', action='store_true', help="Coordinates in RE, in meters by default")
    parser.add_argument('-outfile', default='output.txt', help="output data file")
    args = parser.parse_args()


    if args.var is None:
        #defaults
        varnames=["rho","B.magnitude","B.x","B.y","B.z"]
    else:
        varnames=args.var

    #read in variables and their operatorsx
    operators=[]
    variables=[]
    for i,var in enumerate(varnames):
        varop=var.split(".")
        variables.append(varop[0])
        if len(varop)==1:
            operators.append("pass")
        else:
            operators.append(varop[1])

    #read in coordinates
    if args.c is None:
        coords = np.loadtxt(sys.stdin, dtype=np.float)
    else:
        coords = np.loadtxt(args.c, dtype=np.float)

    #if just single point make it into array with 1 row
    coords = np.atleast_2d(coords)

    with open(args.outfile, 'w') as f:
    #with open(outfile, 'w') as f:
        if(args.re):
            f.write("#t X_RE Y_RE Z_RE CELLID " + " ".join(varnames)+"\n")
        else:
            f.write("#t X Y Z CELLID " + " ".join(varnames)+"\n")

    return args.i,coords,variables,operators,args.re, args.outfile
    #return args.i,coords,variables,operators,args.re


#init mpi
comm = MPI.COMM_WORLD
comSize = comm.Get_size()
rank = comm.Get_rank()
MASTER=0


if (rank==MASTER):

    #parse input
    #files,coords,variables,operators,inRE = parse()
    files,coords,variables,operators,inRE,outfile = parse()
    numfiles=len(files)

    #some data wrangling
    files=chunk(files,comSize)
    print("Size= ",comSize, " Files=",numfiles)
    if (numfiles<comSize):
        print("Exiting:Number of files(%i) less than MPI tasks(%i). Aborting. Rerun with less tasks or more files. :)"%(len(files),comSize),file=sys.stderr)
        sys.exit() # this is not a graceful exit 

else:

    files=None
    operators=None
    variables=None
    coords=None
    inRE=None
    outfile=None


#comms
files = comm.scatter(files, root=MASTER)
operators = comm.bcast(operators, root=MASTER)
variables = comm.bcast(variables, root=MASTER)
coords = comm.bcast(coords, root=MASTER)
inRE = comm.bcast(inRE, root=MASTER)


for taskID in range(comSize):
    if (rank==taskID):
        print("Task %i has %i file(s)"%(rank,len(files)))

#main 
data=[]
for file  in files:
    print(rank,file)
    data.append(extract_file(file))


#//gather back
finalData = []
finalData = comm.gather(data,root=0)
#write
if (rank==MASTER):
    with open(outfile, 'a') as f:
        for i in finalData:
            for j in sorted(i):
                for k in j:
                 f.write(k[1]+"\n")
    
