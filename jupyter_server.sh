#!/bin/bash -l
#SBATCH --job-name=jupyter_notebook
#SBATCH --time=12:00:00
#SBATCH --output=jupyter_notebook_README.txt
#SBATCH --error=jupyter_notebook_IPINFO.err
##SBATCH --nodes=1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -M carrington
#SBATCH -p short
#SBATCH --mem-per-cpu=8000

# get tunneling info
XDG_RUNTIME_DIR=""
node=$(hostname -s)
user=$(whoami)
cluster="turso02"
port=50002 #<SET A PORT HERE, PREFERABLY BETWEEN 50002 AND 60000>

# print tunneling instructions jupyter-log
echo -e "
# Note: below ${port} is used to signify the port.
#       However, it may be another number if ${port} is in use.
#       Check jupyter_notebook_IPINFO.err to find the port.

# Command to create SSH tunnel:
ssh -q -f -L ${port}:${node}:${port} ${user}@${cluster}.cs.helsinki.fi sleep 600
# You will have 10 minutes to connect from VSCode to the jupyter server, otherwise the tunnel closes

# Find the link with "127.0.0.1:${port}/?token=..." from the file "jupyter_notebook_IPINFO.err"
# Go to VSCode with remote connection to ${cluster}, open a jupyter notebook.
# Click the top right button that says "No kernel" or some python version.
# Click "Select another kernel"->"Existing jupyter server"->"Enter URL..."
# Paste in the link you got from the "...IPINFO.err" file, press Enter, and Enter again.
# Run the notebook normally

"

# ACTIVATE YOUR PYTHON VIRTUAL ENVIRONMENT HERE IF NEEDED
#conda activate <my env with jupyter installed>

# Run jupyter server
jupyter-notebook --no-browser --ip=${node} --port=${port}

# keep it alive
sleep 43200




