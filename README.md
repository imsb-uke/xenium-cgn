# xenium-cgn
Subcellular resolution immune landscape of glomerular tranjectories in auto immune kidney diseases


Please use environment.yml to create environment. Additional incompatible packages should have a Dockerfile

```conda env create --name <environment_name> -f environment.yml```

or with the environment path 

```conda env create --prefix <environment_path> -f environment.yml```

environment path can be similar to the following: ```~/envs/<environment_name>```

## create python kernel for jupyter lab

```bash
conda activate <environment_name/path>
ipython kernel install --user --name=xen-cgn-py
conda deactivate
```
Restart jupyterlab to have it available

## create R kernel for jupyter lab
```bash
conda activate <environment_name/path>
R
```
An R shell is started
```R
install.packages('IRkernel')
IRkernel::installspec(name = 'xen-cgn-r', displayname = 'xen-cgn-r')
quit()
```



###
ZS: I had to install shapely to define periglom regions. Making note here to keep track of additionally needed packages.

#pip install shapely
Output : 
Collecting shapely
  Downloading shapely-2.0.5-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl.metadata (7.0 kB)
Requirement already satisfied: numpy<3,>=1.14 in /opt/conda/envs/xen-cgn/lib/python3.11/site-packages (from shapely) (1.24.4)
Downloading shapely-2.0.5-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (2.5 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 2.5/2.5 MB 7.8 MB/s eta 0:00:0000:0100:01
Installing collected packages: shapely
Successfully installed shapely-2.0.5
Note: you may need to restart the kernel to use updated packages.
#####
