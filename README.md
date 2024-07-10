# xenium-cgn
Subcellular resolution immune landscape of glomerular tranjectories in auto immune kidney diseases


Please use environment.yml to create environment. Additional incompatible packages should have a Dockerfile

```conda env create --name <environment_name> -f environment.yml```
or
```conda env create --prefix /home/envs/<environment_name> -f environment.yml```

## create kernel for jupyter lab

```bash
conda activate <environment_name>
ipython kernel install --user --name=<kernel_name>
conda deactivate
```
Restart jupyterlab to have it available
