# RNAit
Reimplementation of RNAit RNAi primer selection tool

RNAit is a tool for selecting RNAi primers for Trypanosome genomes, originally
implemented in Perl (see [Redmond et al. (2003) RNAit: an automated web-based
tool for the selection of RNAi targets in Trypanosoma brucei, Mol. Biochem.
Parasitol.](https://www.sciencedirect.com/science/article/pii/S0166685103000458?via%3Dihub))

This is a reimplementation to replace the now offline original deployment

## Setting up a development instance (MacOS)

bash 4 is required for helper scripts, whereas the MacOS default is 3.2, so grab a newer version with brew...

```bash
brew install bash
```

Create a conda environment, and add the RNAIT_ROOT variable to point to the location of the RNAit software....

```bash
conda create -n RNAit
mkdir -p ~/miniconda3/envs/RNAit/etc/conda/activate.d
mkdir -p ~/miniconda3/envs/RNAit/etc/conda/deactivate.d
echo '#!/bin/env bash' > ~/miniconda3/envs/RNAit/etc/conda/activate.d/vars.sh
echo "export RNAIT_ROOT=/path/to/RNAIT/" >>  ~/miniconda3/envs/RNAit/etc/conda/activate.d/vars.sh
echo '#!/bin/env bash' > ~/miniconda3/envs/RNAit/etc/conda/deactivate.d/vars.sh
echo "unset RNAIT_ROOT" >>  ~/miniconda3/envs/RNAit/etc/conda/deactivate.d/vars.sh
source activate RNAit
```

We need to install nginx and uwsgi from conda, and symlink the nginx config in place...
Nodejs is also required from conda...

```bash
conda install nginx uwsgi nodejs
mv $CONDA_PREFIX/etc/nginx/sites.d/default-site.conf $CONDA_PREFIX/etc/nginx/sites.d/default-site.conf.hiding
ln -s ${RNAIT_ROOT}/etc/nginx-site.conf $CONDA_PREFIX/etc/nginx/sites.d/
```

The node.js broken-link-checker is used for verifying HTML pages, while vnu (available via brew) is jused for HTML syntax validation
```bash
npm install -g broken-link-checker
brew install vnu
brew install java
```

Syntax and link checking is carried out using bin/check_pages.sh. Ensure the correct path to the vnu installation is set in the script, since this will vary with brew versions...


The nginx and uwsgi instances can be started and stopped using:
```bash
bin/start_servers.sh
bin/stop_servers.sh
```

UWSGI can be made to reload the python scripts when these are modified by touching the 'reload' file within the uwsgi directory

