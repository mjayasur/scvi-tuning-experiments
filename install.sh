rm -rf scvi-tools
git clone https://github.com/YosefLab/scvi-tools.git
cd scvi-tools && git checkout michael/autotune && pip install .
cd ..
pip install scanpy
pip install ray
pip install tensorboardX
wget https://ndownloader.figshare.com/files/24539828 -O pancreas.h5ad