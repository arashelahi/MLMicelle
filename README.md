# MLMicelle
The MLMicelle package enables the detection of separate micelles, using an unsupervised machine learning-based algorithm. This model is particularly developed for the detection of separate polymeric micelles, but can be utilized for other molecular aggregates.

## Reference
Please cite the following paper if you use the MLMicelle package:

Elahi, A.; Bidault, X.; Chaudhuri, S. Temperature-transferable coarse-grained model for poly(propylene oxide) to study thermo-responsive behavior of triblock copolymers. J. Phys. Chem. B 2022, 126, 292– 307,  DOI: 10.1021/acs.jpcb.1c06318


## The funtionalities of the MLMicelle:
1. Identification of the separate micelles.
2. Quantification of the structural properties of sphericacl micelles, such as micelle size, micelle core size, etc.
3. Quantification of the diffusion coefficients of micelles.
4. Detection of the shape of the micelles, e.g., spherical, cylinderical, etc.

<img width="700" alt="Fig1" src="https://github.com/user-attachments/assets/52c20287-a657-415b-8fe4-feaadb5fd4d4">
<img width="700" alt="Fig2" src="https://github.com/user-attachments/assets/33e645c9-bd44-436e-be0c-7f36470d4135">



# Installation
Cloning the repository and creating a conda environment.
``` 
git clone https://github.com/arashelahi/MLMicelle.git
cd MLMicelle
conda create -n MLMicelle
```
Installing the prerequisite packages
```
pip install numpy
pip install scikit-learn
pip install pandas
pip install matplotlib
pip install seaborn
pip install --upgrade MDAnalysis
```

Installing the package
```
pip install -e .

```
# Using the Model

Loading the clustering model
```
from   MLMicelle.utils import traj_reader as trj
import MLMicelle.clustering as clst
import MLMicelle.analysis as mic_anal

```

Loading the GROMACS generated files readers.

```
import MDAnalysis as mda
from   MLMicelle.utils import traj_reader as trj
```

Detect the micelles using the spatial coordinates.

```
# Path to the trajectory and topology files
traj_dir = './data/input/' 
# load the trajectory file using MDAnalysis
u = mda.Universe(traj_dir + '/md.tpr', traj_dir + '/md_clus.gro')
# Process the Trajectory Information to a pd.DataFrame
clustered_data = clst.clustering(u, bead_type='PO55') 
```
Measuring the structural properties of micelles, i.e., micele size, micelle core size, aggregation number, and anisotropic shape

```
mice_size, core_size, agg_numb, anisotropy =  mic_anal.analysis(u, clustered_data) 


```


