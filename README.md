# MLMicel
The MLMicel package enables the detection of separate micelles, using an unsupervised machine learning-based algorithm. This model is particularly developed for the detection of separate polymeric micelles, but can be utilized for other molecular aggregates.

## Reference

## The funtionalities of the MLMicel:
1. Identification of the separate micelles.
2. Quantification of the structural properties of sphericacl micelles, such as micelle size, micelle core size, etc.
3. Quantification of the diffusion coefficients of micelles.
4. 4. Detection of the shape of the micelles, e.g., spherical, cylinderical, etc.

<img width="700" alt="Fig1" src="https://github.com/user-attachments/assets/52c20287-a657-415b-8fe4-feaadb5fd4d4">
<img width="700" alt="Fig2" src="https://github.com/user-attachments/assets/33e645c9-bd44-436e-be0c-7f36470d4135">



# Installation
Cloning the repository and creating a conda environment.
``` 
git clone https://github.com/arashelahi/MLMicel.git
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
```

Installing the package
```
pip install -e

```
# Using the Model

Loading the clustering model
```
from   MLMicelle.utils import traj_reader as trj
import MLMicelle.clustering as clst
```

Loading the GROMACS generated files readers.

```
from   MLMicelle.utils import traj_reader as trj
```

Detect the micelles using the spatial coordinates.

```
traj_file='coord.gro' ## trajectory file
coord_data= trj.gro_reader(traj_file)
clustered_data = clst.clustering(pd_data)
```
Measuring the structural properties of micelles, i.e., micele size, micelle core size, and aggregation number.

```
import MLMicelle.analysis as mic_anal
mice_size, core_size, agg_numb mic_anal.analysis(clustered_data)

```