# RetroBioCat
RetroBioCat is a web-based tool for designing biocatalytic cascades and reactions.  

We recommend using retrobiocat through the online version hosted at https://retrobiocat.com  

You may run your own instance of RetroBioCat by following the installation instructions below.
However currently only exemplar data-sets for the reaction rules and for the substrate specificity database are provided.

For more information, please see our preprint:  
[Finnigan, William; Hepworth, Lorna J.; Turner, Nicholas J.; Flitsch, Sabine (2020): RetroBioCat: Computer-Aided Synthesis Planning for Biocatalytic Reactions and Cascades. ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv.12571235.v1](https://chemrxiv.org/articles/preprint/RetroBioCat_Computer-Aided_Synthesis_Planning_for_Biocatalytic_Reactions_and_Cascades/12571235?fbclid=IwAR3PBVXF-MGavQ2ejq3gdVQhdRxUYxNLPvI-EozTnqfm1Ut9R2eDJOD6i4I)  



## Installation
RetroBioCat requires anaconda or miniconda with python 3.7 or later

* First install the following conda packages
```
conda install -c rdkit rdkit -y
```

* Next, clone this repository and change your working directory to that of the repository  
```
git clone https://github.com/willfinnigan/retrobiocat.git 
cd retrobiocat
```

* Install the requirements
```
pip install -r requirements.txt
```

## Usage  
RetroBioCat requires access to a redis server and a mongo database on the default ports.  
We recommend using docker to run redis and mongodb.  


