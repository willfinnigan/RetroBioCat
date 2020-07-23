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

## Running redis and mongodb 
RetroBioCat requires access to a redis server and a mongo database on the default ports.  
We recommend using docker to run redis and mongodb.  

To run redis using docker:
```
docker run -d -p 6379:6379 redis
```

To run mongodb using docker:
```
docker run -d -p 27017-27019:27017-27019 mongo:4.0.4
```

## Running RetroBioCat 
To run the RetroBioCat website, two python scripts are required.  
From the retrobiocat directory, run (in separate terminals):

```
python retrobiocat_web/main.py
python retrobiocat_web/worker.py
```

RetroBioCat should now be available locally at http://127.0.0.1:5000

## Initialiase the database
Before your local version of RetroBioCat can be used, the databases it relies on must be set up.  

To do this, first login using the default admin account:  
- username: admin  
- password: password  

Navigate to the Initialise Database page in the admin menu.

Initialise the database by uploading the required files.  This can be done one at a time (recommended) or all together.

Files are available at:  
Currently only example sets of reaction rules and substrate specificity information are provided, pending future publications.



