# RetroBioCat
RetroBioCat is a web-based tool for designing biocatalytic cascades and reactions.  

We recommend using retrobiocat through the online version hosted at https://retrobiocat.com  

However, you may run your own instance of RetroBioCat by following the installation instructions below.
Currently only exemplar data-sets for the reaction rules and for the substrate specificity database are provided.

For more information, please see our preprint:  
[Finnigan, William; Hepworth, Lorna J.; Turner, Nicholas J.; Flitsch, Sabine (2020): RetroBioCat: Computer-Aided Synthesis Planning for Biocatalytic Reactions and Cascades. ChemRxiv. Preprint. https://doi.org/10.26434/chemrxiv.12571235.v1](https://chemrxiv.org/articles/preprint/RetroBioCat_Computer-Aided_Synthesis_Planning_for_Biocatalytic_Reactions_and_Cascades/12571235?fbclid=IwAR3PBVXF-MGavQ2ejq3gdVQhdRxUYxNLPvI-EozTnqfm1Ut9R2eDJOD6i4I)  

## Option 1 - Quick start using docker
Option 1 requires that you have docker, docker-compose and git installed on your computer.  

 - Clone this repository and move directory to /retrobiocat/docker/
```
git clone https://github.com/willfinnigan/retrobiocat.git 
cd retrobiocat/docker/
```

- Build the docker containers
```
docker-compose build --no-cache
```

- Run
```
docker-compose up
```

RetroBioCat should now be available locally at http://127.0.0.1:5000  
Databases must now be initalised (see below)

## Option 2 - Manual Installation
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

* Install retrobiocat_web along with requirements
```
pip install -e .
```

### 2b. Running redis and mongodb 
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

### 2c. Running RetroBioCat 
To run the RetroBioCat website, two python scripts are required.  
From the retrobiocat directory, run (in separate terminals):

```
python retrobiocat_web/main.py
python retrobiocat_web/worker.py
```

RetroBioCat should now be available locally at http://127.0.0.1:5000



## Initialiase the database (required for both methods of installation)
Before your local version of RetroBioCat can be used, the databases it relies on must be set up.  

To do this, first login using the default admin account:  
- username: admin  
- password: password  

Navigate to the Initialise Database page in the admin menu.

Initialise the database by uploading the required files.  This can be done one at a time (recommended) or all together.

Files are available at:   
https://figshare.com/articles/software/RetroBioCat_database_files/12696482  

Reaction rules: trial_rule_set.yaml  
Activity: trial_biocatdb_will_and_lorna.xlsx  
Building blocks: building_blocks.db

Currently only example sets of reaction rules and substrate specificity information are provided, pending future publications.

Once the databases are initialised RetroBioCat is ready to use.

## Automated testing of pathway test-set  
Our publication on RetroBioCat features an evaluation on a test-set of 52 pathways.  
We automated this evaluation using a script available in the /scripts/pathway_testing/ folder.  

To run the pathway_eval.py script, install retrobiocat via option 2 (above) and ensure that your mongodb instance is running and that the databases have been initialised as described above.  

Move directories to /scripts/pathway_testing/ , and run python pathway_eval.py  

Note this script takes a long time to run.  Results are saved by default to test_pathways.xlsx

(Note, replication of the results in the paper requires the complete set of reaction rules and database file, which are not yet publicly available)  







