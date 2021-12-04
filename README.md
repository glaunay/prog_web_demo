# TP Prog Web -- Corrections

## Setup
1. Install [poetry](https://python-poetry.org/docs/)

2. Install project dependencies, by typing in the cloned repository
`poetry install`
3. get the database [sql file](https://perso.liris.cnrs.fr/pierre-antoine.champin/2019/progweb-python/_static/ensembl_hs63_simple.sqlite)
4. run flask through poetry virtual env by typing `poetry run python server.py`

## Test
The POST API can be tested with curl

```sh
curl -d 'Ensembl_Gene_ID=MY_GENE_ID&Chromosome_Name=Y&Band="q.12&Gene_Start=10&Gene_End=200' -H "Content-Type: application/x-www-form-urlencoded" -X POST http://localhost:5000/api/genes/
```