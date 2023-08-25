# Thesis_MaStat

This repository contains the code for the masterthesis titled: Improving interim decisions in a Simon’s two-stage design for single arm trials. 

## Practical
To execute the simulations under settings you will need to install [R](https://www.r-project.org) with the [goldilocks](https://cran.r-project.org/web/packages/goldilocks/index.html) package, [python 3](https://www.python.org/downloads/), [parrallel](https://www.gnu.org/software/parallel/) and use the following command to run the simulations:
```bash
python3 name_of_python_file.py | parrallel -j 16 {}
```
with 16 the number of threads.
