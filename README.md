# SIMSEARCH
Current tool is aimed for similarity searching between two databases, one query and another target databases,through different fingerprints with the help of Tanimoto similarity. Query database: The experimental property of the database compounds is unknown, e.g., molecules obtained from virtual screening. Target database: The experimental property of the database compounds are known, e.g., compounds obtained from ChEMBL or DrugBank. The tool is dependent on RDKit and Pandas. Follow the instruction provided in the below lik to install RDKit. https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md Windows users may follow the following steps:

    Download and install Anaconda
    Open Anaconda prompt and type 'conda create -c conda-forge -n my-rdkit-env rdkit'
    The type 'conda activate my-rdkit-env'
    Under 'my-rdkit-env', 'conda install -c conda-forge rdkit'
    Under same environment install pandas by the command 'pip install pandas'
    Next, move to the directory where the downloded files of SIMSEARCH are placed and type 'python similarity_search.py'
