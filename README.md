# MLSF-protocol

![Protocol-Workflow](https://github.com/vktrannguyen/MLSF-protocol/blob/main/Figure-1_Nat-Protoc.png)

You will find herein the code (bash scripts, Python scripts, Jupyter notebooks), the input and output files related to our Nature Protocols paper:

**Tran-Nguyen, V. K., Junaid, M., Simeon, S. & Ballester, P. J. A practical guide to machine-learning scoring for structure-based virtual screening. *Nat. Protoc.* **18**, 3460–3511 (2023)**

Here we provide examples for three targets: **ACHE** (acetylcholinesterase), **HMGR** (HMG-CoA reductase), and **PPARA** (peroxisome proliferator-activated receptor alpha). The input and output files for a target are found in the corresponding folder.

Inside each of the three folders **ACHE**, **HMGR**, **PPARA**, you will find the following sub-folders:

- **DEKOIS2.0** sub-folder: input and output files for Section A of this protocol.
- **Own_data** sub-folder: input and output files for Sections B, C, D of this protocol. This sub-folder contains:
  - **SMILES** sub-folder: SMILES strings of the users' own true actives, true inactives and decoys.
  - **ChemAxon** sub-folder: raw data for Figure 4 in the manuscript and Figures S1-S2 in Supporting Information.
  - **MLSF_PETS** sub-folder: input and output files related to training-test partitions obtained from the "Pre-Existing Test Set" (PETS) option.
  - **MLSF_OTS** sub-folder: input and output files related to training-test partitions obtained from the "Own Test Set" (OTS) option.
  - *data.xlsx*: the master Excel file where all necessary information on the users' own true actives and true inactives is stored.
- **D-test-sets** sub-folder: files related to the "Dissimilar" (D) test sets of each target's PETS and OTS (Tables 2, 3).

The code is found in the **Protocol_Code** folder:

- *bash-commands.zip*: zip file containing all bash commands for the compilation of scores from Smina, CNN-Score, RF-Score-VS, IFP (Steps 14, 21, 30, 41).
- *Precision-Recall-curve.ipynb*: Jupyter notebook for plotting the precision-recall (PR) curve.
- *EF1-NEF1.sh*: bash script for computing the EF1% and NEF1% values.
- *vlookup.awk* and *potency.sh*: necessary files for extracting the potency of true hits among the top 1%-ranked molecules.
- *Morgan-fp-simil.ipynb*: Jupyter notebook for computing the similarity of Morgan fingerprints.
- *Compound-clustering_Morgan-fp.ipynb*: Jupyter notebook for clustering molecules based on their Morgan fingerprint similarity.
- *Remove_AVE_Python3.py*: Python code for splitting a data set into four subsets (training actives, training inactives, test actives, and test inactives) in an unbiased manner.
- *MLSFs.ipynb*: Jupyter notebook for training and evaluating target-specific SFs.

**ATTENTION**: if you use the scripts *EF1-NEF1.sh* and *potency.sh* to process the csv hit lists issued by the *MLSFs.ipynb* Jupyter notebook (in **Section D**):

- You must remove the "Predicted_Class" column (while keeping the "Real_Class" column as is) from the hit lists beforehand.
- The *EF1-NEF1.sh* code must be slightly modified as follows (you can either edit it using the *vi*/*vim* command or open it in a text editor, e.g. Notepad++, and save it after modifying):

```
#Count the number of true active molecules (true hits) in the whole test set:
A=$(grep -c 'Active' $hitlistname)-1
A=$((A))
```

All Supporting Information files are found in the **Supporting_Information** folder:

- *Supporting-Information_MLSFs-SBVS.docx*: Tables S1-S5, Figures S1, S2 and other supplementary information as indicated in the manuscript.
- *Supporting-Information_MLSFs-SBVS_DEKOIS-retrieved-actives.xlsx*: raw data for Figure 3 in the manuscript.
- *Supporting-Information_MLSFs-SBVS_10-runs.xlsx*: virtual screening performance of five learning algorithms across 10 training-test runs on five training-test partitions (full partitions and "Dissimilar" test sets), raw data for Figures 5A, 5C in the manuscript.

Two environments have to be set up in order to run the code of this protocol: **DeepCoy-env** and **protocol-env**. The **DeepCoy-env** environment can be installed according to DeepCoy authors (https://github.com/fimrie/DeepCoy). The **protocol-env** environment can be set up by using the *protocol-env.yml* file provided here as follows:

```
conda env create -f protocol-env.yml
conda activate protocol-env
```

A part of the code for Step 60 (to retrieve SMILES strings from PubChem) is provided here (you can simply copy and do not have to retype):

```
ipython
import pandas as pd
import pubchempy as pcp
```

```
cid_list = df['cid']
smiles_list = []
def get_smiles(input_cid):
    mol = pcp.Compound.from_cid(int(input_cid))
    smiles_list.append(mol.canonical_smiles)
    return smiles_list
```

```
[get_smiles(mol) for mol in cid_list]
df['smiles']= smiles_list
```

The Anaconda installer v4.13.0 is needed to set up these environments. Installation instructions for Anaconda can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html.

Several parts of the code/scripts/Jupyter notebooks used in this protocol were developed from the original code accessible in the following github repositories:

- https://github.com/sawsimeon/MLSF-PDL1
- https://github.com/sawsimeon/PDL1_Generic

Programming languages used in this protocol: Python v3.7, Bash.

All other necessary information is included in our protocol.

For further queries, please contact **Dr. Viet-Khoa Tran-Nguyen** (khoatnv1993@gmail.com) or **Dr. Pedro J. Ballester** (p.ballester@imperial.ac.uk).
