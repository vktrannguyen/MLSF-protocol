# MLSF-protocol

![Protocol-Workflow](https://github.com/vktrannguyen/MLSF-protocol/blob/main/Figure1-protoc.png)

You will find herein the code (Python scripts, Jupyter notebooks), the input and output files related to our Nature Protocols paper:

**Tran-Nguyen, V. K., Junaid, M., Simeon, S. & Ballester, P. J. A practical guide to machine-learning scoring for structure-based virtual screening. *Nat. Protoc.* (2022)**

The code is found in the **Protocol_Code** folder.

Here we provide examples for three targets: **ACHE** (acetylcholinesterase), **HMGR** (HMG-CoA reductase), and **PPARA** (peroxisome proliferator-activated receptor alpha). The input and output files for a target are found in the corresponding folder.

All Supporting Information files are found in the **Supporting_Information** folder.

Two environments have to be set up in order to run the code of this protocol: **DeepCoy-env** and **protocol-env**. The **DeepCoy-env** environment can be installed according to DeepCoy authors (https://github.com/fimrie/DeepCoy). The **protocol-env** environment can be set up by using the **protocol-env.yml** file provided here as follows:

*conda env create -f protocol-env.yml;*

*conda activate protocol-env*

The Anaconda installer v4.13.0 is needed to set up these environments. Installation instructions for Anaconda can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html.

Languages: Python v3.7, Jupyter notebook, Bash.

All other necessary information is included in our protocol.

For further queries, please contact **Dr. Viet-Khoa Tran-Nguyen** (khoatnv1993@gmail.com) or **Dr. Pedro J. Ballester** (p.ballester@imperial.ac.uk).
