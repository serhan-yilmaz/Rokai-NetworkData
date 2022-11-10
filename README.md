# Rokai-NetworkData
A repository for the data preprocessing steps to prepare the network data for RokaiApp.

### Step 1 - PhosphoSitePlus Data
- Download the ``Kinase_Substrate_Dataset.gz`` and ``Phosphorylation_site_dataset.gz`` files from [PhosphoSitePlus website](https://www.phosphosite.org/staticDownloads) and extract its contents to the ``/data/big/phosphositeplus/`` folder. 
- Take note of the last modified dates of the extracted files and set the ``psp_version`` variable accordingly.
- Run the ``demo_01_parse_phosphositeplus.m`` script on Matlab. 
