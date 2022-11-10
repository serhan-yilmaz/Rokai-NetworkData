# Rokai-NetworkData
A repository for the data preprocessing steps to prepare the network data for RokaiApp.

### Step 1 - PhosphoSitePlus Data
- Download the ``Kinase_Substrate_Dataset.gz`` and ``Phosphorylation_site_dataset.gz`` files from [PhosphoSitePlus website](https://www.phosphosite.org/staticDownloads) and extract its contents to the ``/data/big/phosphositeplus/`` folder. 
- Take note of the last modified dates of the extracted files and set the ``psp_version`` variable accordingly.
- Run the ``demo_01_parse_psp_kinase_substrates.m`` script on Matlab. 

### Step 2 - Signor Data
- Visit the [Signor website](https://signor.uniroma2.it/downloads.php), select ``file xls (slower)`` option and download the datasets for human/mouse/rat species. Place the contents into the ``/data/big/signor/`` folder. 
- Take note of the dates provided in the downloaded file names and set the ``signor_date`` variable accordingly.
- Run the ``demo_02_parse_signor_kinase_substrates.m`` script on Matlab. 
