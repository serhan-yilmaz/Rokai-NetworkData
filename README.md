# Rokai-NetworkData
This repository explains how to update the network data used in RokaiApp.

### Step 1 - PhosphoSitePlus Data
- Download the ``Kinase_Substrate_Dataset.gz`` and ``Phosphorylation_site_dataset.gz`` files from [PhosphoSitePlus website](https://www.phosphosite.org/staticDownloads) and extract their contents to the ``/data/input/phosphositeplus/`` folder. 
- Take note of the last modified dates of the extracted files and set the ``psp_version`` variable accordingly.
- Run the ``/src/demo_01_parse_psp_kinase_substrates.m`` script on Matlab. 

### Step 2 - Signor Data
- Visit the [Signor website](https://signor.uniroma2.it/downloads.php), select ``file xls (slower)`` option and download the datasets for ``human/mouse/rat`` species. Place the contents into the ``/data/input/signor/`` folder. 
- Take note of the dates provided in the downloaded file names and set the ``signor_date`` and ``signor_version`` variables accordingly.
- Run the ``/src/demo_02_parse_signor_kinase_substrates.m`` script on Matlab. 

### Step 3 - Uniprot ID Mapping Data
- Visit the [Uniprot data download website](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/) and download the ``HUMAN_9606_idmapping.dat.gz``, ``MOUSE_10090_idmapping.dat.gz``, ``RAT_10116_idmapping.dat.gz`` files. Extract their contents to the ``/data/input/uniprot/`` folder. 
- Take note of the last modified dates of the extracted files and set the ``uniprot_version`` variable accordingly.
- Run the ``/src/demo_03_uniprot_id_mapping.m`` script on Matlab. 

### Step 4 - Depod Data
- Visit [Depod download page](http://depod.bioss.uni-freiburg.de/download.php) and download the file under ``Interactions between human phosphatases and their protein substrates``. Place the contents into the ``/data/input/depod/`` folder. 
- Take note of the date provided in the downloaded file name and set the ``depod_date`` and ``depod_version`` variables accordingly.
- Run the ``/src/demo_04_parse_depod_phosphatase_substrates.m`` script on Matlab. 

### Step 5 - STRING Data
- Visit the [STRING download page](https://string-db.org/cgi/download), select ``Homo sapiens``, ``Mus musculus`` and ``Rattus norvegicus`` species, and download ``9606.protein.links.v11.5.txt.gz``, ``10090.protein.links.v11.5.txt.gz``, ``10116.protein.links.v11.5.txt.gz`` files. Extract their contents to the ``/data/input/string/`` folder. 
- Take note of the last modified dates and file names of the extracted files and set the ``string_version_number`` and ``string_version_date`` variables accordingly.
- Run the ``/src/demo_05_parse_string_ppi_network.m`` script on Matlab. 

### Step 6 - PTMcode Data
- Visit the [PTMcode download page](https://ptmcode.embl.de/data.cgi) and download ``PTMcode2_associations_within_proteins.txt.gz``, ``PTMcode2_associations_between_proteins.txt.gz`` files.  Extract their contents to the ``/data/input/ptmcode/`` folder. 
- Take note of the last modified dates of the extracted files and set the ``ptmcode_version`` variable accordingly.
- Run the ``/src/demo_06_parse_ptmcode_networks.m`` script on Matlab. 
- This step can be skipped entirely since PTMcode seems to stop the updates (last update date in 2014). 

### Step 7 - Preparing the NetworkData
- Run the ``demo_07_prepare_combined_network_data.m`` script on Matlab to combine all datasets into a single file. 
- This creates ``rokai_network_data_uniprotkb_human.mat``, ``rokai_network_data_uniprotkb_mouse.mat``, ``rokai_network_data_uniprotkb_rat.mat`` files in the ``/data/`` folder. 

### Step 8 - Transferring the files to R
- First, run the ``demo_transfer_to_R.m`` script on Matlab. This will prepare the files for transfer to R in the ``data/r/`` folder. 
- Next, open the ``demo_prepare_rokai_files.R`` script on R and set the ``src`` variable to the working directory (e.g., ``src = "C:/Users/Admin/Documents/MATLAB/Projects/Rokai_NetworkData/"``). 
- Run the ``demo_prepare_rokai_files.R`` script in R. 
- This creates ``rokai_network_data_uniprotkb_human.rds``, ``rokai_network_data_uniprotkb_mouse.rds``, ``rokai_network_data_uniprotkb_rat.rds`` files in the ``/data/r/`` folder. 



