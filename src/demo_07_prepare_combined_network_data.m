warning('off', 'MATLAB:table:ModifiedVarnames');
%%
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};
    fprintf('[Running] Preparing NetworkData - %s\n', species)
    % Load networks
    dataFolder = '../data/processed/';
    PSP = load([dataFolder, 'signor/psp_signor_sites_and_ks_network_', species, '.mat']);
    
    % Protein-Protein Interaction (PPI) network obtained from STRING.
    % See 'load_string_ppi_network.m' for more information. 
    STRING = load([dataFolder, 'string/string_ppi_network_', species, '.mat']);
    
    % Phosphosite-level interaction network based on co-evolution and/or 
    % structure distance evidence obtained from PTMcode. 
    % See 'load_ptmcode_networks.m' for more information. 
    PTMcode = load([dataFolder, 'ptmcode/ptmcode_networks_', species, '.mat']);
    
    % Phosphatase-Substrate network from Depod. 
    % See 'load_depod_phosphothase_substrates.m' for more information. 
    DEPOD = load([dataFolder, 'depod/depod_phosphatase_targets_', species, '.mat']);
    
    % ------------------------------------------------------------------------
    % --- Mapping Ensembl protein identifiers to UniprotKb ---
    % In order to create the kinase-kinase interaction network, we map the
    % string protein identifiers (ENSP) to psp kinase identifiers (UniprotKB).
    % We obtain the ensembl protein to uniprotkb id mappings using Uniprot
    % mapping tool: https://www.uniprot.org/uploadlists/
    % The 'string_proteins_uniprotkb_mouse.tab' contains the results of the query
    % run on 2021-02-03, using the ids given in 'STRING.Proteins' as input. 
    
    dataFolder = '../data/processed/uniprot/';
    filename = [dataFolder, 'uniprot_id_mapping_', species, '_Ensembl_PRO.csv'];
    ds = datastore(filename, 'ReadVariableNames', true);
    Mapping = ds.readall();
    
    % Filter for kinase identifiers
    validRows = ismember(Mapping.Entry, PSP.Kinase.ID);
    M = Mapping(validRows, :);
    
    % Map to primary identifier
    M.PrimaryEnsp = regexprep(M.Ensembl_Protein, ',[^.]*', '');
    
    validRows = ismember(M.PrimaryEnsp, STRING.Proteins);
    M = M(validRows, :);
    
    [~, kin_indices] = ismember(M.Entry, PSP.Kinase.ID);
    [~, ppi_protein_indices] = ismember(M.PrimaryEnsp, STRING.Proteins);
    Mkin2ppi = sparse(kin_indices, ppi_protein_indices, ...
            true, height(PSP.Kinase), length(STRING.Proteins));
    
    % Create kinase-kinase interaction network, a subset of STRING ppi network
    % Edges are weighted according to the combined score provided by STRING.
    Kinase_PPI = Mkin2ppi * STRING.PPI * Mkin2ppi';
    
    % - Kinase+Phosphatase PPI - 
    % Filter for kinase+phosphatase identifiers
    kinIds = [PSP.Kinase.ID; DEPOD.Phosphatase.ID];
    
    validRows = ismember(Mapping.Entry, kinIds);
    M = Mapping(validRows, :);
    
    % Map to primary identifier
    M.PrimaryEnsp = regexprep(M.Ensembl_Protein, ',[^.]*', '');
    
    validRows = ismember(M.PrimaryEnsp, STRING.Proteins);
    M = M(validRows, :);
    
    [~, kin_indices] = ismember(M.Entry, kinIds);
    [~, ppi_protein_indices] = ismember(M.PrimaryEnsp, STRING.Proteins);
    Mkin2ppi_combined = sparse(kin_indices, ppi_protein_indices, ...
            true, length(kinIds), length(STRING.Proteins));
    
    % Create kinase-kinase interaction network, a subset of STRING ppi network
    % Edges are weighted according to the combined score provided by STRING.
    Kinase_PPI_combined = Mkin2ppi_combined * STRING.PPI * Mkin2ppi_combined';
    
    % ------------------------------------------------------------------------
    % --- Mapping PTMcode Gene Identifiers to UniprotKb ---
    % In order to create the phosphosite level interaction network, we map the
    % PTMcode gene identifiers to psp protein identifiers (UniprotKB).
    % We obtain the gene name and ensg to uniprotkb id mappings using 
    % Uniprot mapping tool: https://www.uniprot.org/uploadlists/
    % The 'ptmcode_genes_uniprotkb.tab' contains the results of the query
    % run on 2021-02-03, using the ids given in 'PTMcode.Gene' as input. 
    
    dataFolder = '../data/processed/uniprot/';
    filename = [dataFolder, 'uniprot_id_mapping_', species, '_Gene_Name_Combined.csv'];
    ds = datastore(filename, 'ReadVariableNames', true);
    Mapping = ds.readall();
    Mapping.Gene = Mapping.Gene_Name;
    
    % dataFolder = '../data/raw__/';
    % filePath = [dataFolder, 'ptmcode_genes_uniprotkb_', species, '.tab'];
    % ds = tabularTextDatastore(filePath, 'FileExtensions', '.tab', 'Delimiter', '\t');
    % ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
    % Mapping = ds.readall();
    
    % Map to primary identifiers 
    Mapping.Gene = regexprep(Mapping.Gene, ',[^.]*', '');
    
    [b, indices] = ismember(Mapping.Entry, PSP.Protein);
    [b2, indices2] = ismember(lower(Mapping.Gene), lower(regexprep(PTMcode.Gene, ' ', '')));
    M = Mapping(b&b2, :);
    Mmapping2pspprotein = sparse((1:height(M))', indices(b&b2), ...
            1, height(M), length(PSP.Protein));
    Mmapping2ptmcodegene = sparse((1:height(M))', indices2(b&b2), ...
            1, height(M), length(PTMcode.Gene));
    
    % Mapping of PSP proteins to PTM Gene Identifiers
    Mpspprotein2ptmcodegene = Mmapping2pspprotein' * Mmapping2ptmcodegene;
    
    % Map the ptmcode sites to psp sites according to protein/gene ids
    [~, pspProteinIndices] = ismember(PSP.Site.Protein, PSP.Protein);
    Mpsp2pspprotein = sparse((1:height(PSP.Site))', ...
        pspProteinIndices, 1, height(PSP.Site), length(PSP.Protein));
    
    Mptmcode2ptmcodegene = sparse((1:height(PTMcode.Site))', ...
        PTMcode.Site.GeneIndex, 1, height(PTMcode.Site), length(PTMcode.Gene));
    
    % Mapping of PSP sites to PTMcode sites according to protein/gene ids
    Mprotein_psp2ptmcode = Mpsp2pspprotein ...
                         * Mpspprotein2ptmcodegene ...
                         * Mptmcode2ptmcodegene';
    
    % Map the ptmcode sites to psp sites according to residue
    residues = [PSP.Site.Residue; PTMcode.Site.Residue];
    [~, pspResidueIndices] = ismember(PSP.Site.Residue, residues);
    [~, ptmcodeResidueIndices] = ismember(PTMcode.Site.Residue, residues);
    
    Mpsp2residue = sparse((1:height(PSP.Site))', ...
        pspResidueIndices, 1, height(PSP.Site), length(residues));
    Mptmcode2residue = sparse((1:height(PTMcode.Site))', ...
        ptmcodeResidueIndices, 1, height(PTMcode.Site), length(residues));
    
    % Mapping of P  SP sites to PTMcode sites according to residue
    Mresidue_psp2ptmcode = Mpsp2residue * Mptmcode2residue';
    
    % Combined mapping of PSP sites to PTMcode sites
    Mpsp2ptmcode = double(Mprotein_psp2ptmcode & Mresidue_psp2ptmcode);
    
    % Map the co-evolution network to PSP sites
    Wcoev = logical(Mpsp2ptmcode * PTMcode.Wcoev * Mpsp2ptmcode');
    
    % Map the structure distance network to PSP sites
    Wsd = logical(Mpsp2ptmcode * PTMcode.Wsd * Mpsp2ptmcode');
    
    % Cleanup
    clear Mmapping2pspprotein Mmapping2ptmcodegene Mpspprotein2ptmcodegene
    clear Mpsp2pspprotein Mptmcode2ptmcodegene
    clear Mpsp2residue Mptmcode2residue
    clear Mprotein_psp2ptmcode Mresidue_psp2ptmcode
    
    % Prepare output 
    Kinase = table();
    Kinase.KinaseID = PSP.Kinase.ID;
    Kinase.KinaseName = PSP.Kinase.Name;
    Kinase.Gene = PSP.Kinase.Gene;
    
    position = regexprep(PSP.Site.Residue, '[^0-9]', '');
    
    Site = table();
    Site.Identifier = cellstr(join([PSP.Site.Protein, position], '_'));
    Site.Gene = PSP.Site.Gene;
    Site.Protein = PSP.Site.Protein;
    Site.Position = PSP.Site.Residue;
    Site.Flanking = PSP.Site.Flanking;
    
    [genes, ~, ic] = unique(Site.Gene);
    Gene = table();
    Gene.Name = genes;
    Wgene2site = sparse(ic, (1:height(Site))', true, height(Gene), height(Site));
    
    % Uniprot
    dataFolder = '../data/processed/uniprot/';
    filename = [dataFolder, 'uniprot_id_mapping_', species, '_Gene_Name.csv'];
    ds = datastore(filename, 'ReadVariableNames', true);
    UniprotProt2Gene = ds.readall();
    
    UniprotGene = table();
    UniprotGene.ID = UniprotProt2Gene.Entry;
    UniprotGene.Gene = UniprotProt2Gene.Gene_Name;

    dataFolder = '../data/processed/uniprot/';
    load([dataFolder, 'uniprot_version_', species, '.mat']);
    
    % 
    NetworkData = struct();
    NetworkData.Site = Site;
    NetworkData.Kinase = Kinase;
    NetworkData.Phosphatase = DEPOD.Phosphatase;
    NetworkData.Gene = Gene;
    NetworkData.UniprotGene = UniprotGene;
    NetworkData.Wkin2site = PSP.KS;
    NetworkData.Wkin2kin = Kinase_PPI;
    NetworkData.Wkin2kin_phospha = Kinase_PPI_combined;
    NetworkData.Wsite2site_coev = Wcoev;
    NetworkData.Wsite2site_sd = Wsd;
    NetworkData.Wgene2site = Wgene2site;
    NetworkData.Wphospha2site = DEPOD.Mphospha2pspsite;
    NetworkData.KS = struct();
    NetworkData.KS.Wkin2site_psp = PSP.KSpsp;
    NetworkData.KS.Wkin2site_psp_base = PSP.KSbase;
    NetworkData.KS.Wkin2site_signor = PSP.KSsignor;
    NetworkData.nSite = height(Site);
    NetworkData.nKinase = height(Kinase);
    NetworkData.nPhosphatase = DEPOD.Stats.nPhosphatase;
    NetworkData.nGene = height(Gene);
    NetworkData.Stats = PSP.Stats;
    NetworkData.Stats.nProtein = length(unique(Site.Protein));
    NetworkData.Stats.nSiteWithKinase = nnz(sum(NetworkData.Wkin2site, 1)>0);
    NetworkData.Stats.nProteinWithKinase = length(unique(NetworkData.Site.Protein(sum(NetworkData.Wkin2site, 1)>0)));
    NetworkData.Stats.nnzPPI = nnz(NetworkData.Wkin2kin);
    NetworkData.Stats.nnzPPI_phospha = nnz(NetworkData.Wkin2kin_phospha);
    NetworkData.Stats.nnzPTMcodeCoev = nnz(NetworkData.Wsite2site_coev);
    NetworkData.Stats.nnzPTMcodeSD = nnz(NetworkData.Wsite2site_sd);
    NetworkData.Stats.nPhosphatase = DEPOD.Stats.nPhosphatase;
    NetworkData.Stats.nSiteWithPhospha = DEPOD.Stats.nSite;
    NetworkData.Stats.nProteinWithPhospha = length(unique(NetworkData.Site.Protein(sum(NetworkData.Wphospha2site, 1)>0)));
    NetworkData.Stats.nnzDEPOD = DEPOD.Stats.nnzTotal;
    NetworkData.Versions = struct();
    NetworkData.Versions.version_uniprot = uniprot_version;
    NetworkData.Versions.version_psp = PSP.psp_version;
    NetworkData.Versions.version_signor = PSP.signor_version;
    NetworkData.Versions.version_ptmcode = PTMcode.ptmcode_version;
    NetworkData.Versions.version_string = STRING.string_version_date;
    NetworkData.Versions.version_depod = DEPOD.depod_version;
%     NetworkData.Versions.version_uniprot = '2022-10-04';
%     NetworkData.Versions.version_psp = '2021-04-19';
%     NetworkData.Versions.version_signor = '2021-05-21';
%     NetworkData.Versions.version_ptmcode = '2014-09-17';
%     NetworkData.Versions.version_string = '2018-12-20';
%     NetworkData.Versions.version_depod = '2019-03-01';
    
    % Save the results
    fprintf('[Running] Saving the output files - %s\n', species)
    outputFolder =  '../data/';
    save([outputFolder, 'rokai_network_data_uniprotkb_', species, '.mat'], 'NetworkData');
end

fprintf('[Done] Preparing NetworkData\n');







