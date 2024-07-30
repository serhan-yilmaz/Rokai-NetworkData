warning('off', 'MATLAB:table:ModifiedVarnames');
%%
signor_date = '30_07_24';
signor_version = '2024-07-30';
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};
    fprintf('[Running] Parsing Signor Kinase Substrates - %s\n', species)

    switch(species)
        case 'mouse'
            species_signor = '_M_musculus';
        case 'human'
            species_signor = '';
        case 'rat'
            species_signor = '_R_norvegicus';
        otherwise
            species_signor = ['_', species];
    end
    
    dataFolder = '../data/input/signor/';
    filename = [dataFolder, 'all_data', species_signor, '_', signor_date, '.xlsx'];
    ds = datastore(filename);
    ds.VariableTypes = repmat({'char'}, 1, length(ds.VariableNames));
    KStable = ds.readall();
    
    KStable = KStable(strcmpi(KStable.TYPEA, 'protein'), :);
    KStable = KStable(strcmpi(KStable.TYPEB, 'protein'), :);
    KStable = KStable(strcmpi(KStable.MECHANISM, 'phosphorylation'), :);
    KStable = KStable(~strcmpi(KStable.RESIDUE, ''), :);
    
    KStable.Residue = regexprep(KStable.RESIDUE, '[a-z]', '');
    KStable.PositionStr = regexprep(KStable.Residue, '[^0-9]', '');
    
    % Create unique substrate identifier
    KStable.SubstrateProtein = KStable.IDB;
    KStable.SubstrateID = cellstr(join([KStable.IDB, KStable.PositionStr], '_'));
    % KStable.SubstrateRID = cellstr(join([KStable.IDB, KStable.Residue], '_'));
    KStable.KinaseID = KStable.IDA;
    KStable.RowID = cellstr(join([KStable.KinaseID, KStable.SubstrateID], '_'));
    
    [~, ib] = unique(KStable.RowID);
    KStable = KStable(ib, :);
    
    dataFolder = '../data/processed/phosphositeplus/';
    load([dataFolder, 'psp_sites_and_kinase_substrate_network_', species, '.mat']);
    
    validKins = ismember(KStable.KinaseID, Kinase.ID);
    KStable = KStable(validKins, :);
    validProteins = ismember(KStable.SubstrateProtein, Site.Protein);
    KStable = KStable(validProteins, :);
    
    Site.PosStr = regexprep(Site.Residue, '[^0-9]', '');
    Site.Identifier = cellstr(join([Site.Protein, Site.PosStr], '_'));
    validSites = ismember(KStable.SubstrateID, Site.Identifier);
    KStable2 = KStable(validSites, :);
    KStable2b = KStable(~validSites, :);
    
    [~, ikins] = ismember(KStable2.KinaseID, Kinase.ID);
    [~, isites] = ismember(KStable2.SubstrateID, Site.Identifier);
    Wks = logical(sparse(ikins, isites, 1, height(Kinase), height(Site)));
    
    KSpsp = KS;
    KS = KSpsp | Wks;
    KSsignor_base = Wks;
    
    Site.PosStr = [];
    Site.Identifier = [];
    
    [proteins, ib] = unique(Site.Protein);
    genes = Site.Gene(ib);
    
    [~, ib] = ismember(KStable2b.SubstrateProtein, proteins);
    KStable2b.ProteinIndex = ib;
    KStable2b.Gene = genes(ib);
    
    [~, ib] = unique(KStable2b.SubstrateID);
    
    Site2 = table();
    Site2.Gene = KStable2b.Gene(ib);
    Site2.Protein = KStable2b.SubstrateProtein(ib);
    Site2.Residue = KStable2b.Residue(ib);
    Site2.Flanking = KStable2b.SEQUENCE(ib);
    % Site2.PosStr = 
    
    nSitePSP = height(Site);
    Site = [Site; Site2];
    nKin = size(KS, 1);
    KS = [KS sparse(nKin, height(Site2))];
    KSpsp = [KSpsp sparse(nKin, height(Site2))];
    KSbase = [KSbase sparse(nKin, height(Site2))];
    KSsignor_base = [KSsignor_base sparse(nKin, height(Site2))];
    
    Site.PosStr = regexprep(Site.Residue, '[^0-9]', '');
    Site.Identifier = cellstr(join([Site.Protein, Site.PosStr], '_'));
    [~, ikins] = ismember(KStable2b.KinaseID, Kinase.ID);
    [~, isites] = ismember(KStable2b.SubstrateID, Site.Identifier);
    KSaddition = logical(sparse(ikins, isites, 1, height(Kinase), height(Site)));
    
    Site.PosStr = [];
    Site.Identifier = [];
    
    KS1 = KS;
    KS = KS | KSaddition;
    
    Protein = unique(Site.Protein);
    KSsignor = KSaddition | KSsignor_base;
    
    Stats = struct();
    Stats.nSiteTotal = height(Site);
    Stats.nSitePSP = nSitePSP;
    Stats.nKinase = height(Kinase);
    Stats.nnzTotal = nnz(KS);
    Stats.nnzPSP = nnz(KSpsp);
    Stats.nnzPSPbase = nnz(KSbase);
    Stats.nnzSignor = nnz(KSsignor);
    Stats.nnzSignorBase = nnz(KSsignor_base);
    
    fprintf('[Running] Saving the output files - %s\n', species)
    outputFolder =  '../data/processed/signor/';
    if(~exist(outputFolder, 'dir')); mkdir(outputFolder); end
    save([outputFolder, 'psp_signor_sites_and_ks_network_', species, '.mat'], ...
        'Protein', 'Site', 'Kinase', 'KS', 'KSbase', 'KSpsp', 'KSsignor', ...
        'nSitePSP', 'Stats', 'signor_date', 'signor_version', 'psp_version');
end
fprintf('[Done] Parsing Signor Kinase Substrates\n')



