warning('off', 'MATLAB:table:ModifiedVarnames');
%%
psp_version = '2022-10-22';
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};

    fprintf('[Running] Parsing PSP Kinase Substrates - %s\n', species)

    % - Load Kinase-Substrate network - 
    % We obtain the kinase-substrate annotations from PhosphositePlus at:
    % https://www.phosphosite.org/staticDownloads
    dataFolder = '../data/input/phosphositeplus/';
    filename = [dataFolder, 'Kinase_Substrate_Dataset'];
    ds = datastore(filename);
    ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
    KStable = ds.readall();
    % Filter for KIN_ORGANISM = 'mouse, human' and SUB_ORGANISM = 'mouse'
    % validRows = ismember(KStable.KIN_ORGANISM, {'mouse', 'human'});
    % validRows = validRows & strcmpi(KStable.SUB_ORGANISM, 'mouse');
    validRows = strcmpi(KStable.SUB_ORGANISM, species);
    KStable = KStable(validRows, :);
    
    % Remove isoform identifiers
    KStable.KIN_ACC_ID = regexprep(KStable.KIN_ACC_ID, '-[^.]', '');
    KStable.SUB_ACC_ID = regexprep(KStable.SUB_ACC_ID, '-[^.]', '');
    
    % Extract Kinases
    [kins, ib, ic] = unique(KStable.KIN_ACC_ID);
    KStable.KinaseIndex = ic;
    
    Kinase = table();
    Kinase.ID = kins;
    Kinase.Gene = KStable.GENE(ib);
    Kinase.Name = KStable.KINASE(ib);
    Kinase.Organism = KStable.KIN_ORGANISM(ib);
    
    % Create unique substrate identifier
    KStable.SubstrateID = cellstr(join([KStable.SUB_ACC_ID, KStable.SUB_MOD_RSD], '_'));
    
    % Extract Substrates
    [subst, ib, ic] = unique(KStable.SubstrateID);
    KStable.SubstrateIndex = ic;
    
    Substrate = table();
    Substrate.ID = subst;
    Substrate.Gene = KStable.SUB_GENE(ib);
    Substrate.Name = KStable.SUBSTRATE(ib);
    Substrate.Protein = KStable.SUB_ACC_ID(ib);
    Substrate.Residue = KStable.SUB_MOD_RSD(ib);
    Substrate.Flanking = KStable.SITE____7_AA(ib);
    
    % Create Kinase-Substrate network
    KS = logical(sparse(KStable.KinaseIndex, KStable.SubstrateIndex, ...
            1, height(Kinase), height(Substrate)));

    % - Load Known Phosphorylation Sites -
    % We obtain a list of known phosphosites from PhosphositePlus.
    % Due to its excessive file size, 'Phosphorylation_site_dataset'
    % is not included in this repository, but is available at: 
    % https://www.phosphosite.org/staticDownloads
    dataFolder = '../data/input/phosphositeplus/';
    filename = [dataFolder, 'Phosphorylation_site_dataset'];
    ds = datastore(filename);
    ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
    Sitetable = ds.readall();
    
    % Filter for ORGANISM = 'mouse'
    validRows = strcmpi(Sitetable.ORGANISM, species);
    Sitetable = Sitetable(validRows, :);
    
    % Create a unique identifier
    Sitetable.Residue = regexprep(Sitetable.MOD_RSD, '-[^.]*', '');
    Sitetable.ID = join([Sitetable.ACC_ID Sitetable.Residue], '_');
    
    % Combine the site table with substrate table
    siteIds = unique([Sitetable.ID; Substrate.ID]);
    Site = table();
    Site.Gene = cell(length(siteIds), 1);
    Site.Protein = cell(length(siteIds), 1);
    Site.Residue = cell(length(siteIds), 1);
    Site.Flanking = cell(length(siteIds), 1);
    
    % Map the sitetable to combined 'Site' table
    [b, indices] = ismember(siteIds, Sitetable.ID);
    Site.Gene(b) = Sitetable.GENE(indices(b));
    Site.Protein(b) = Sitetable.ACC_ID(indices(b));
    Site.Residue(b) = Sitetable.Residue(indices(b));
    Site.Flanking(b) = Sitetable.SITE____7_AA(indices(b));
    
    % Map the substrates to combined 'Site' table
    [b, indices] = ismember(siteIds, Substrate.ID);
    Site.Gene(b) = Substrate.Gene(indices(b));
    Site.Protein(b) = Substrate.Protein(indices(b));
    Site.Residue(b) = Substrate.Residue(indices(b));
    Site.Flanking(b) = Substrate.Flanking(indices(b));
    
    % Map the kinase-substate network to combined site table
    [~, indices] = ismember(Substrate.ID, siteIds);
    Msubstrate2site = sparse((1:height(Substrate))', ... 
        indices, 1, height(Substrate), height(Site));
    KS = logical(KS * Msubstrate2site);
    
    % Add kinase-substrate links from other organisms
    mouseKinases = strcmpi(Kinase.Organism, species);
    otherKinases = ~mouseKinases;
    
    KinaseOther = Kinase(otherKinases, :);
    KSother = KS(otherKinases, :);
    
    Kinase = Kinase(mouseKinases, :);
    KSbase = KS(mouseKinases, :);
    
    [b, ib] = ismember(KinaseOther.Name, Kinase.Name);
    KSother = KSother(b, :);
    Mkinother2kin = sparse((1:size(KSother, 1))', ib(b), ...
        true, size(KSother, 1), height(Kinase));
    
    KSmapped = logical(double(Mkinother2kin') * KSother);
    
    KS = KSbase | KSmapped;
    
    % Create a list of unique proteins
    Protein = unique(Site.Protein);
    
    % % Remove redundant fields
    clear KStable Sitetable
    
    fprintf('[Running] Saving the output files - %s\n', species)

    % Save the results
    outputFolder =  '../data/processed/phosphositeplus/';
    if(~exist(outputFolder, 'dir')); mkdir(outputFolder); end
    save([outputFolder, 'psp_sites_and_kinase_substrate_network_', species, '.mat'], ...
        'Protein', 'Site', 'Kinase', 'KS', 'KSbase', 'psp_version');
end
fprintf('[Done] Parsing PSP Kinase Substrates\n')







