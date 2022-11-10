warning('off', 'MATLAB:table:ModifiedVarnames');
%%
depod_date = '201903';
depod_version = '2019-03-01';
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};
    fprintf('[Running] Parsing Depod Phosphatase Substrates - %s\n', species)
    
    dataFolder = '../data/input/depod/';
    filename = [dataFolder, 'PPase_protSubtrates_', depod_date, '.xls'];
    ds = datastore(filename);
    ds.VariableTypes = repmat({'char'}, 1, length(ds.VariableNames));
    PTtable = ds.readall();
    PTtable.Flanking5Org = PTtable{:, 6};
    
    validRows = ~ismember(lower(PTtable.Flanking5Org), {'', 'n/a', 'na', 'null', 'multiple'});
    PTtable = PTtable(validRows, :);
    PTtable.Flanking5Org = regexprep(PTtable.Flanking5Org, '-', '_');
    
    delims = {',', ';', ' '};
    seqs = cellfun(@(x) split(x, delims), PTtable.Flanking5Org, 'UniformOutput', false);
    
    Aflanking = DynamicArray(2000, 'cell');
    Arows = DynamicArray(2000, 'numeric');
    for iRow = 1:length(seqs)
        seq_list = seqs{iRow};
        Aflanking.insert(seq_list);
        Arows.insert(repmat(iRow, length(seq_list), 1));
    end
    Aflanking = Aflanking.finalize();
    Arows = Arows.finalize();
    
    PTtable2 = PTtable(Arows, :);
    PTtable2.Flanking5 = Aflanking;
    PTtable2.FlankLength = cellfun(@length, PTtable2.Flanking5);
    validRows = PTtable2.FlankLength == 11;
    irregularRows = ~validRows & (PTtable2.FlankLength > 0);
    PTtable2irregular = PTtable2(irregularRows, :);
    PTtable2 = PTtable2(validRows, :);
    
    dataFolder = '../data/processed/signor/';
    load([dataFolder, 'psp_signor_sites_and_ks_network_', species, '.mat']);
    Site.Flanking5 = cellfun(@(x) x(3:(end-2)), Site.Flanking, 'UniformOutput', false);
    [b, ib] = ismember(lower(PTtable2.Flanking5), lower(Site.Flanking5));
    PTtable2.PSP_SiteIndex = ib;
    fprintf('[Info] %s - %d/%d DEPOD interactions are mapped\n', species, nnz(b), length(b));
    
    PTtable3 = PTtable2(b, :);
    PTtable3.PSP_SiteGene = Site.Gene(ib(b));
    
    uFlank5 = unique(lower(Site.Flanking5));
    [b, ib] = ismember(lower(PTtable3.Flanking5), uFlank5);
    Mptsite2usite = sparse(1:height(PTtable3), ib, true, height(PTtable3), length(uFlank5));
    
    [b, ib] = ismember(lower(Site.Flanking5), uFlank5);
    Mpspsite2usite = sparse(1:height(Site), ib, true, height(Site), length(uFlank5));
    
    Mptsite2pspsite = logical(double(Mptsite2usite) * Mpspsite2usite');
    [i1, i2] = find(Mptsite2pspsite);
    
    PTtable4 = PTtable3(i1, :);
    PTtable4.PSP_SiteIndex = i2;
    PTtable4.PSP_SiteGene = Site.Gene(PTtable4.PSP_SiteIndex);
    PTtable4.PSP_SiteProtein = Site.Protein(PTtable4.PSP_SiteIndex);
    PTtable4.PSP_SiteResidue = Site.Residue(PTtable4.PSP_SiteIndex);
    
    % Map the phosphothases
    dataFolder = '../data/processed/uniprot/';
    filename = [dataFolder, 'uniprot_id_mapping_', species, '_Gene_Name_Combined.csv'];
    ds = datastore(filename, 'ReadVariableNames', true);
    Mapping = ds.readall();
    
    [b, ib] = ismember(lower(PTtable4.PhosphataseEntryNames), lower(Mapping.Gene_Name));
    PTtable5 = PTtable4(b, :);
    PTtable5.PhosphaID = Mapping.Entry(ib(b));
    PTtable5.PhosphaGene = Mapping.Gene_Name(ib(b));
    
    [phosphathases, ib, ic] = unique(PTtable5.PhosphaGene);
    PTtable5.PhosphaIndex = ic;
    
    Phosphatase = table();
    Phosphatase.ID = PTtable5.PhosphaID(ib);
    Phosphatase.Gene = PTtable5.PhosphaGene(ib);
    
    Mphospha2pspsite = sparse(PTtable5.PhosphaIndex, PTtable5.PSP_SiteIndex, true, ...
            height(Phosphatase), height(Site));
    
    Stats = struct();
    Stats.nSite = nnz(sum(Mphospha2pspsite, 1) > 0);
    Stats.nPhosphatase = nnz(sum(Mphospha2pspsite, 2) > 0);
    Stats.nnzTotal = nnz(Mphospha2pspsite);
%     Stats
    
    fprintf('[Running] Saving the output files - %s\n', species)
    outputFolder =  '../data/processed/depod/';
    if(~exist(outputFolder, 'dir')); mkdir(outputFolder); end
    save([outputFolder, 'depod_phosphatase_targets_', species, '.mat'], ...
        'Phosphatase', 'Mphospha2pspsite', 'Stats', 'depod_date', 'depod_version');
end
fprintf('[Done] Parsing Depod Phosphatase Substrates\n')


