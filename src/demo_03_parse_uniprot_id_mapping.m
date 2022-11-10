warning('off', 'MATLAB:table:ModifiedVarnames');
%%
uniprot_version = '2022-10-04';
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};

    fprintf('[Running] Parsing Uniprot ID Mapping - %s\n', species)

    switch(species)
        case 'human'
            species_uniprot = 'HUMAN_9606';
        case 'mouse'
            species_uniprot = 'MOUSE_10090';
        case 'rat'
            species_uniprot = 'RAT_10116';
        otherwise
            species_uniprot = ['_', species];
    end
    
    % Needed Mappings: 
    % - STRING.Proteins -> Ensembl_PRO
    % - PTMcode.Gene -> Gene_Name
    % - Depod.Gene -> Gene_Name
    
    % Relevant Fields:
    % - Ensembl_PRO
    % - Gene_Name
    % - UniProtKB-ID
    
    dataFolder = '../data/input/uniprot/';
    filename = [dataFolder, species_uniprot, '_idmapping.dat'];
    ds = datastore(filename, 'ReadVariableNames', false);
    %ds.VariableTypes = repmat({'char'}, 1, length(ds.VariableNames));
    UniprotTable = ds.readall();
    
    [types, ~, ic] = unique(UniprotTable.Var2);
    UniprotTable.TypeIndex = ic;
    
    [~, si] = sort(UniprotTable.TypeIndex, 'ascend');
    UniprotTable = UniprotTable(si, :);
    
    outputFolder = '../data/processed/uniprot/';
    if(~exist(outputFolder, 'dir')); mkdir(outputFolder); end
    
    [~, indices_last] = unique(UniprotTable.TypeIndex, 'last');
    indices_last = [0; indices_last];
    
    [~, selectedTypes] = ismember({'Ensembl_PRO', 'Gene_Name', 'Gene_Synonym', 'UniProtKB-ID'}, types);
    % selectedTypes = 1:length(types); % For all types
    

    S = struct();
    for iType = selectedTypes
        type = types{iType};
        first_index = indices_last(iType) + 1;
        last_index = indices_last(iType+1);
        ind = first_index:last_index;
        T = UniprotTable(ind, :);
        
        a = cellfun(@(x) split(x, '-'), T.Var1, 'UniformOutput', false);
        b = cellfun(@(y) y{1}, a, 'UniformOutput', false);

        Tx = table();
        Tx.Entry = b;
        Tx.EntryWithIsoform = T.Var1;
        Tx.(type) = T.Var3;
        if(strcmpi(type, 'Ensembl_PRO'))
            a = cellfun(@(x) split(x, '.'), T.Var3, 'UniformOutput', false);
            b = cellfun(@(y) y{1}, a, 'UniformOutput', false);
            Tx.Ensembl_Protein = b;
        end
        switch(type)
            case 'Gene_Name'
                S.TGeneName = Tx;
            case 'Gene_Synonym'
                S.TGeneSynonym = Tx;
        end
        fprintf('[Running] Saving the output - %s - %s\n', species, type)
        writetable(Tx, [outputFolder, 'uniprot_id_mapping_', species, '_', type, '.csv']);
    end
    T1 = table();
    T1.Entry = S.TGeneName.Entry;
    T1.EntryWithIsoform = S.TGeneName.EntryWithIsoform;
    T1.Gene_Name = S.TGeneName.Gene_Name;
    T2 = table();
    T2.Entry = S.TGeneSynonym.Entry;
    T2.EntryWithIsoform = S.TGeneSynonym.EntryWithIsoform;
    T2.Gene_Name = S.TGeneSynonym.Gene_Synonym;
    Tx = [T1; T2];
    type = 'Gene_Name_Combined';
    fprintf('[Running] Saving the output - %s - %s\n', species, type)
    writetable(Tx, [outputFolder, 'uniprot_id_mapping_', species, '_', type, '.csv']);

    save([outputFolder, 'uniprot_version_', species, '.mat'], 'uniprot_version'); 
end
fprintf('[Done] Uniprot ID Mapping\n')




