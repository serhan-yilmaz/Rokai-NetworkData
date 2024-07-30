warning('off', 'MATLAB:table:ModifiedVarnames');
%%
string_version_number = 'v12.0';
string_version_date = '2023-05-15';
species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};
    fprintf('[Running] Parsing STRING PPI Network - %s\n', species)

    switch(species)
        case 'human'
            species_code = '9606.';
        case 'mouse'
            species_code = '10090.';
        case 'rat'
            species_code = '10116.';
        otherwise
            error('Unknown species code.');
    end
    
    dataFolder = '../data/input/string/';
    filename = [dataFolder, species_code, 'protein.links.', string_version_number, '.txt'];
    ds = datastore(filename);
    STRING_PPI = ds.readall();
    
    Proteins = unique([STRING_PPI.protein1; STRING_PPI.protein2]);
    [~, proteinIndices1] = ismember(STRING_PPI.protein1, Proteins);
    [~, proteinIndices2] = ismember(STRING_PPI.protein2, Proteins);
    Proteins = regexprep(Proteins, species_code, '');
    
    PPI = sparse(proteinIndices1, proteinIndices2, ...
        STRING_PPI.combined_score, length(Proteins), length(Proteins));
    clear STRING_PPI proteinIndices1 proteinIndices2
    
    % Save the results
    fprintf('[Running] Saving the output files - %s\n', species)
    outputFolder =  '../data/processed/string/';
    if(~exist(outputFolder, 'dir')); mkdir(outputFolder); end
    save([outputFolder, 'string_ppi_network_', species, '.mat'], 'Proteins', 'PPI', 'string_version_number', 'string_version_date');
end
fprintf('[Done] Parsing STRING PPI Network\n');



