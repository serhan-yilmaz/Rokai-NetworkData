species_list = {'human', 'mouse', 'rat'};
for iSpecies = 1:length(species_list)
    species = species_list{iSpecies};

    fprintf('[Running] Preparing the files for transfer to R - %s\n', species)
    fname = ['rokai_network_data_uniprotkb_', species];
    load(['data/', fname, '.mat']);
    
    Wkin2site = NetworkData.Wkin2site;
    Wkin2kin = NetworkData.Wkin2kin;
    Wkin2kin_phospha = NetworkData.Wkin2kin_phospha;
    Wsite2site_coev = NetworkData.Wsite2site_coev;
    Wsite2site_sd = NetworkData.Wsite2site_sd;
    Wkin2site_psp = logical(NetworkData.KS.Wkin2site_psp);
    Wkin2site_psp_base = logical(NetworkData.KS.Wkin2site_psp_base);
    Wkin2site_signor = NetworkData.KS.Wkin2site_signor;
    Wphospha2site = NetworkData.Wphospha2site;
    
    version_psp = NetworkData.Versions.version_psp;
    version_signor = NetworkData.Versions.version_signor;
    version_string = NetworkData.Versions.version_string;
    version_ptmcode = NetworkData.Versions.version_ptmcode;
    version_depod = NetworkData.Versions.version_depod;
    version_uniprot = NetworkData.Versions.version_uniprot;
    
    folder = ['data/r/',fname, '/'];
    if(~exist(folder, 'dir')); mkdir(folder); end
    save([folder, 'rokai_networks_r.mat'], ...
        'Wkin2site', 'Wkin2kin', 'Wsite2site_coev', 'Wsite2site_sd', ...
        'Wkin2site_psp', 'Wkin2site_psp_base', 'Wkin2site_signor', ...
        'Wphospha2site', 'Wkin2kin_phospha', ...
        'version_psp', 'version_signor', 'version_string', 'version_ptmcode', ...
        'version_depod', 'version_uniprot');
    writetable(NetworkData.Site, [folder, 'site.csv']);
    writetable(NetworkData.Kinase, [folder, 'kinase.csv']);
    writetable(NetworkData.Gene, [folder, 'gene.csv']);
    writetable(NetworkData.Phosphatase, [folder, 'phosphatase.csv']);
    writetable(NetworkData.UniprotGene, [folder, 'uniprot_gene.csv']);
end
fprintf('[Done] Preparing the files for transfer to R\n');


