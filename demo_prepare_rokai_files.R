library(R.matlab)
#src = "D:/Archive/Matlab/Projects/Rokai_NetworkData/rokai_to_shiny/"
src = "C:/Users/Admin/Documents/MATLAB/Projects/Rokai_NetworkData/rokai_uniprot_allspecies/"

species_list = c('human', 'mouse', 'rat');

for (iSpecies in 1:length(species_list)){
  species = species_list[iSpecies]
  fname = paste('rokai_network_data_uniprotkb', species, sep = "_");
  folder = paste(src, "data/r/", fname, '/', sep = "");
  outfile = paste(src, "data/r/", fname, '.rds', sep = "");
  
  sites <- read.csv(paste(folder, "site.csv", sep=""))
  kinases <- read.csv(paste(folder, "kinase.csv", sep=""))
  phosphatases <- read.csv(paste(folder, "phosphatase.csv", sep=""))
  genes <- read.csv(paste(folder, "gene.csv", sep=""))
  uniprot.genes <- read.csv(paste(folder, "uniprot_gene.csv", sep=""))
  net <- readMat(paste(folder, "rokai_networks_r.mat", sep=""))
  
  NetworkData = list("Site" = sites, "Kinase" = kinases, "Gene" = genes,  
                     "Phosphatase" = phosphatases, "UniprotGene" = uniprot.genes, 
                     "Wkin2site" = net$Wkin2site, "net" = net)
  saveRDS(NetworkData, outfile)
}
#readRDS(outfile)
