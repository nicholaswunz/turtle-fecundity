# Dataset for turtle reproductive output
## data/reproductive output.csv  
### Variables descriptions
*author.name:* Surname of the first author of the paper;  
*paper.ID:* Unique identifiers for each paper;  
*species:* Binomial nomenclature for species (as *Genus_latin*);  
*red.list:* Endangered status classification defined by the International Union for Conservation of Nature [IUCN](https://www.iucnredlist.org/);  
*published.year:* Year of publication;  
*title:* Title of the paper extracted;  
*year:* Year of the study conducted;  
*locality:* General location of the study site;  
*country:* Country of the study;  
*lat:* Latitude decimal geographic coordinate in decimals (north–south position);  
*lon:* Longitude decimal georeferenced coordinate in decimals (east–west position);  
*MPA:* Logical binary (YES/NO) whether the coordinates overlay marine protected areas;  
*MPA2:* Logical binary (YES/NO) whether the coordinates overlay marine protected areas +5 lm buffer zone;  
*air.temp:* Mean air temperature of the study site during the nesting period (°C);  
*sand.temp:* Mean nest temperature of the study site during the nesting period (°C);  
*egg.diameter:* Mean diameter of the egg in millimetres (mm);  
*egg.volume:* Estimated egg volume in cm³;  
*egg.mass:* Mean egg mass in grams (g);  
*hatchling.length:* Mean hatchling length in millimetre (mm);  
*hatchling.mass:* Mean hatchling mass in grams (g);  
*clutch.size:* Mean number of eggs laid per female;  
*clutch.mass:* Mean clutch mass laid per female in grams (g);  
*clutch.freq:* Mean number of clutches laid during the breeding season;  
*pred.freq:* Estimated clutch frequency based on linear relationship of species-specific clutch size and female CCL;  
*tot.clutch:* Total clutch laid per female during the breeding season;  
*pred.clutch:* Estimated total clutch size (clutch.size x pred.freq);  
*length:* Mean curved carapace length in centimetre (cm);  
*mass:* Mean body mass in kilogram (kg);  
*pred.mass:* Estimated body mass (kg) based on the Mass-CCL relationship in the Supplementary file;  
*notes:* General comments;  
*og.author:* Review papers only – Surname of the first author of the paper extracted from the review paper;  
*pred.mass:* Review papers only – Title of the paper extracted from the review paper.


# Dataset for turtle egg energy content  
## data/egg energy.csv
### Variables descriptions  
*first.author:* Surname of the first author of the paper;  
*published.year:* Year of publication;  
*paper.ID:* Unique study identifier;  
*title:* Title of the paper extracted;  
*locality:* General location of the study site;  
*country:* Country of the study;  
*species:* Binomial nomenclature for species (as *Genus_latin*);  
*egg.energy:* Energy content of the egg analysed in joules (J);  
*egg.mass:* Wet mass in grams (g) of the egg analysed;  
*egg.mass:* Dry mass in grams (g) of the egg analysed;  
*energy.g:* energy content per gram of tissue (J g);  
*notes:* General comments;  
*link:* Link to the paper extracted.  


# Dataset for Terengganu population survey 
## data/nesting year.csv  
### Variables descriptions  
*year:* Year of the nesting survey;  
*ID:* Unique Identifier assigned to each nesting turtle;  
*recruit:* Logical binary (First time/Recurring) of whether the turtle surveyed is a first time nester on the island, or a recurring nester from the previous years;  
*Count of Year:* Number of a times the female arrived on the island (note not all females arriving on the beach did lay eggs);  
*sum_n of eggs:* Total number of eggs laid during the nesting season;  
*mean_n of eggs	:* Mean number of eggs laid per nesting event (clutch size);  
*SD_ n of eggs:* Variation in clutch size as standard deviation;  
*mean_length:* Mean curved carapace length in centimetre (cm);  
*SD_length:* Variation in curved carapace length as standard deviation;  
*mean_width:* Mean curved carapace width in centimetre (cm);  
*SD_width:* Variation in curved carapace width as standard deviation;  
*mass:* Estimated body mass in kilogram (kg) from the length-mass conversion in Ganyai (2017). 


# Dataset for global protected beaches by country
Dataset for the global percentage of protected beaches extracted from [Mazaris et al. 2014](sciencedirect.com/science/article/pii/S000632071400113X?casa_token=VW7hq97bFSIAAAAA:0gV4jzKB9krEjDhuJTkdlpP2NEAWxjb9YgxuK4MB0pxfTscTD0hUX90WlVCMMFbflTFiGSE)
## data/world turtle.csv
### Variables descriptions  
*name:* Country name;  
*protected site:* Percentage of nesting sites protected at species level;  
*temp:* Mean nest temperature by each country (°C);  
*CC%:* Percentage of *Caretta caretta* nesting sites protected;  
*CM%:* Percentage of *Chelonia mydas* nesting sites protected;  
*DC%:* Percentage of *Dermochelys coriacea* nesting sites protected;  
*EI%:* Percentage of *Eretmochelys imbricata* nesting sites protected;  
*LK%:* Percentage of *Lepidochelys kempii* nesting sites protected;  
*LO%:* Percentage of *Lepidochelys olivacea* nesting sites protected;  
*ND%:* Percentage of *Natator depressus* nesting sites protectedy;  
*CC_mass:* Mean *Caretta caretta* body mass in kg;  
*CM_mass:* Mean *Chelonia mydas* body mass in kg;  
*DC_mass:* Mean *Dermochelys coriacea* body mass in kg;  
*EI_mass:* Mean *Eretmochelys imbricata* body mass in kg;  
*LK_mass:* Mean *Lepidochelys kempii* body mass in kg;  
*LO_mass:* Mean *Lepidochelys olivacea* body mass in kg;  
*ND_mass:* Mean *Natator depressus* body mass in kg. 


# Turtle phylogeny data 
Phylogenetic matrix of all seven sea tuertle species used for the analysis `data/phylo_cor_all.Rdata`.


    
  


