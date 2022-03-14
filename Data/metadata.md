# Dataset for turtle reproductive output
## data/meta-raw_data.csv  
### Variables descriptions
*author:* Surname of the first author of the paper;  
*study_ID:* Unique identifiers for each paper;  
*species:* Binomial nomenclature for species (as *Genus_latin*);  
*red_list:* Endangered status classification defined by the International Union for Conservation of Nature [IUCN](https://www.iucnredlist.org/);  
*sample_year:* Year of sampling for biometric parameters; 
*location:* General location of the study site;  
*pop_ID:* Unique identifier for populations;  
*country:* Country of the study;  
*lat:* Latitude decimal geographic coordinate in decimals (north–south position);  
*lon:* Longitude decimal georeferenced coordinate in decimals (east–west position);  
*MPA:* Logical binary (YES/NO) whether the coordinates overlay marine protected areas;  
*MPA2:* Logical binary (YES/NO) whether the coordinates overlay marine protected areas +5 km buffer zone;  
*air_temp:* Mean air temperature of the study site during the nesting period (°C);  
*sand_temp:* Mean nest temperature of the study site during the nesting period (°C);  
*egg_diameter_mm:* Mean diameter of the egg in millimetres (mm);  
*egg_volume_cm3:* Estimated egg volume in cm³;  
*egg_mass_g:* Mean egg mass in grams (g);  
*hatchling_length_mm:* Mean hatchling length in millimetre (mm);  
*hatchling_mass_g:* Mean hatchling mass in grams (g);  
*clutch_size:* Mean number of eggs laid per female;  
*clutch_mass_g:* Mean clutch mass laid per female in grams (g);  
*clutch_freq:* Mean number of clutches laid during the breeding season;  
*pred_freq:* Estimated clutch frequency based on linear relationship of species-specific clutch size and female CCL;  
*total_clutch:* Total clutch laid per female during the breeding season;  
*pred_clutch:* Estimated total clutch size (clutch_size x pred_freq);  
*adult_length_cm:* Mean adult curved carapace length in centimetre (cm);  
*adult_mass_kg:* Mean adult body mass in kilogram (kg);  
*pred_mass:* Estimated body mass (kg) based on the Mass-CCL relationship in the Supplementary file;  
*pub_year:* Year of publication; 
*notes:* General comments;  
*title:* Title of the paper extracted;  
*orig_author:* Review papers only – Surname of the first author of the paper extracted from the review paper;  
*orig_title:* Review papers only – Title of the paper extracted from the review paper;  
*link:* Link to paper extracted.


# Dataset for turtle egg energy content  
## data/egg energy.csv
### Variables descriptions  
*author:* Surname of the first author of the paper;  
*published_year:* Year of publication;  
*study_ID:* Unique study identifier;  
*title:* Title of the paper extracted;  
*locality:* General location of the study site;  
*country:* Country of the study;  
*species:* Binomial nomenclature for species (as *Genus_latin*);  
*energy_kJ:* Energy content of the egg analysed in kilojoules (kJ);  
*egg_mass_g:* Wet mass in grams (g) of the egg analysed;  
*dry_mass_g:* Dry mass in grams (g) of the egg analysed;  
*energy_KJ_g:* energy content per gram of tissue (J g);  
*notes:* General comments;  
*link:* Link to the paper extracted.  


# Dataset for the Redang Island population survey 
## data/nesting year.csv  
### Variables descriptions  
*year:* Year of the nesting survey;  
*ID:* Unique Identifier assigned to each nesting turtle;  
*recruit:* Logical binary (First time/Recurring) of whether the turtle surveyed is a first time nester on the island, or a recurring nester from the previous years;  
*count_year:* Number of a times the female arrived on the island (note not all females arriving on the beach did lay eggs);  
*sum_eggs:* Total number of eggs laid during the nesting season;  
*mean_eggs:* Mean number of eggs laid per nesting event (clutch size);  
*SD_eggs:* Variation in clutch size as standard deviation;  
*mean_length_cm:* Mean curved carapace length in centimetre (cm);  
*SD_length:* Variation in curved carapace length as standard deviation;  
*mean_width_cm:* Mean curved carapace width in centimetre (cm);  
*SD_width:* Variation in curved carapace width as standard deviation;  
*mass_kg:* Estimated body mass in kilogram (kg) from the length-mass conversion in Ganyai (2017). 


# Dataset for global protected beaches by country
Dataset for the global percentage of protected beaches extracted from Mazaris et al. 2014
## data/world_turtle.csv
### Variables descriptions  
*name:* Country name;  
*all_protected_perc:* Percentage of nesting sites protected at species level;  
*mean_temp:* Mean nest temperature by each country (°C);  
*CC_perc:* Percentage of *Caretta caretta* nesting sites protected;  
*CM_perc:* Percentage of *Chelonia mydas* nesting sites protected;  
*DC_perc:* Percentage of *Dermochelys coriacea* nesting sites protected;  
*EI_perc:* Percentage of *Eretmochelys imbricata* nesting sites protected;  
*LK_perc:* Percentage of *Lepidochelys kempii* nesting sites protected;  
*LO_perc:* Percentage of *Lepidochelys olivacea* nesting sites protected;  
*ND_perc:* Percentage of *Natator depressus* nesting sites protectedy;  
*CC_mass:* Mean *Caretta caretta* body mass in kg;  
*CM_mass:* Mean *Chelonia mydas* body mass in kg;  
*DC_mass:* Mean *Dermochelys coriacea* body mass in kg;  
*EI_mass:* Mean *Eretmochelys imbricata* body mass in kg;  
*LK_mass:* Mean *Lepidochelys kempii* body mass in kg;  
*LO_mass:* Mean *Lepidochelys olivacea* body mass in kg;  
*ND_mass:* Mean *Natator depressus* body mass in kg. 


## data/world_turtle_grouped.csv
### Variables descriptions  
*name:* Country name;  
*species:* Binomial nomenclature for species (as *Genus_latin*);  
*protected_perc:* Percentage of nesting sites protected by country and species;  
*mean_temp:* Mean nest temperature by each country (°C);  
*mass_mean:* Mean body mass in kg;  
*sd_mass:* Standard diviation of body mass;  
*mass_n:* Number of mass measurements extracted by species and country;  


# Turtle phylogeny data
## data/world turtle.csv
Phylogenetic co-variance matrix of all seven sea turtle species used for the analysis `data/phylo_cor.rds`.

## data/world turtle.csv
Phylogenetic tree reconstruction of all seven sea turtle species used for the analysis `data/phylo_tree.rds`.

    
  


