1. From NCBI:
Download of available genomes from https://www.ncbi.nlm.nih.gov/genome/browse/# in file 'genomes_overview.txt' on November 30, 2017. There are 855 listed species (but for some, a proteome has not been predicted). Download of taxonomic information and (when available) proteome sequence for these 855 species:

``for species in `grep -v '^#' genomes_overview.txt | awk -F '\t' '{print $1}' | sed 's| |_|g'`;do wget -O esearch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy\&term=`echo $species | sed 's|_|+|g'`;id=`grep '^<Id>' esearch_$species | sed -e 's|^<Id>||' -e 's|</Id>||'`;wget -O efetch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy\&id=$id;sleep 0.5;done;for species in `grep -v '^#' genomes_overview.txt | awk -F '\t' '{print $1}' | sed 's| |_|g'`;do wget -O $species'_proteome.fa.gz' ftp://ftp.ncbi.nlm.nih.gov/genomes/$species/protein/protein.fa.gz;done``

For 328 species, predicted proteome sequence had been successfully downloaded from NCBI; for the remaining 527 species, it could not:

``for f in `ls *proteome*`;do if test -s $f;then echo $f;fi;done # these 328 species are OK``
``for f in `ls *proteome*`;do if [ ! -s $f ];then echo $f;fi;done # these 527 species are not OK: I delete these empty files:``
``for f in `ls *proteome*`;do if [ ! -s $f ];then rm -f $f;fi;done``

For 33 of the missing species: their predicted proteomes were available from a former download (on August 6, 2015, from NCBI, VectorBase, FlyBase, JGI, Ensembl and WormBase): these 33 old proteomes were included in the list.


2. From VectorBase:
With an interactive web browser on https://www.vectorbase.org/downloads?field_organism_taxonomy_tid%5B%5D=372&field_organism_taxonomy_tid%5B%5D=373&field_organism_taxonomy_tid%5B%5D=519&field_organism_taxonomy_tid%5B%5D=1505&field_organism_taxonomy_tid%5B%5D=364&field_organism_taxonomy_tid%5B%5D=556&field_organism_taxonomy_tid%5B%5D=549&field_organism_taxonomy_tid%5B%5D=1463&field_organism_taxonomy_tid%5B%5D=553&field_organism_taxonomy_tid%5B%5D=365&field_organism_taxonomy_tid%5B%5D=1533&field_organism_taxonomy_tid%5B%5D=550&field_organism_taxonomy_tid%5B%5D=554&field_organism_taxonomy_tid%5B%5D=847&field_organism_taxonomy_tid%5B%5D=366&field_organism_taxonomy_tid%5B%5D=367&field_organism_taxonomy_tid%5B%5D=1518&field_organism_taxonomy_tid%5B%5D=551&field_organism_taxonomy_tid%5B%5D=548&field_organism_taxonomy_tid%5B%5D=547&field_organism_taxonomy_tid%5B%5D=1534&field_organism_taxonomy_tid%5B%5D=1375&field_organism_taxonomy_tid%5B%5D=1521&field_organism_taxonomy_tid%5B%5D=546&field_organism_taxonomy_tid%5B%5D=880&field_organism_taxonomy_tid%5B%5D=1385&field_organism_taxonomy_tid%5B%5D=368&field_organism_taxonomy_tid%5B%5D=1241&field_organism_taxonomy_tid%5B%5D=1531&field_organism_taxonomy_tid%5B%5D=375&field_organism_taxonomy_tid%5B%5D=889&field_organism_taxonomy_tid%5B%5D=913&field_organism_taxonomy_tid%5B%5D=887&field_organism_taxonomy_tid%5B%5D=355&field_organism_taxonomy_tid%5B%5D=888&field_organism_taxonomy_tid%5B%5D=356&field_organism_taxonomy_tid%5B%5D=339&field_organism_taxonomy_tid%5B%5D=340&field_organism_taxonomy_tid%5B%5D=383&field_organism_taxonomy_tid%5B%5D=891&field_organism_taxonomy_tid%5B%5D=398&field_organism_taxonomy_tid%5B%5D=385&field_organism_taxonomy_tid%5B%5D=1343&field_organism_taxonomy_tid%5B%5D=393&field_organism_taxonomy_tid%5B%5D=1527&field_organism_taxonomy_tid%5B%5D=892&field_organism_taxonomy_tid%5B%5D=1557&field_organism_taxonomy_tid%5B%5D=1506&field_organism_taxonomy_tid%5B%5D=1556&field_download_file_type_tid%5B%5D=464&field_download_file_format_tid=All&field_status_value=Current

When several proteome files are available for a given species (e.g., from various isolates), pick the largest file (done on November 30, 2017); donwloaded proteomes: in the archive 'vbo_archive_20171130_0.zip':

``unzip vbo_archive_20171130_0.zip``

Renaming VectorBase files to fit with my nomenclature for other species (and in case I already had the proteome from the NCBI, simply delete that from VectorBase):

``for f in `ls *.fa.gz | grep -v '_proteome'`;do name=`echo $f | sed -e 's|^\([A-Z][a-z]*\)\-\([a-z]*\)-.*|\1_\2_proteome.fa.gz|'`;if [ ! -f $name ];then mv $f $name;else rm -f $f;fi;done``

(addition from VectorBase added another 32 proteomes)

3. From Uniprot:
On November 30, 2017: download of the list of eukaryotic proteomes available from Uniprot (http://www.uniprot.org/) in file 'proteomes-all.tab.gz'. Comparing that list to the list of remaining genomes for which the predicted proteome has not yet been downloaded, then (on December 1, 2017), download of the missing proteomes that Uniprot was offering:

``gunzip proteomes-all.tab.gz;for species in `grep -v '^#' genomes_overview.txt | awk -F '\t' '{print $1}' | sed 's| |_|g'`;do if [ ! -f $species'_proteome.fa'* ];then species_genus=`echo $species | sed 's|_.*||'`;species_species=`echo $species | sed 's|.*_||'`;if test `grep -w $species_genus proteomes-all.tab | grep -wc $species_species` -ne 0;then echo $species;fi;fi;done # There are 102 of my missing proteomes in Uniprot!``

``for species in `grep -v '^#' genomes_overview.txt | awk -F '\t' '{print $1}' | sed 's| |_|g'`;do if [ ! -f $species'_proteome.fa'* ];then species_genus=`echo $species | sed 's|_.*||'`;species_species=`echo $species | sed 's|.*_||'`;ID=`grep -w $species_genus proteomes-all.tab | grep -w $species_species | awk '{print $1}'`;if test "$ID" != "";then wget -O $species'_proteome.fa.gz' http://www.uniprot.org/uniprot/?query=proteome:$ID'&format=fasta&compress=yes';sleep 0.2;fi;fi;done``


4. From FlyBase:
On December 1, 2017: download of the Drosophila predicted proteomes from FlyBase:

``wget -O Drosophila_ananassae_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dana_r1.05/fasta/dana-all-translation-r1.05.fasta.gz;wget -O Drosophila_erecta_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dere_r1.05/fasta/dere-all-translation-r1.05.fasta.gz;wget -O Drosophila_grimshawi_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dgri_r1.05/fasta/dgri-all-translation-r1.05.fasta.gz;wget -O Drosophila_melanogaster_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dmel_r6.18/fasta/dmel-all-translation-r6.18.fasta.gz;wget -O Drosophila_mojavensis_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dmoj_r1.04/fasta/dmoj-all-translation-r1.04.fasta.gz;wget -O Drosophila_persimilis_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dper_r1.3/fasta/dper-all-translation-r1.3.fasta.gz;wget -O Drosophila_pseudoobscura_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dpse_r3.04/fasta/dpse-all-translation-r3.04.fasta.gz;wget -O Drosophila_sechellia_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dsec_r1.3/fasta/dsec-all-translation-r1.3.fasta.gz;wget -O Drosophila_simulans_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dsim_r2.02/fasta/dsim-all-translation-r2.02.fasta.gz;wget -O Drosophila_virilis_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dvir_r1.06/fasta/dvir-all-translation-r1.06.fasta.gz;wget -O Drosophila_willistoni_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dwil_r1.05/fasta/dwil-all-translation-r1.05.fasta.gz;wget -O Drosophila_yakuba_proteome.fa.gz ftp://ftp.flybase.net/releases/FB2017_05/dyak_r1.05/fasta/dyak-all-translation-r1.05.fasta.gz``


5. From JGI:
On December 1, 2017: comparison of the list in ftp://ftp.jgi-psf.org/pub/JGI_data/ with the list of proteomes already downloaded.Two additional proteomes are offered by JGI:

``wget -O Capitella_teleta_proteome.fa.gz ftp://ftp.jgi-psf.org/pub/JGI_data/Capitella/v1.0/allModels.aa.fasta.gz;wget -O Helobdella_robusta_proteome.fa.gz ftp://ftp.jgi-psf.org/pub/JGI_data/Helobdella_robusta/v1.0/proteins.Helro1_FilteredModels3.fasta.gz``

6. From Ensembl:
On December 1, 2017: download of the list of available proteomes from http://www.ensembl.org/index.html in file 'Available_at_Ensembl.txt'. Comparison with the list of already downloaded proteomes:

``for s in `sed -e 's|.*(||' -e 's|)$||' -e 's| |_|g' Available_at_Ensembl.txt`;do if test `ls *proteome.fa* | grep -c $s` -eq 0;then echo $s;fi;done``

Result: 91 of these were already downloaded (+ 4 with a variable name: "Gorilla gorilla gorilla" in Ensembl, "Gorilla gorilla" in NCBI; and 3 subspecies of "Mus musculus", which will be ignored because the Mus musculus proteome was already downloaded). But Ensembl offers 4 additional metazoan proteomes. Download:

``wget -O Mus_spretus_proteome.fa.gz ftp://ftp.ensembl.org/pub/release-90/fasta/mus_spretus_spreteij/pep/Mus_spretus_spreteij.SPRET_EiJ_v1.pep.all.fa.gz;wget -O Cavia_aperea_proteome.fa.gz ftp://ftp.ensembl.org/pub/release-90/fasta/cavia_aperea/pep/Cavia_aperea.CavAp1.0.pep.all.fa.gz;wget -O Peromyscus_maniculatus_proteome.fa.gz ftp://ftp.ensembl.org/pub/release-90/fasta/peromyscus_maniculatus_bairdii/pep/Peromyscus_maniculatus_bairdii.Pman_1.0.pep.all.fa.gz;wget -O Notamacropus_eugenii_proteome.fa.gz ftp://ftp.ensembl.org/pub/release-90/fasta/notamacropus_eugenii/pep/Notamacropus_eugenii.Meug_1.0.pep.all.fa.gz``


7. From WormBase:
On December 1, 2017: download of the list of available proteomes from http://www.wormbase.org/rest/widget/index/all/all/downloads?download=1&content-type=text%2Fhtml (file saved as 'From_wormbase.html').
Formatting of that list:

``tail -n +2 From_wormbase.html | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N 
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed '/^[A-Z][a-z]/ {
N
s|\(.*\)\n\(.*\)|\1 \2|
}' | sed -e 's|^\([A-Za-z]*\) \([a-z]*\).*protein.fa.gz <|\1_\2 |' -e 's|>$||' | grep -v '^$' | grep -v '^Link' > cmd_wormbase_dwnld``

Note that the 'cmd_wormbase_dwnld' file has to be corrected by hand for Caehnorhabditis sinica (weird formatting issues with the WormBase files caused it to point to the C. tropicalis proteome). Then: download (on December 1, 2017):

``for s in `awk '{print $1}' cmd_wormbase_dwnld | sort | uniq`;do url=`grep '^'$s' ' cmd_wormbase_dwnld | awk '{print $2}'`;wget -O $s'_proteome.fa.gz' $url;sleep 0.5;done``

For proteomes already downloaded from other sources: remove the WormBase version:

``for s in `ls *_proteome.fa.gz | sed 's|_proteome.fa.gz$||'`;do if [ -f $s'_proteome.fa' ];then echo "I already have "$s", I delete the WormBase file";rm -f $s'_proteome.fa.gz';fi;done``

For each species where it was not yet downloaded: download its taxonomic information:

``for s in `ls *_proteome.fa* | sed 's|_proteome.fa.*||'`;do if [ ! -f efetch_$s ];then echo $s;fi;done > Missing_taxonomy;for species in `cat Missing_taxonomy`;do wget -O esearch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy\&term=`echo $species | sed 's|_|+|g'`;id=`grep '^<Id>' esearch_$species | sed -e 's|^<Id>||' -e 's|</Id>||'`;wget -O efetch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy\&id=$id;sleep 0.5;done``

For two incompletely identified species ("Rhabditophanes_sp" and "Trichinella_sp"), manual correction is needed:

``rm -f efetch_Rhabditophanes_sp efetch_Trichinella_sp;for species in Rhabditophanes Trichinella;do wget -O esearch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy\&term=`echo $species | sed 's|_|+|g'`;id=`grep '^<Id>' esearch_$species | sed -e 's|^<Id>||' -e 's|</Id>||'`;wget -O efetch_$species http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy\&id=$id;sleep 0.5;done;rm -f Trichinella_proteome.fa.gz efetch_Trichinella``


8. Uncompressing (some files were gz-compressed):

``for f in `ls *.gz`;do gunzip $f;done``


9. Final list:
A predicted proteome was downloaded from 538 species in total. They are available in the archive:
https://www.igh.cnrs.fr/images/microsite/herve-seitz/files/pinzon-et-al-2018/00_Proteome_download//All_predicted_proteomes.tar.bz2


10. Measuring protein length distribution in each predicted proteome:

``for s in `ls *_proteome.fa | sed 's|_proteome\.fa$||'`;do ./Module_measures_proteome_size_dist.pl $s'_proteome.fa' | sed 's|.* ||' > Predicted_protein_lengths_$s'.dat';done``


11. Selection of proteomes with large numbers of long proteins (either at least 5,000 proteins of at least 500 amino acids; or at least 1,000 proteins of at least 1,000 amino acids):

``R CMD BATCH R_commands_statistics_on_protein_lengths;awk -F ',' '$5>=5000 {print $2}' Length_numbers.csv | sed 's|"||g' > Species_with_at_least_5000_proteins_of_more_than_500aa.txt;awk -F ',' '$7>=1000 {print $2}' Length_numbers.csv | sed 's|"||g' > Species_with_at_least_1000_proteins_of_more_than_1000aa.txt``
