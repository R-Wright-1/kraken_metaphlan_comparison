---
title: "Making a custom database for Kraken2"
---

If you aren't interested in any of the rest of this file (why it's necessary, an explainer of the files used etc.), feel free to skip ahead to "Run the scripts" where you can either: (1) make the files ready to run Kraken2 with the [Genome Taxonomy Database]() (with the other domains from NCBI RefSeq if you wish); or (2) make the files ready to run Kraken2 with any other custom database.

## When is this useful?

For most applications, I recommend following the instructions on database construction that are on the [Kraken2 wiki page](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases), using the NCBI taxonomy. The aim here is to make a database using genomes for which you want to use a custom taxonomy. In my case, I wanted to use the [GTDB taxonomy](https://gtdb.ecogenomic.org/), and I tried doing this in a few ways before giving up and writing some scripts to do this for myself. The approaches that I tried first are listed here in case they are of use to anyone else:

1. Firstly, I downloaded the GTDB genomes (latest files [here](https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/)) and slightly modified the script that I have [here](https://github.com/R-Wright-1/peptides) (that downloads files from NCBI and modifies the fasta file headers to match the expected input for Kraken2 database building). The aim here was to modify all of the fasta files for the GTDB genomes to match what Kraken2 expected them to be. The issue that I ended up running into here was that the taxid's for about 5000 of these genomes weren't in the kraken taxonomy mapping file. These taxonomy ID's did generally exist, but weren't included in Kraken2 previously, presumably because they generally (for the ones that I looked at) weren't taxonomically classified. 

2. [GTDB_Kraken](https://github.com/hcdenbakker/GTDB_Kraken) (listed on the GTDB website) - following the build from source instructions. For this one, I modified the GTDB metadata file that I'd downloaded to have a .tsv file with the genome ID in the first column and the taxonomy (i.e., colon separated d__Kingdom;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species) in the second column, and then ran their script. However, I quickly realised that this was getting the NCBI taxonomy for the genomes, which clearly wasn't going to work for the same reason as the first approach. Maybe this worked when this github was last updated in 2018, or this wasn't an issue for what they wanted it for.... but this wasn't going to work either.

3. [Struo](https://github.com/leylabmpi/Struo) (interestingly, this was previously listed on the GTDB website and isn't any longer) - I tried to follow the instructions for using their gtdb_to_taxdump script, but ultimately kept running into errors that I was having problems trouble-shooting.

## Files needed to start

### So in order to build a database for Kraken2, we will need:

1. **Genomes that we want to add**</br>
   In my case, I first downloaded all GTDB representative genomes: https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz

2. **Full taxonomy information for each of these genomes**</br>
   For GTDB, this was both the bacterial and archael metadata files downloaded from: https://data.gtdb.ecogenomic.org/releases/latest/

### And what we want to get out of these is information on fields in each file from NCBI):

1. **names.dmp** - a tab delimited file with four columns: </br>
  *tax_id*					-- the id of node associated with this name</br>
  *name_txt*				-- name itself</br>
  *unique name*			-- the unique variant of this name if name not unique</br>
  *name class*			-- (synonym, common name, ...)</br>

2. **nodes.dmp** - a file with 13 fields that are both tab-delimited and separated with "|":</br>
  *tax_id*					                    -- node id in GenBank taxonomy database</br>
  *parent tax_id*		                    -- parent node id in GenBank taxonomy database</br>
  *rank*					                      -- rank of this node (superkingdom, kingdom, ...) </br>
  *embl code*				                    -- locus-name prefix; not unique</br>
  *division id*				                  -- see division.dmp file</br>
  *inherited div flag  (1 or 0)*		    -- 1 if node inherits division from parent</br>
  *genetic code id*				              -- see gencode.dmp file</br>
  *inherited GC flag  (1 or 0)*		      -- 1 if node inherits genetic code from parent</br>
  *mitochondrial genetic code id*		    -- see gencode.dmp file</br>
  *inherited MGC flag  (1 or 0)*		    -- 1 if node inherits mitochondrial gencode from parent</br>
  *GenBank hidden flag (1 or 0)*        -- 1 if name is suppressed in GenBank entry lineage</br>
  *hidden subtree root flag (1 or 0)*   -- 1 if this subtree has no sequence data yet</br>
  *comments*				                    -- free-text comments and citations</br>
  **And what is necessary from this is**: tax_id, parent tax_id, rank (phylum, class, ....), embl code=XX, division id=0, inherited div flag=0, genetic code id=11, inherited GC flag=1, mitochondrial genetic code id=1, inherited MGC flag=0, GenBank hidden flag=0, hidden subtree root flag=0

3. **nucl.taxid2accession** - a tab-delimited file with four columns:</br>
  *accession*           -- the genome accession number</br>
  *accession version*   -- the genome accession number with version</br>
  *taxid*               -- the taxid for that genome</br>
  *gi*                  -- this was in NCBI but these are being phased out, so while this column is required, it can be filled with "NA"</br>
  **As an example**:  accession       accession.version       taxid   gi</br>
                      A00001          A00001.1                10641   NA</br>
  So in the case of the GTDB genomes, the accession version is given in the first column (e.g. GB_GCA_000006155.2), and removing the final number (i.e. .2) will give the accession (e.g. GB_GCA_000006155). I imagine that if you only had a list of genome names, it would be fine to just add ".1" to the end of each to make this compatible. 

### So what we actually need to put together is:

- taxid - these are actually pretty simple and just need a unique number for each taxonomy and level of the taxonomy, and we need to keep a track of what the parent node is for each ID</br>
- name </br>
- rank </br>
- parent taxid (this can be itself in the case of root only)</br>

and we can then use these to make all of the other files. 

## Run the script

There are four scripts that can be run here, and you should run them in this order, but you don't need to run all of them depending on what you want to include in your database:

1. `get_ncbi_other_domains.py` - this can be used to download all genomes from any domain from the most recent NCBI RefSeq release. By default, this will download all *but* the bacteria and archaea. You can check available domains [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/).  It requires the following to be installed:
  - Python + default packages os, pandas, argparse, multiprocessing and sys
  - wget
  
  - You can run it as follows:
    `python get_ncbi_other_domains.py` - this will download all fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other and viral genomes (regardless of whether they are complete or not) to the folder ncbi_genomes/
    **OR** *e.g.* `python get_ncbi_other_domains.py --domain protozoa,vertebrate_mammalian --complete True --folder $FOLDER --download_genomes False --log_file logfile.txt --processors 12` - with any or all of these options for specific domains to download (as a comma separated list with no spaces), whether to download only complete genomes (True/False, False by default), the folder to download the genomes to (if this already exists then the additional genomes will be added to it), and whether to download the genomes or not (if False, then the information on the genomes within the domains given will still be downloaded and saves as summary_to_download.csv). The processors option determines how many processors to use for downloading genomes (default is 1). Any genomes that can't be downloaded will be written to the logfile (default logfile_download_genomes.txt)
    This will save summary_to_download.csv (containing assembly accessions, NCBI taxonomy IDs, organism name, assembly level, ftp path and domain/taxonomy information) as well as genomes.tsv (a tab-delimited file containing assembly accessions and taxonomy information needed for script 2)

2. `download_GTDB_latest.py` - this will download the latest GTDB representative genomes and will make the taxonomy structure needed for Kraken2 and output the nodes.dmp, names.dmp and nucl.taxid2accession files. You can also add your own custom list of genomes. It requires the following to be installed:
  - Python + default packages pandas, numpy, os, argparse, sys, subprocess
  - wget
  
  - You can run it as follows:
    `python download_GTDB_latest.py` OR
    `python download_GTDB_latest.py --download_genomes False` - if, for whatever reason, you don't want to actually download all of the GTDB genomes 
  - Note that you can also add your own genomes to this by having a file called genomes.tsv that is a tab-delimited file containing two columns - one containing the genome accession and the second containing the taxonomy in the format d__ ;k__ ;p__ ;c__ ;o__ ;f__ ;g__ ;s__ ;s1__ for domain, kingdom, phylum, class, order, family, genus, species and strain. You don't need to give this as an option, but if it exists in the folder this script is in then they will be added. You should then put the fasta files for these in the same folder as the GTDB genomes that have been downloaded. 
  This will save the files: genomes_and_taxonomy.tsv (a tab-delimited file containing the accession and taxonomy for all GTDB genomes, not needed for anything further), db_samples.tsv (a tab-delimited file containing four named columns, accession, gtdb_taxonomy, taxid and parent - this will be needed for script 4), db_taxid_info.tsv (a tab-delimited file containing four names columns, taxid, parent, name and rank - not needed for anything further), nodes.dmp, names.dmp and nucl.taxid2accession (files in the format described above, neede for construction of the Kraken2 database), folder containing all genomes, gtdb_genomes_reps/ (currently - as of May 2021 - there are 47,915 GTDB genomes). 

3. `get_taxid_custom_files.py` - this requires you to input a list of genomes and their taxonomy and will output the nodes.dmp, names.dmp and nucl.taxid2accession files needed to build the . It requires the following to be installed:
  - Python + default packages numpy, sys, argparse and pandas
  
  - You can run it as follows:
    `python get_taxid_custom_files.py --genome_taxonomy $FILE` where $FILE is a tab delimited .tsv file with two columns - one containing a list of genomes, and the second with their taxonomy (as is output by download_GTDB_latest.py - genomes_and_taxonomy.tsv)
    This will save the files: db_samples.tsv, db_taxid_info.tsv, names.dmp, nodes.dmp and nucl.taxid2accession (each as described for script 2)

4. `rename_fasta_headers.py` - after running either of script 2 or 3, you should have the nodes.dmp, names.dmp and nucl.taxid2accession files that kraken2 requires, as well as some additional files (db_samples.tsv and db_taxid_info.tsv) and a folder containing all of the genomes that you'd like to add. To be compatible with Kraken2 you now need to add the taxid to the fasta headers of each of your genome files. It requires the following to be installed:
  - Python + default packages os, argparse, multiprocessing and pandas
  - Biopython - [install instructions here](https://biopython.org/wiki/Download), or using conda:`conda install -c conda-forge biopython`
  
  - You can run it as follows:
    `python rename_fasta_headers.py --genome_folder $FOLDER --genome_list $FILE_ACC_TAXID --log_file logfile.txt --processors 12` where $FOLDER is a folder containing all of the genomes that you'd like to rename prior to adding to kraken2 and $FILE_ACC_TAXID is a tab-delimited file with a list of genome accessions in the first column (db_samples.tsv, output by the previous step) and has at least one other column called 'taxid' that contains the taxid that you'd like to add to this genome. The logfile will note any genomes that can't be renamed and the processors option allows you to rename multiple genomes at the same time.
  - If you just ran the `download_GTDB_latest.py` script then you don't need to give any options  and can simply run this with `python rename_fasta_headers.py`
  
### Which scripts to run

For example, to get a database containing all GTDB genomes and the NCBI RefSeq genomes for the other domains (fungi, invertebrate, plant, protozoa, vertebrate_mammalian, vertebrate_other and viral), I ran script 1 followed by script 2, moved the genomes from ncbi_genomes into gtdb_genomes_reps and then ran script 4. If you only wanted to download the bacteria and archaeal genomes from GTDB, you could just run scripts 2 and 4, and if you wanted to use only a custom set of genomes, you could run scripts 3 and 4. 

## Make the database

Assuming that these ran fine then you can make your kraken2 database by (see the [Kraken documentation on building custom databases](https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases) for more information):

1. Making a directory for your database: `mkdir kraken2_database`

2. Moving the files to the directory:   `mv nodes.dmp kraken2_database/`
                                      `mv names.dmp kraken2_database/`
                                      `mv nucl.taxid2accession kraken2_database/`

3. Adding all fasta files to your library: `for i in fasta_renamed/* ; do kraken2-build --add-to-library $i --db kraken2_database ; done`
   Note that files need to be unzipped for this to work. If you want to do this a bit faster then you can use the `unzip_add_library.py` script to gunzip, add files to the library and zip again, using as many processors as you like.

4. Building the database: `kraken2-build --build --db kraken2_database --threads 12 --max-db-size XX` where the maximum size is given in bytes. Kraken will downsample minimizers to fit this size. 

5. Building the files for Bracken: `bracken-build -d kraken2_database/ -t 12 -l 150` (If you leave the -l parameter then by default this will be 150)
