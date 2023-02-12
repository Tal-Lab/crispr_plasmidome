
# Identification of potential hosts and the host range of plasmids based on alignment to CRISPR spacers

This computational method is used to identify potential hosts for newly discovered plasmids from metagenomic data and alternative hosts for plasmids identified in isolates. It is based on the alignment of CRISPR spacers and uses the BLAST+ tool and Python modules pandas and biopython. The method provides a prediction for the host range of the plasmids and shows the network of connection between different bacterial taxa.

## Requirements

- Python
- Pandas
- Biopython
- BLAST+ (makeblastdb, blastn)
- Unix tools
- COPLA (https://github.com/santirdnd/COPLA)
- Plasmid dataset in FASTA format (PLSDB used in this study: https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2)
- CRISPRCasdb spacers database (last updated version: https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=spacer_34.zip)
## Installation

1. Clone the source code from GitHub:

```bash
git clone https://github.com/Tal-Lab/crispr_plasmidome.git
cd crispr_plasmidome
```
2. Change directories in all coding files to the relevant ones.
    
## Usage and Examples

1. Finding matches:

```css
python path/to/code/Reaserch_project.py path/to/spacers.fasta path/to/plasmids.fasta

```
This script uses the BLAST+ tool to align sequences and outputs a 'Blast.csv' file in the working directory. Make sure to replace the path in the 'db' argument with the relevant directory before running the script.

2. Calculating alignment ratio:
```css
python path/to/code/alignment_ratio.py  path/to/Blast.csv path/to/Blast_ratio.csv
```
3. Finding the closest host match:

Replace directories in the script metadata.py and run the script search_taxonomy_match.py.

4. Calculating host range grade:
```css
python host_range.py /path/to/file/BLAST_database.csv
```

5. Preparing a table for network:
Change paths in network_table.py and uncomment necessary functions.

6. Analyze host range grades with regard to plasmid mobility:
- Obtain fasta with plasmids ORFs: run ORFeome_unique.py and replace paths as necessary.
- Get fasta sequences of unique plasmids with hits to spacers: run plsdb_match_fasta.py and replace paths as necessary.
- Run COPLA with copla_orfeome_crispr.csh and check input and output file paths. This will generate plasmid PTUs and MOB families.
- Run mobility_analysis.py to separate plasmids into mobilizable and non-mobilizable, add mobility


## Feedback

If you have any feedback, please reach out to us at lucyandrosyuk@gmail.com or leave a comment in Discussions.

Please, include “CRISPR github” in the subject and specify your issue. 



## Acknowledgements
- We thank Maria Pilar Garcillan Barcia for discussions and great ideas  
- Tal Shay for resources, discussions and advice 
- This study was supported (in part) by grant no. 3-17700 from the Office of the Chief Scientist, Israel Ministry of Health, as a part of the MAPMAR project (a part of the ERA-NET Cofund AquaticPollutants project). 
- Lucy Androsiuk is the recipient of a Hi-Tech, Bio-Tech, and Chemo-tech fellowship of Ben363 Gurion University of the Negev.
- The README was generated with the help of https://readme.so/editor
