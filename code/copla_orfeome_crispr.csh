#!/bin/tcsh
#$ -o /gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR/output_copla
#$ -e /gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR/error_copla

/gpfs0/tals/projects/software/copla/bin/copla.py /gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR/res/plsdb_unique.fasta \
    -a /gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR/data_calculations/plsdb_orfeome.fasta \
    /gpfs0/tals/projects/software/copla/databases/Copla_RS84/RS84f_sHSBM.pickle \
    /gpfs0/tals/projects/software/copla/databases/Copla_RS84/CoplaDB.lst \
    /gpfs0/tals/projects/Analysis/Lucy_plasmidome/Plasmidome/CRISPR/data_calculations/copla_CRISPR_orfeome