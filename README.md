# Development and calculation of the structural descriptors for nuclear receptors  
The project is aimed at exploring the possibility of using structural descriptors to distinguish patterns in protein binding methods, in particular nuclear receptors (NRs). The study is conducted on the example of Vitamin D Receptor (VDR), as a nuclear receptor with known structure and binding methods.  
### Project goal  
The aim is to create a library of structural descriptors and calculate them on the example of nuclear receptors.  
### Tasks of the project  
- To make a preliminary list of structural descriptors based on literature references and own suggestions;
- To form the final list of descriptors;
- To download protein structures of NRs from PDB;
- To create the calculation tools of structural descriptors for the PDB structures.
### Methods
Data for the project were obtained from ProteinDataBank (www.rcsb.org). It is 151 PDB-files, that contain information about three-dimensional structures of VDR.
Data preprocessing was performed by Schrodinger 2017-1. Library for calculation of descriptors was written using Python 3 and Bio.PDB package; information about secondary structures was gained using DSSP algorithm. Validation was performed by PyMOL.
Then the implemented descriptors' calculations were calculated for preprocessed and unprocessed structures and k-means clustering of the results was performed, visualisation was made by tSNE method.
Statistical analysis of obtained data was made on R, using such methods as nonparametric Wilcoxon test and Kruskal test both with Bonferroni correction.
### Usage
To calculate most descriptors (except ligand descriptors and sandwitch layer descriptor, user need to launch them separately) it's convenient to use ```Calculate_Descriptors.py```. For reproducibility of the analysis results we can offer the user to try following steps. 
1. Clone repository to personal mashine using ```git clone```

git clone https://github.com/prefabsprout/NuclReceptors_Descriptors.git

2. Then run the script from terminal. In general case command for running looks like that:

python3 ./NuclReceptors_Descriptors/src/Calculate_Descriptors.py {folder with PDB-files for species of interest} {name of output csv-file} —clamp-resid {three integers with - charge clamp residue numbers} —species {name of species: zebrafish/human/rat} —prepared {0 - unprocessed, 1 - preprocessed}

Examples:

python3 ./NuclReceptors_Descriptors/src/Calculate_Descriptors.py ./NuclReceptors_Descriptors/src/raw_data/VDR_PDB/Danio_rerio calc_results_danio.csv --clamp-resid 274 292 446 --species zebrafish --prepared 0

python3 ./NuclReceptors_Descriptors/src/Calculate_Descriptors.py ./NuclReceptors_Descriptors/src/raw_data/VDR_PDB/Homo_sapiens calc_results_homosap_prep.csv --clamp-resid 246 264 420 --species human --prepared 1

python3 ./NuclReceptors_Descriptors/src/Calculate_Descriptors.py ./NuclReceptors_Descriptors/src/raw_data/VDR_PDB/Rattus_norvegicus calc_results_rattus_prep.csv --clamp-resid 242 260 416 --species rat --prepared 1

Currently, the script for calculating the angle between sandwich layers and for calculation of solvent-accessibility area per helix need to be improved in order to facilitate user interaction with the script.
### Recomended system requirements for developed software
- Ubuntu 18.04
- Python 3.6
- Biopython 1.76
- Numpy 1.18.4
- Sympy 1.6
- Pandas 0.25.2
### Example of output file
The result of running is csv-table with names of descriptors in columns and PDB-structures in rows.
Resulting csv-tables can be found in ```Analysis/data```.
### References
1. rcsb.org. H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne.
(2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.
2. Stephen K Burley, Helen M. Berman, et al.
RCSB Protein Data Bank: biological macromolecular structures enabling research and education in fundamental biology, biomedicine, biotechnology and energy (2019) Nucleic Acids Research 47: D464–D474. doi: 10.1093/nar/gky1004.
3. Jain P., Hirst J. D. Study of protein structural descriptors: towards similarity and classification //Proceedings of Computational Biophysics to Systems Biology. – 2007. – С. 165-167.
4. Yaghmaei S. et al. Agonist and antagonist binding to the nuclear vitamin D receptor: dynamics, mutation effects and functional implications //In silico pharmacology. – 2013. – Т. 1. – №. 1. – С. 2.
5. Bourguet W., Germain P., Gronemeyer H. Nuclear receptor ligand-binding domains: three-dimensional structures, molecular interactions and pharmacological implications //Trends in pharmacological sciences. – 2000. – Т. 21. – №. 10. – С. 381-388.
6. Schrödinger Release 2020-2: Maestro, Schrödinger, LLC, New York, NY, 2020.
7. The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.
