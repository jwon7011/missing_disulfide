# Analysis of missing disulfide bonds in PDB files

The scripts in this repository are used to identify disulfide bonds that are missing in certain PDB chains either in within a PDB structure or in a different PDB structure. Currently the analysis is limited to the analysis of XRay structures.

Only the installation of Python is necessary to use the scripts. However, additional analysis files required include:

1. Directory containing all PDB files in *.ent.gz (ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/)
2. Chain level map file from PDBSWS (http://www.bioinf.org.uk/pdbsws/pdbsws_chain.txt.gz)
3. Residue level map file from PDBSWS (http://www.bioinf.org.uk/pdbsws/pdbsws_res.txt.gz)
4. DSB analysis saved in tab-delimited format (http://bioinf.med.unsw.edu.au/python/disulfideanalysis/download.html)

To run the analysis use the following commands as an example:

  python pdbDisulfide_interXray_gz -i ~/DB/pdb/  -m ~/DB/pdb/pdb_uniprot_chain_map.lst.2 -b ~/DB/pdb_uniprot_map.lst.2 -o interXRay_all_out.txt
  
  python annotatePDB.py interXRay_all_out.txt ~/DB/AllXray.txt interXray_annotated.txt
  
The output can be filtered based on column 8 which lists whether a bond is missing (SSMissing), mutated (SSMutated) or simply due to a chain sequence not covering the disulfide bond (NonOverlap). The status of whether a bond is missing in the same or different structure can be further filtered based on column 10 (0 = different PDB, 1 = same PDB).
