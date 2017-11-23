# Name: R_Get_EEnergy.R
# Author: Alvin Farrel
# Description: Calculates static electrostatic Energy in the system
# Usage: Rscript R_Get_EEnergy.R Your_File.pdb

args <- commandArgs(TRUE)
PDB_filename = args

PDB = readLines(PDB_filename)

#Get all Atom Records from PDB file
PDB_Atoms = subset(PDB,substr(PDB,1,4)=="ATOM")

#Remove extra white space from PDB records
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


#Get all Amino acid atoms
PDB_Atoms_AA = subset(PDB_Atoms,nchar(trim(substr(PDB_Atoms,18,20)))==3)
#Get all DNA atoms
PDB_Atoms_DNA = subset(PDB_Atoms,nchar(trim(substr(PDB_Atoms,18,20)))==2)

# Group all Histidine atoms
PDB_Atoms_AA_His_ND1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="HIS" & trim(substr(PDB_Atoms_AA,13,16))=="ND1") #-0.255
PDB_Atoms_AA_His_HD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="HIS" & trim(substr(PDB_Atoms_AA,13,16))=="HD2") #0.068
PDB_Atoms_AA_His_HE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="HIS" & trim(substr(PDB_Atoms_AA,13,16))=="HE1") #0.107
PDB_Atoms_AA_His_HE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="HIS" & trim(substr(PDB_Atoms_AA,13,16))=="HE2") #0.217

# Assign Histidne charges to atoms an create Histidine charge data.frame from the TF-DNA complex
max.len = max(length(PDB_Atoms_AA_His_ND1),length(PDB_Atoms_AA_His_HD2),length(PDB_Atoms_AA_His_HE1),length(PDB_Atoms_AA_His_HE2))
PDB_Atoms_AA_His = data.frame(PDB_Atoms_AA_His_ND1[1:max.len],PDB_Atoms_AA_His_HD2[1:max.len],PDB_Atoms_AA_His_HE1[1:max.len],PDB_Atoms_AA_His_HE2[1:max.len])
His_charges=c(-0.255,0.068,0.107,0.217)

# Group all Tyrosine atoms
PDB_Atoms_AA_TYR_CZ = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CZ") #-0.249
PDB_Atoms_AA_TYR_CG = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CG") #-0.249
PDB_Atoms_AA_TYR_CD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CD1") #0.062
PDB_Atoms_AA_TYR_CD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CD2") #0.062
PDB_Atoms_AA_TYR_CE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CE1") #0.065
PDB_Atoms_AA_TYR_CE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="CE2") #0.065
PDB_Atoms_AA_TYR_HD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="HD1") #0.062
PDB_Atoms_AA_TYR_HD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="HD2") #0.062
PDB_Atoms_AA_TYR_HE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="HE1") #0.065
PDB_Atoms_AA_TYR_HE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="HE2") #0.065
PDB_Atoms_AA_TYR_HH = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TYR" & trim(substr(PDB_Atoms_AA,13,16))=="HH") #0.218

#Assign Tyrosine atomic charges and create Tyrosine atomic charge dataframe
max.len = max(length(PDB_Atoms_AA_TYR_CZ),length(PDB_Atoms_AA_TYR_CG),length(PDB_Atoms_AA_TYR_CD1),length(PDB_Atoms_AA_TYR_CD2),length(PDB_Atoms_AA_TYR_CE1),length(PDB_Atoms_AA_TYR_CE2),length(PDB_Atoms_AA_TYR_HD1),length(PDB_Atoms_AA_TYR_HD2),length(PDB_Atoms_AA_TYR_HE1),length(PDB_Atoms_AA_TYR_HE2),length(PDB_Atoms_AA_TYR_HH))
PDB_Atoms_AA_TYR = data.frame(PDB_Atoms_AA_TYR_CZ[1:max.len],PDB_Atoms_AA_TYR_CG[1:max.len],PDB_Atoms_AA_TYR_CD1[1:max.len],PDB_Atoms_AA_TYR_CD2[1:max.len],PDB_Atoms_AA_TYR_CE1[1:max.len],PDB_Atoms_AA_TYR_CE2[1:max.len],PDB_Atoms_AA_TYR_HD1[1:max.len],PDB_Atoms_AA_TYR_HD2[1:max.len],PDB_Atoms_AA_TYR_HE1[1:max.len],PDB_Atoms_AA_TYR_HE2[1:max.len],PDB_Atoms_AA_TYR_HH[1:max.len])
Tyr_charges = c((-0.249*6),(-0.249*6),(-0.249*6),(-0.249*6),(-0.249*6),(-0.249*6),0.062,0.062,0.065,0.065,0.218)

# Group all Phenylalanine atoms
PDB_Atoms_AA_PHE_CZ = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CZ") #-0.372
PDB_Atoms_AA_PHE_CG = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CG") #-0.372
PDB_Atoms_AA_PHE_CD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CD1") #0.062
PDB_Atoms_AA_PHE_CD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CD2") #0.062
PDB_Atoms_AA_PHE_CE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CE1") #0.062
PDB_Atoms_AA_PHE_CE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="CE2") #0.062
PDB_Atoms_AA_PHE_HD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="HD1") #0.062
PDB_Atoms_AA_PHE_HD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="HD2") #0.062
PDB_Atoms_AA_PHE_HE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="HE1") #0.062
PDB_Atoms_AA_PHE_HE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="HE2") #0.062
PDB_Atoms_AA_PHE_HZ = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="PHE" & trim(substr(PDB_Atoms_AA,13,16))=="HZ") #0.062

#Assign Phenylalanine atomic charges and create Phenylalanine atomic charge dataframe
max.len = max(length(PDB_Atoms_AA_PHE_CZ),length(PDB_Atoms_AA_PHE_CG),length(PDB_Atoms_AA_PHE_CD1),length(PDB_Atoms_AA_PHE_CD2),length(PDB_Atoms_AA_PHE_CE1),length(PDB_Atoms_AA_PHE_CE2),length(PDB_Atoms_AA_PHE_HD1),length(PDB_Atoms_AA_PHE_HD2),length(PDB_Atoms_AA_PHE_HE1),length(PDB_Atoms_AA_PHE_HE2),length(PDB_Atoms_AA_PHE_HZ))
PDB_Atoms_AA_PHE = data.frame(PDB_Atoms_AA_PHE_CZ[1:max.len],PDB_Atoms_AA_PHE_CG[1:max.len],PDB_Atoms_AA_PHE_CD1[1:max.len],PDB_Atoms_AA_PHE_CD2[1:max.len],PDB_Atoms_AA_PHE_CE1[1:max.len],PDB_Atoms_AA_PHE_CE2[1:max.len],PDB_Atoms_AA_PHE_HD1[1:max.len],PDB_Atoms_AA_PHE_HD2[1:max.len],PDB_Atoms_AA_PHE_HE1[1:max.len],PDB_Atoms_AA_PHE_HE2[1:max.len],PDB_Atoms_AA_PHE_HZ[1:max.len])
Phe_charges = c(-0.372,-0.372,-0.372,-0.372,-0.372,-0.372,0.062,0.062,0.062,0.062,0.062)

# Group all Tryptophan atoms
PDB_Atoms_AA_TRP_CG = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TRP" & trim(substr(PDB_Atoms_AA,13,16))=="CG") #-0.685
PDB_Atoms_AA_TRP_CH2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TRP" & trim(substr(PDB_Atoms_AA,13,16))=="CH2") #-0.685
PDB_Atoms_AA_TRP_HD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TRP" & trim(substr(PDB_Atoms_AA,13,16))=="HD1") #0.104
PDB_Atoms_AA_TRP_HE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TRP" & trim(substr(PDB_Atoms_AA,13,16))=="HE1") #0.252
PDB_Atoms_AA_TRP_HZ2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="TRP" & trim(substr(PDB_Atoms_AA,13,16))=="HZ2") #0.054

#Assign Tryptophan atomic charges and create Tryptophan atomic charge dataframe
max.len = max(length(PDB_Atoms_AA_TRP_CG),length(PDB_Atoms_AA_TRP_CH2),length(PDB_Atoms_AA_TRP_HD1),length(PDB_Atoms_AA_TRP_HE1),length(PDB_Atoms_AA_TRP_HZ2))
PDB_Atoms_AA_TRP = data.frame(PDB_Atoms_AA_TRP_CG[1:max.len],PDB_Atoms_AA_TRP_CH2[1:max.len],PDB_Atoms_AA_TRP_HD1[1:max.len],PDB_Atoms_AA_TRP_HE1[1:max.len],PDB_Atoms_AA_TRP_HZ2[1:max.len])
Trp_charges = c(-0.685,-0.685,0.104,0.252,0.054)

# Group all Aspartate atoms
PDB_Atoms_AA_ASP_OD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ASP" & trim(substr(PDB_Atoms_AA,13,16))=="OD1") #-1.0 #-0.323
PDB_Atoms_AA_ASP_OD2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ASP" & trim(substr(PDB_Atoms_AA,13,16))=="OD2") #-1.0 #-0.323

# Assign Aspartate charges create Aspartate  dataframe
max.len = max(length(PDB_Atoms_AA_ASP_OD1),length(PDB_Atoms_AA_ASP_OD2))
PDB_Atoms_AA_ASP = data.frame(PDB_Atoms_AA_ASP_OD1[1:max.len],PDB_Atoms_AA_ASP_OD2[1:max.len])
Asp_charges = c(((-0.323 -1.0)/2),((-0.323 -1.0)/2))

# Group all Glutamate atoms
PDB_Atoms_AA_GLU_OE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="GLU" & trim(substr(PDB_Atoms_AA,13,16))=="OE1") #-1.0 #-0.322
PDB_Atoms_AA_GLU_OE2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="GLU" & trim(substr(PDB_Atoms_AA,13,16))=="OE2") #-1.0 #-0.322 

# Assign Glutamate charges create Glutamate  dataframe
max.len = max(length(PDB_Atoms_AA_GLU_OE1),length(PDB_Atoms_AA_GLU_OE2))
PDB_Atoms_AA_GLU = data.frame(PDB_Atoms_AA_GLU_OE1[1:max.len],PDB_Atoms_AA_GLU_OE2[1:max.len])
Glu_charges = c(((-0.322 -1.0)/2),((-0.322 -1.0)/2))

# Group all Glutamine atoms
PDB_Atoms_AA_GLN_OE1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="GLN" & trim(substr(PDB_Atoms_AA,13,16))=="OE1") #-0.466
PDB_Atoms_AA_GLN_HE21 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="GLN" & (trim(substr(PDB_Atoms_AA,13,16))=="HE21" | trim(substr(PDB_Atoms_AA,13,16))== "1HE2")) #0.157
PDB_Atoms_AA_GLN_HE22 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="GLN" & (trim(substr(PDB_Atoms_AA,13,16))=="HE22" | trim(substr(PDB_Atoms_AA,13,16))== "2HE2")) #0.157

# Assign Glutamate charges create Glutamine  dataframe
max.len = max(length(PDB_Atoms_AA_GLN_OE1),length(PDB_Atoms_AA_GLN_HE21),length(PDB_Atoms_AA_GLN_HE22))
PDB_Atoms_AA_GLN = data.frame(PDB_Atoms_AA_GLN_OE1[1:max.len],PDB_Atoms_AA_GLN_HE21[1:max.len],PDB_Atoms_AA_GLN_HE22[1:max.len])
Gln_charges = c(-0.466,0.157,0.157)

# Group all Asparagine atoms
PDB_Atoms_AA_ASN_OD1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ASN" & trim(substr(PDB_Atoms_AA,13,16))=="OD1") #-0.466
PDB_Atoms_AA_ASN_HD21 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ASN" & (trim(substr(PDB_Atoms_AA,13,16))=="HD21" | trim(substr(PDB_Atoms_AA,13,16))=="1HD2")) #0.157
PDB_Atoms_AA_ASN_HD22 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ASN" & (trim(substr(PDB_Atoms_AA,13,16))=="HD22" | trim(substr(PDB_Atoms_AA,13,16))== "2HD2")) #0.157
# Assign Glutamate charges create Asparagine  dataframe
max.len = max(length(PDB_Atoms_AA_ASN_OD1),length(PDB_Atoms_AA_ASN_HD21),length(PDB_Atoms_AA_ASN_HD22))
PDB_Atoms_AA_ASN = data.frame(PDB_Atoms_AA_ASN_OD1[1:max.len],PDB_Atoms_AA_ASN_HD21[1:max.len],PDB_Atoms_AA_ASN_HD22[1:max.len])
Asn_charges = c(-0.466,0.157,0.157)

# Group all Serine atoms
PDB_Atoms_AA_SER_HG = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="SER" & trim(substr(PDB_Atoms_AA,13,16))=="HG") #0.21
PDB_Atoms_AA_SER_HB2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="SER" & trim(substr(PDB_Atoms_AA,13,16))=="HB2") #0.056
PDB_Atoms_AA_SER_HB3 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="SER" & trim(substr(PDB_Atoms_AA,13,16))=="HB3") #0.056

# Assign Serine charges create Serine  dataframe
max.len = max(length(PDB_Atoms_AA_SER_HG),length(PDB_Atoms_AA_SER_HB2),length(PDB_Atoms_AA_SER_HB3))
PDB_Atoms_AA_SER = data.frame(PDB_Atoms_AA_SER_HG[1:max.len],PDB_Atoms_AA_SER_HB2[1:max.len],PDB_Atoms_AA_SER_HB3[1:max.len])
Ser_charges = c(0.21,0.056,0.056)

# Group all Threonine atoms
PDB_Atoms_AA_THR_HG1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="THR" & trim(substr(PDB_Atoms_AA,13,16))=="HG1") #0.21
PDB_Atoms_AA_THR_HB = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="THR" & trim(substr(PDB_Atoms_AA,13,16))=="HB") #0.059

# Assign Threonine charges create Threonine  dataframe
max.len = max(length(PDB_Atoms_AA_THR_HG1),length(PDB_Atoms_AA_THR_HB))
PDB_Atoms_AA_THR = data.frame(PDB_Atoms_AA_THR_HG1[1:max.len],PDB_Atoms_AA_THR_HB[1:max.len])
Thr_charges = c(0.21,0.059)

# Group all Lysine atoms
PDB_Atoms_AA_LYS_HZ1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="LYS" & trim(substr(PDB_Atoms_AA,13,16))=="HZ1") #1.0 #0.18
PDB_Atoms_AA_LYS_HZ2 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="LYS" & trim(substr(PDB_Atoms_AA,13,16))=="HZ2") #1.0 #0.18
PDB_Atoms_AA_LYS_HZ3 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="LYS" & trim(substr(PDB_Atoms_AA,13,16))=="HZ3") #1.0

# Assign Lysine charges create Lysine  dataframe
max.len = max(length(PDB_Atoms_AA_LYS_HZ1),length(PDB_Atoms_AA_LYS_HZ2),length(PDB_Atoms_AA_LYS_HZ3))
PDB_Atoms_AA_LYS = data.frame(PDB_Atoms_AA_LYS_HZ1[1:max.len],PDB_Atoms_AA_LYS_HZ2[1:max.len],PDB_Atoms_AA_LYS_HZ3[1:max.len])
Lys_charges = c((0.18 + 1.0/3),(0.18 + 1.0/3),(0.18 + 1.0/3))
#Lys_charges = c(0.18,0.18,0.18)

# Group all Arginine atoms
PDB_Atoms_AA_ARG_HE = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ARG" & trim(substr(PDB_Atoms_AA,13,16))=="HE") #1.0 #0.152
PDB_Atoms_AA_ARG_HH11 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ARG" & (trim(substr(PDB_Atoms_AA,13,16))=="HH11" | trim(substr(PDB_Atoms_AA,13,16))== "1HH1")) #1.0 #0.148
PDB_Atoms_AA_ARG_HH12 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ARG" & (trim(substr(PDB_Atoms_AA,13,16))=="HH12" | trim(substr(PDB_Atoms_AA,13,16))== "2HH1"))#1.0 #0.148
PDB_Atoms_AA_ARG_HH21 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ARG" & (trim(substr(PDB_Atoms_AA,13,16))=="HH21" | trim(substr(PDB_Atoms_AA,13,16))== "1HH2"))#1.0 #0.148
PDB_Atoms_AA_ARG_HH22 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="ARG" & (trim(substr(PDB_Atoms_AA,13,16))=="HH22" | trim(substr(PDB_Atoms_AA,13,16))== "2HH2"))#1.0 #0.148

# Assign charges create Arginine  dataframe
max.len = max(length(PDB_Atoms_AA_ARG_HE),length(PDB_Atoms_AA_ARG_HH11),length(PDB_Atoms_AA_ARG_HH12),length(PDB_Atoms_AA_ARG_HH21),length(PDB_Atoms_AA_ARG_HH22))
PDB_Atoms_AA_ARG = data.frame(PDB_Atoms_AA_ARG_HE[1:max.len],PDB_Atoms_AA_ARG_HH11[1:max.len],PDB_Atoms_AA_ARG_HH12[1:max.len],PDB_Atoms_AA_ARG_HH21[1:max.len],PDB_Atoms_AA_ARG_HH22[1:max.len])
Arg_charges = c((0.152 + 1/5),(0.148 + 1/5),(0.148 + 1/5),(0.148 + 1/5),(0.148 + 1/5))
#Arg_charges = c(0.152,0.148,0.148,0.148,0.148)

# Group all Cysteine atoms
PDB_Atoms_AA_CYS_HG1 = subset(PDB_Atoms_AA,trim(substr(PDB_Atoms_AA,18,20))=="CYS" & trim(substr(PDB_Atoms_AA,13,16))=="HG") #0.102

# Assign charges create Cysteine  dataframe
max.len = max(length(PDB_Atoms_AA_CYS_HG1))
PDB_Atoms_AA_CYS = data.frame(PDB_Atoms_AA_CYS_HG1[1:max.len])
Cys_charges = c(0.102)

#Group Guanine atoms and assign charges
PDB_Atoms_DNA_DG_N3 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="N3") #-0.206
PDB_Atoms_DNA_DG_H21 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="H21") #0.152
PDB_Atoms_DNA_DG_H22 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="H22") #0.152
PDB_Atoms_DNA_DG_H1 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="H1") #0.174
PDB_Atoms_DNA_DG_O6 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="O6") #-0.441
PDB_Atoms_DNA_DG_N7 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="N7") #-0.215
PDB_Atoms_DNA_DG_H8 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DG" & trim(substr(PDB_Atoms_DNA,13,16))=="H8") #0.107

max.len = max(length(PDB_Atoms_DNA_DG_N3),length(PDB_Atoms_DNA_DG_H21),length(PDB_Atoms_DNA_DG_H22),length(PDB_Atoms_DNA_DG_H1),length(PDB_Atoms_DNA_DG_O6),length(PDB_Atoms_DNA_DG_N7),length(PDB_Atoms_DNA_DG_H8))
PDB_Atoms_DNA_DG = data.frame(PDB_Atoms_DNA_DG_N3[1:max.len],PDB_Atoms_DNA_DG_H21[1:max.len],PDB_Atoms_DNA_DG_H22[1:max.len],PDB_Atoms_DNA_DG_H1[1:max.len],PDB_Atoms_DNA_DG_O6[1:max.len],PDB_Atoms_DNA_DG_N7[1:max.len],PDB_Atoms_DNA_DG_H8[1:max.len])
DG_charges = c(-0.206,0.152,0.152,0.174,-0.441,-0.215,0.107)

#Group Adenine atoms and assign charges
PDB_Atoms_DNA_DA_N3 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="N3") #-0.239
PDB_Atoms_DNA_DA_N1 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="N1") #-0.241
PDB_Atoms_DNA_DA_H2 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="H2") #0.089
PDB_Atoms_DNA_DA_H61 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="H61") #0.157
PDB_Atoms_DNA_DA_H62 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="H62") #0.157
PDB_Atoms_DNA_DA_N7 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="N7") #-0.21
PDB_Atoms_DNA_DA_H8 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DA" & trim(substr(PDB_Atoms_DNA,13,16))=="H8") #0.115

max.len = max(length(PDB_Atoms_DNA_DA_N3),length(PDB_Atoms_DNA_DA_N1),length(PDB_Atoms_DNA_DA_H2),length(PDB_Atoms_DNA_DA_H61),length(PDB_Atoms_DNA_DA_H62),length(PDB_Atoms_DNA_DA_N7),length(PDB_Atoms_DNA_DA_H8))
PDB_Atoms_DNA_DA = data.frame(PDB_Atoms_DNA_DA_N3[1:max.len],PDB_Atoms_DNA_DA_N1[1:max.len],PDB_Atoms_DNA_DA_H2[1:max.len],PDB_Atoms_DNA_DA_H61[1:max.len],PDB_Atoms_DNA_DA_H62[1:max.len],PDB_Atoms_DNA_DA_N7[1:max.len],PDB_Atoms_DNA_DA_H8[1:max.len])
DA_charges = c(-0.239,-0.241,0.089,0.157,0.157,-0.21,0.115)

#Group Cytosine atoms and assign charges
PDB_Atoms_DNA_DC_O2 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="O2") #-0.468
PDB_Atoms_DNA_DC_N3 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="N3") #-0.22
PDB_Atoms_DNA_DC_H41 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="H41") #0.157
PDB_Atoms_DNA_DC_H42 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="H42") #0.157
PDB_Atoms_DNA_DC_H5 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="H5") #0.066
PDB_Atoms_DNA_DC_H6 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DC" & trim(substr(PDB_Atoms_DNA,13,16))=="H6") #0.085

max.len = max(length(PDB_Atoms_DNA_DC_O2),length(PDB_Atoms_DNA_DC_N3),length(PDB_Atoms_DNA_DC_H41),length(PDB_Atoms_DNA_DC_H42),length(PDB_Atoms_DNA_DC_H5),length(PDB_Atoms_DNA_DC_H6))
PDB_Atoms_DNA_DC = data.frame(PDB_Atoms_DNA_DC_O2[1:max.len],PDB_Atoms_DNA_DC_N3[1:max.len],PDB_Atoms_DNA_DC_H41[1:max.len],PDB_Atoms_DNA_DC_H42[1:max.len],PDB_Atoms_DNA_DC_H5[1:max.len],PDB_Atoms_DNA_DC_H6[1:max.len])
DC_charges = c(-0.468,-0.22,0.157,0.157,0.066,0.085)

#Group Thymine atoms and assign charges
PDB_Atoms_DNA_DT_O2 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="O2") #-0.373
PDB_Atoms_DNA_DT_H3 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="H3") #0.194
PDB_Atoms_DNA_DT_O4 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="O4") #-0.478
PDB_Atoms_DNA_DT_H6 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="H6") #0.096

max.len = max(length(PDB_Atoms_DNA_DT_O2),length(PDB_Atoms_DNA_DT_H3),length(PDB_Atoms_DNA_DT_O4),length(PDB_Atoms_DNA_DT_H6))
PDB_Atoms_DNA_DT = data.frame(PDB_Atoms_DNA_DT_O2[1:max.len],PDB_Atoms_DNA_DT_H3[1:max.len],PDB_Atoms_DNA_DT_O4[1:max.len],PDB_Atoms_DNA_DT_H6[1:max.len])
DT_charges = c(-0.373,0.194,-0.478,0.096)

#Group Thymine methyl group atoms
PDB_Atoms_DNA_DTMe_H71 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="H71") #-0.373
PDB_Atoms_DNA_DTMe_H72 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="H72") #0.194
PDB_Atoms_DNA_DTMe_H73 = subset(PDB_Atoms_DNA,trim(substr(PDB_Atoms_DNA,18,20))=="DT" & trim(substr(PDB_Atoms_DNA,13,16))=="H73") #-0.478

max.len = max(length(PDB_Atoms_DNA_DTMe_H71),length(PDB_Atoms_DNA_DTMe_H72),length(PDB_Atoms_DNA_DTMe_H73))
PDB_Atoms_DNA_DTMe = data.frame(PDB_Atoms_DNA_DTMe_H71[1:max.len],PDB_Atoms_DNA_DTMe_H72[1:max.len],PDB_Atoms_DNA_DTMe_H73[1:max.len])





#########################################################
### This function creates a vector for all DNA Charges ##
#########################################################
Create.DNA.charge.vector <- function(PDB_Atoms_DNA_DA,DA_charges,PDB_Atoms_DNA_DC,DC_charges,PDB_Atoms_DNA_DG,DG_charges,PDB_Atoms_DNA_DT,DT_charges)
{ 
  DNA.charge.vector = c()
  
  for(j in 1:length(PDB_Atoms_DNA_DA[1,]))
    for(k in 1:length(PDB_Atoms_DNA_DA[,1])) DNA.charge.vector[(length(DNA.charge.vector)+1)] = DA_charges[j]              
  
  for(j in 1:length(PDB_Atoms_DNA_DC[1,])) 
    for(k in 1:length(PDB_Atoms_DNA_DC[,1])) DNA.charge.vector[(length(DNA.charge.vector)+1)] = DC_charges[j]               
  
  for(j in 1:length(PDB_Atoms_DNA_DG[1,]))
    for(k in 1:length(PDB_Atoms_DNA_DG[,1])) DNA.charge.vector[(length(DNA.charge.vector)+1)] = DG_charges[j]
  
  for(j in 1:length(PDB_Atoms_DNA_DT[1,]))
    for(k in 1:length(PDB_Atoms_DNA_DT[,1])) DNA.charge.vector[(length(DNA.charge.vector)+1)] = DT_charges[j]      
  
  return(DNA.charge.vector)
}

################################################################
### This function creates a vector for all Amino acid Charges ##
################################################################
ADD.AA.charge.vector <- function(PDB_Atoms_AA_NAME,AA_charges)
{
  
  AA.charge.vector = c()
  
  for(h in 1:length(PDB_Atoms_AA_NAME[1,]))
    for(i in 1:length(PDB_Atoms_AA_NAME[,1])) 
    {
      AA.charge.vector[(length(AA.charge.vector)+1)] = AA_charges[h]
     
    }
  
  return(AA.charge.vector)
}


#####################################################################
#### This function creates an XYZ matrix for all DNA coordiantes ####
## This Matrix will be used to determine the final Distance Matrix ##
#####################################################################
Get.DNA.coord.matrix <- function(PDB_Atoms_DNA_DA,PDB_Atoms_DNA_DC,PDB_Atoms_DNA_DG,PDB_Atoms_DNA_DT)
{
  
  DNA.coord.matrix = data.frame()
  temp.DNA.coord.matrix = data.frame()
         
      for(j in 1:length(PDB_Atoms_DNA_DA[1,]))
      {
        DXcoords =as.numeric(substr(PDB_Atoms_DNA_DA[,j],31,38)) 
        DYcoords =as.numeric(substr(PDB_Atoms_DNA_DA[,j],39,46)) 
        DZcoords =as.numeric(substr(PDB_Atoms_DNA_DA[,j],47,54)) 
        
        temp.DNA.coord.matrix=cbind(DXcoords,DYcoords,DZcoords)
        
        DNA.coord.matrix = rbind(DNA.coord.matrix,temp.DNA.coord.matrix)
                
      }
      for(j in 1:length(PDB_Atoms_DNA_DC[1,]))
      {
        
        DXcoords =as.numeric(substr(PDB_Atoms_DNA_DC[,j],31,38)) 
        DYcoords =as.numeric(substr(PDB_Atoms_DNA_DC[,j],39,46)) 
        DZcoords =as.numeric(substr(PDB_Atoms_DNA_DC[,j],47,54)) 
        
        temp.DNA.coord.matrix=cbind(DXcoords,DYcoords,DZcoords)
        
        DNA.coord.matrix = rbind(DNA.coord.matrix,temp.DNA.coord.matrix)
      }
      for(j in 1:length(PDB_Atoms_DNA_DG[1,]))
      {
              
        DXcoords =as.numeric(substr(PDB_Atoms_DNA_DG[,j],31,38)) 
        DYcoords =as.numeric(substr(PDB_Atoms_DNA_DG[,j],39,46)) 
        DZcoords =as.numeric(substr(PDB_Atoms_DNA_DG[,j],47,54)) 
        
        temp.DNA.coord.matrix=cbind(DXcoords,DYcoords,DZcoords)
        
        DNA.coord.matrix = rbind(DNA.coord.matrix,temp.DNA.coord.matrix)
        
      }
      for(j in 1:length(PDB_Atoms_DNA_DT[1,]))
      {
        DXcoords =as.numeric(substr(PDB_Atoms_DNA_DT[,j],31,38)) 
        DYcoords =as.numeric(substr(PDB_Atoms_DNA_DT[,j],39,46)) 
        DZcoords =as.numeric(substr(PDB_Atoms_DNA_DT[,j],47,54)) 
        
        temp.DNA.coord.matrix=cbind(DXcoords,DYcoords,DZcoords)
        
        DNA.coord.matrix = rbind(DNA.coord.matrix,temp.DNA.coord.matrix)
        
      }       
  return(DNA.coord.matrix)
}


########################################################
# Determine the XYZ matrix of the Thymine Mythle atoms #
########################################################
DTMe.coord.matrix = data.frame()
for(j in 1:length(PDB_Atoms_DNA_DTMe[1,]))
{
  DXcoords =as.numeric(substr(PDB_Atoms_DNA_DTMe[,j],31,38)) 
  DYcoords =as.numeric(substr(PDB_Atoms_DNA_DTMe[,j],39,46)) 
  DZcoords =as.numeric(substr(PDB_Atoms_DNA_DTMe[,j],47,54)) 
  
  temp.DTMe.coord.matrix=cbind(DXcoords,DYcoords,DZcoords)
  
  DTMe.coord.matrix = rbind(DTMe.coord.matrix,temp.DTMe.coord.matrix)
  
}  


####################################################################################################
# This function creates a Distance matrix between the AA record atoms ant an XYZ coordiante matrix #
####################################################################################################
Create.Dist.Matrix <- function(PDB_Atoms_AA_NAME,DNA.coords)
{
  
  Charge.Dist.matrix = data.frame()
  temp.Charge.Dist.vector = c()
  
  for(h in 1:length(PDB_Atoms_AA_NAME[1,]))
  {
    for(i in 1:length(PDB_Atoms_AA_NAME[,1]))
    {  
      temp.Charge.Dist.vector = c()
                      
        
        
        Xcoords =as.numeric(substr(PDB_Atoms_AA_NAME[i,h],31,38)) 
        Ycoords =as.numeric(substr(PDB_Atoms_AA_NAME[i,h],39,46)) 
        Zcoords =as.numeric(substr(PDB_Atoms_AA_NAME[i,h],47,54)) 
        
      temp.Charge.Dist.vector = sqrt((Xcoords-DNA.coords[,1])^2+(Ycoords-DNA.coords[,2])^2+(Zcoords-DNA.coords[,3])^2)

      #Check aromatic compounds as they have an electron cloud giving them a longer distance range of interactions.for this function
      if((grepl("TYR",PDB_Atoms_AA_NAME[i,h]) | grepl("PHE",PDB_Atoms_AA_NAME[i,h])) & grepl("C",substr(PDB_Atoms_AA_NAME[i,h],1,16))) 
      {
        temp.Charge.Dist.vector[temp.Charge.Dist.vector > 4.5] = Inf #If outside if distance cutoff it =infinity
      } else {temp.Charge.Dist.vector[temp.Charge.Dist.vector > 2.9] = Inf} #If outside if distance cutoff it = infinity
      
      Charge.Dist.matrix = rbind(Charge.Dist.matrix,temp.Charge.Dist.vector)  
    }        
  }
    
  
  colnames(Charge.Dist.matrix) <- c(1:length(Charge.Dist.matrix[1,]))
  return(Charge.Dist.matrix)
}


#Create DNA XYZ coordiante Matrix for all DNA atoms)
DNA.coords = Get.DNA.coord.matrix(PDB_Atoms_DNA_DA,PDB_Atoms_DNA_DC,PDB_Atoms_DNA_DG,PDB_Atoms_DNA_DT)

#Create DNA charge vector for all DNA atoms
DNA.charge.vector = Create.DNA.charge.vector(PDB_Atoms_DNA_DA,DA_charges,PDB_Atoms_DNA_DC,DC_charges,PDB_Atoms_DNA_DG,DG_charges,PDB_Atoms_DNA_DT,DT_charges)

DTMe.Dist.Matrix = matrix()
Dist.Matrix = matrix()
Dist.arom.Matrix = data.frame()
AA.charge.vector = c()
AA.arom.charge.vector = c()


#######################################################################
# Create Amino acid-DNA Distance matrix and amino acid charge vectors #
#######################################################################
Dist.Matrix = Create.Dist.Matrix(PDB_Atoms_AA_ASN,DNA.coords)
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_ASN,Asn_charges))  
DTMe.Dist.Matrix = Create.Dist.Matrix(PDB_Atoms_AA_ASN,DTMe.coord.matrix)

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_GLN,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_GLN,Gln_charges))  
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_GLN,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_ARG,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_ARG,Arg_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_ARG,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Miatrix,Create.Dist.Matrix(PDB_Atoms_AA_LYS,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_LYS,Lys_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_LYS,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_ASP,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_ASP,Asp_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_ASP,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_GLU,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_GLU,Glu_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_GLU,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_His,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_His,His_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_His,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_SER,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_SER,Ser_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_SER,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_THR,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_THR,Thr_charges))
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_THR,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_TYR,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_TYR,Tyr_charges)) 
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_TYR,DTMe.coord.matrix))

Dist.Matrix = rbind(Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_PHE,DNA.coords))
AA.charge.vector = c(AA.charge.vector,ADD.AA.charge.vector(PDB_Atoms_AA_PHE,Phe_charges)) 
DTMe.Dist.Matrix = rbind(DTMe.Dist.Matrix,Create.Dist.Matrix(PDB_Atoms_AA_PHE,DTMe.coord.matrix))
############################################################
############################################################
############################################################

##########################
#Check for steric Clashes#
Clash.rows = row(Dist.Matrix)[which(Dist.Matrix < 1.5)]
DTMe.Clash.rows = row(DTMe.Dist.Matrix)[which(DTMe.Dist.Matrix < 2.5)]

AA.charge.vector[Clash.rows] <- 0         #Assign charged of Clashing amino acids to 0
AA.charge.vector[DTMe.Clash.rows] <- 0    #so energy won't be skewed. Amino acids are dynamic so theortrically the AA can move
###########################               #and not cause an unfavorable energy so it is simple assigned to 0


Avogadro = 6.022140857*10^23
ke =  8.9875517873681764*10^9
dielectric = 4
electron = 1.602176565*10^-19
Lcutoff =1.5
cutoff = 2.9



temp.AA.charge.vector = AA.charge.vector
temp.DNA.charge.vector = DNA.charge.vector

EE_matrix = data.frame()
temp.Dist.Matrix = Dist.Matrix
temp.Dist.Matrix[is.na(temp.Dist.Matrix)] <- Inf

######################################################
## Use Distance Matrix and charge vectors to calculate electrostatic energy between each Amino Acid and Base atom using Coloumb's law
######################################################
for(i in 1:length(temp.Dist.Matrix[,1]))
{
  #Coloumb's law with charges and Constants#
  EE_matrix =rbind(EE_matrix,(ke*Avogadro*temp.AA.charge.vector[i]*electron*temp.DNA.charge.vector*electron)/(dielectric*1*10^-10)*2.39*10^-4)
}

clash_penalty = abs((ke*Avogadro*max(AA.charge.vector)*electron*min(DNA.charge.vector)*electron)/(dielectric*1.7*10^-10)*2.39*10^-4)

#Divide my distances to get Electostatic energy#
temp.EE_matrix = EE_matrix/temp.Dist.Matrix

#Add in Clash penalty for clashing atoms#
temp.EE_matrix[temp.Dist.Matrix < 1.5] <- clash_penalty

############################################
############################################

#Total Energy is the sum of all the Energies in the matrix#

temp.EEnergy = sum(temp.EE_matrix)
print(temp.EEnergy)


write(temp.EEnergy, file = "EEnergy.txt")

