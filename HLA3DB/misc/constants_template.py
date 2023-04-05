#       Sgourakis Lab
#   Author: Sagar Gupta
#   Date: April 28, 2022
#   Email: sagarg@sas.upenn.edu


# parameters
MUTATION_THRESH = 5 # number of allowed mutations in the HLA (inclusive)

LONG_PEP_LENGTH = 20 # Upper bound on what we consider a "peptide" chain
PEP_LOWER = 8 # Smallest peptide length to consider
PEP_UPPER = 10 # Longest peptide length to consider

HC_CUTOFF = 90 # Sequence identity score cutoff for HLA heavy chain (out of 180)
B2M_CUTOFF = 90 # Sequence identity score cutoff for HLA B2M chain (out of 100)

RECEPTOR_DIST = 5 # Distance (Å) threshold between peptide and receptor

DISTANCE_PEP_MHC = 15.0 # 65th residue of MHC Ca atom and 2nd residue for peptide Ca atom (from observation)

MHC_TRIM_LENGTH = 181 # Number of residues in MHC chain + 1 (so really it is 180)

FIRST_MHC_STRUC_YEAR = 1988

# reference values
A0201_SEQ = "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQ"
B2M_SEQ = "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"

DEFAULT_MHC_CHAIN = 'A'
DEFAULT_PEP_CHAIN = 'B'

MHC_SEQ_MOTIF_SHS = 'SHS'
MHC_SEQ_MOTIF_HS = 'HS'

DATABASE_SCRIPT_PATH = "/database/"
ROSETTA_INSTALL_DIR = "/rosetta/"

SINGLE_DIST_THRESHOLD = 1.0 # Distance score threshold for every single dihedral angle pair
SUM_DIST_THRESHOLD = 1.5 # Distance score threshold for every single dihedral angle pair

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
