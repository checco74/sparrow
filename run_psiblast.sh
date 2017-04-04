#!/usr/bin/tcsh
# runs psiblast (blastpgp) for each of the input files provided
# writes one file with extension '.pssm' for each of the input files provided

# set the path to the blastpgp executable here
set BLAST_PATH = $PWD
# set number of blast iterations here
set ROUNDS = 3

foreach ifile ( $* )
	${BLAST_PATH} -j ${ROUNDS} -d ${BLAST_PATH}/nr -i ${ifile} -Q ${ifile}.psiblast
end
