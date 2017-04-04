#!/bin/tcsh
# runs sparrow on the given input file

set BDIR = "."
set WDIR = "."

set psiblast = "${BDIR}/run_psiblast.sh"

set globalstatus = 0

foreach arg ( $* )
	if ( "${arg}" == "-h" ) then 
		set helpme = 1
	endif
	if ( "${arg}" == "-fasta" ) then 
		set fastamode = 1
	endif
	set argvar = `echo "${arg}" | grep "^\-"`
	if ( "${argvar}" == '' ) then
		set inputfile = ${arg}
	endif
end
if ( ${?helpme} ) then
	echo "\n  usage:"
	echo "    run_sparrow-nn.sh [options] <input file>"
	echo "\n  options:"
	echo "    -fasta|-psiblast        specifies the type of input to be expected:"
	echo "                            if '-fasta' is specified psiblast will be run first;"
	echo "                            make sure it's installed and the nr database is in the right path."
	echo "    -width=<n>              <n> is the number of amino acids per line of prediction output."
	echo "                            [default <n>=40]"
	echo "    -dsspmode=loose|strict  DSSP reduction scheme. [default 'loose'; 'strict' not implemented.]"
	echo ""
	exit
endif
if ( ${?inputfile} ) then
	if ( "${WDIR}" != "${BDIR}" ) then
		echo "creating directory '${WDIR}'...."
		mkdir -p ${WDIR}/.prm
		ln -vs ${BDIR}/*.prm ${WDIR}
		ln -vs ${BDIR}/.prm/* ${WDIR}/.prm
		if ( -drwx ${WDIR} ) then
			cd ${WDIR}
		else
			echo "\ncouldn't access working directory ${WDIR}."
			exit
		endif
	endif
	if ( ${?fastamode} ) then
		echo "running PSI-BLAST...."
		${psiblast} ${inputfile}
	endif
	echo "running SPARROW stage-I...."
	${BDIR}/sparrow $*
	set globalstatus = $status
	set NDIR = ${BDIR}/annlight
	if ( -drx ${NDIR} && -fr stage2learn.base.cfg && ${globalstatus} == 0 ) then
		cat stage2learn.base.cfg > stage2learn.cfg
		if ( -fr ${inputfile}.nni ) then
			ls ${inputfile}.nni >> stage2learn.cfg
			echo "running SPARROW stage-II...."
			if ( -fr ${inputfile}.out ) then
				echo -n "" > ${inputfile}.nno
				echo "extracting sequence information from stage-I output file...."
				grep "sequence" ${inputfile}.out | sed 's/sequence >/seq >>/g' | split -l 1
				foreach sfile ( x* )
					echo "writing '${inputfile}.seq.${sfile}'...."
					mv ${sfile} ${inputfile}.seq.${sfile}
				end
				echo "invoking the artificial neural network...."
				${NDIR}/ann | grep -A 7 "Prediction for" | tail -1 | split -l 1
				echo "copying neural network output values...."
				if ( -fr scores.pred.csv ) then
					cp scores.pred.csv ${inputfile}.csv
					foreach sfile ( x* )
						echo "writing '${inputfile}.str.${sfile}'...."
						mv ${sfile} ${inputfile}.str.${sfile}
					end
					echo "computing confidence values...."
# 					${BDIR}/computeConfInd -l=15 scores.pred.csv
#					${BDIR}/computeConfInd -l=1 -sort scores.pred.csv
					${BDIR}/computeConfInd -mf scores.pred.csv
					echo "extracting confidence information...."
					grep "conf" scores.pred.csv.cnf | split -l 1
					foreach sfile ( x* )
						mv ${sfile} ${inputfile}.cnf.${sfile}
						cat ${inputfile}.seq.${sfile} >> ${inputfile}.nno
						cat ${inputfile}.str.${sfile} >> ${inputfile}.nno
						cat ${inputfile}.cnf.${sfile} >> ${inputfile}.nno
					end
					foreach junkfile ( \
							ann.data ann.debug \
							confidence-hit_table.dat \
							prediction.raw prediction_wrong.log \
							recall.raw recall_wrong.log \
							scores.csv scores.pred.csv scores.pred.csv.cnf \
							${inputfile}.cnf.* ${inputfile}.seq.* ${inputfile}.str.* ${inputfile}.nni \
						)
						rm -rf ${junkfile}
					end
					echo "finished."
				else
					echo "no neural network scores found." > /dev/stderr
				endif
				set globalstatus = $status
				rm ${inputfile}.out
			else
				echo "file ${inputfile}.out error." > /dev/stderr
			endif
		else
			echo "file '${inputfile}.nni' error." > /dev/stderr
		endif
	else
		echo "no neural network available." > /dev/stderr
	endif
	if ( "${WDIR}" != "${BDIR}" && "${WDIR}" != "." ) then
		rm -r ${WDIR}
	endif
else
	echo "no input file specified." > /dev/stderr
	echo "run with '-h' for help." > /dev/stderr
endif

if ( $globalstatus != 0 ) then
	echo "script finished with status ${globalstatus}." > /dev/stderr
endif
