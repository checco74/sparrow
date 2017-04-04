/* toolz.h functions definition */
#include "toolz.h"

ofstream testDatalogfile;





void SampleType::init(int L)
{
	this->length = L;
	this->key = L/2;
}





resultsTable::resultsTable()
{
	this->cnt = 0;
	this->mycnt = 0;
	this->good = 0;
	this->miss = 0;
	this->confusionVector = NULL;
	this->confusionVector = new unsigned long*[NUM_OF_AMINO_ACIDS+1];
	for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++) this->confusionVector[aa] = NULL;
}



resultsTable::~resultsTable()
{
	for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
		if (this->confusionVector[aa] != NULL) delete[] confusionVector[aa];
	delete[] confusionVector;
}


void resultsTable::init(int classesCount)
{
	for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
	{
		this->confusionVector[aa] = new unsigned long[classesCount];
		for (int p = 0; p<classesCount; p++) this->confusionVector[aa][p] = 0;
	}
}




void eraseStuff()
{
	for (int p = 0; p<patternCount; p++)
	{
		if (neighboursCount[p]!=NULL) free(neighboursCount[p]);
		if (auxnneighboursCount[p]!=NULL) free(auxnneighboursCount[p]);
		if (auxaneighboursCount[p]!=NULL) free(auxaneighboursCount[p]);
	}
	if (neighboursCount!=NULL) free(neighboursCount);
	if (auxnneighboursCount!=NULL) free(auxnneighboursCount);
	if (auxaneighboursCount!=NULL) free(auxaneighboursCount);
	if (minlength!=NULL) free(minlength);
	if (powFactor!=NULL) free(powFactor);
	if (shift!=NULL) free(shift);
}




void checkDirectory()
{
	FILE *fp;
	if ((fp = fopen(".prm/.chk", "ab"))==NULL)
	{
		cout<<"\n.prm/ directory missing"<<endl;
		exit(1);
	}
	else fclose(fp);
}




int getStructCount()
{
	switch(detail)
	{
	case FULL: return totalStructCount;
	case ANY_COIL: return totalStructCount-1;
	case NO_BENDS: return totalStructCount-2;
	case NO_TURNS: return totalStructCount-3;
	case NO_REGIONS: return totalStructCount-8;
	case NO_REGIONS_ANY_COIL: return totalStructCount-9;
	case NO_REGIONS_NO_BENDS: return 11;
	case NO_REGIONS_NO_TURNS: return 10;
	case DSSP_PLUS_SHORTS: return 9;
	case DSSP_STANDARD: return 8;
	case DSSP_NO_BENDS: return 7;
	case DSSP_NO_TURNS: return 6;
	case DSSP_NO_EXOTIC: return 5;
	case Q3_PLUS_TURNS: return 4;
	default: return NUM_OF_MOTIVES;
	}
}




int project(StructureType p, int definitionLevel)
{
	detailLevel local_definitionLevel = DSSP_STANDARD;
	if (definitionLevel == -1) local_definitionLevel = detail;
	switch(local_definitionLevel)
	{
	case FULL: return p;
	case ANY_COIL:
	{
	if (p!=SHORT_COIL) return p;
	else return RANDOM_COIL;
	}
	case NO_BENDS:
	{
	if (p!=SHORT_COIL && p!=BEND) return p;
	else return RANDOM_COIL;
	}
	case NO_TURNS:
	{
	if (p!=SHORT_COIL && p!=BEND && p!=TURN) return p;
	else return RANDOM_COIL;
	}
	case NO_REGIONS:
	case NO_REGIONS_ANY_COIL:
	case NO_REGIONS_NO_BENDS:
	case NO_REGIONS_NO_TURNS:
	{
	switch(p)
	{
	case SHORT_HELIX :
	return SHORT_ALPHA_HELIX;
	case START_HELIX :
	case CENTER_HELIX :
	case END_HELIX :
	return ALPHA_HELIX;
	case SHORT_T_HELIX :
	return SHORT_THIN_HELIX;
	case START_T_HELIX :
	case CENTER_T_HELIX :
	case END_T_HELIX :
	return THIN_HELIX;
	case SHORT_P_HELIX :
	return SHORT_FAT_HELIX;
	case START_P_HELIX :
	case CENTER_P_HELIX :
	case END_P_HELIX :
	return FAT_HELIX;
	case SHORT_STRAND :
	return SHORT_BETA_STRAND;
	case START_STRAND :
	case CENTER_STRAND :
	case END_STRAND :
	return BETA_STRAND;
	case ISOLATED_STRAND :
	return BETA_BRIDGE;
	case TURN :
	if (detail!=NO_REGIONS_NO_TURNS) return TURN_LOOP;
	else return GENERIC_COIL;
	case BEND :
	if (detail==NO_REGIONS or detail==NO_REGIONS_ANY_COIL) return BEND_LOOP;
	else return GENERIC_COIL;
	case SHORT_COIL :
	if (detail==NO_REGIONS) return SHORT_GENERIC_COIL;
	else return GENERIC_COIL;
	default :
	return GENERIC_COIL;
	}
	}
	case DSSP_PLUS_SHORTS:
	case DSSP_STANDARD:
	case DSSP_NO_BENDS:
	case DSSP_NO_TURNS:
	{
	switch(p)
	{
	case SHORT_HELIX :
	case START_HELIX :
	case CENTER_HELIX :
	case END_HELIX :
	return DSSP_HELIX;
	case SHORT_T_HELIX :
	case START_T_HELIX :
	case CENTER_T_HELIX :
	case END_T_HELIX :
	return DSSP_T_HELIX;
	case SHORT_P_HELIX :
	case START_P_HELIX :
	case CENTER_P_HELIX :
	case END_P_HELIX :
	return DSSP_P_HELIX;
	case SHORT_STRAND :
	case START_STRAND :
	case CENTER_STRAND :
	case END_STRAND :
	return DSSP_STRAND;
	case ISOLATED_STRAND :
	return DSSP_BRIDGE;
	case TURN :
	if (detail!=DSSP_NO_TURNS) return DSSP_TURN;
	else return DSSP_COIL;
	case BEND :
	if (detail!=DSSP_NO_BENDS and detail!=DSSP_NO_TURNS) return DSSP_BEND;
	else return DSSP_COIL;
	case SHORT_COIL :
	if (detail==DSSP_PLUS_SHORTS) return DSSP_SHORT;
	else return DSSP_COIL;
	default :
	return DSSP_COIL;
	}
	}
	case DSSP_NO_EXOTIC:
	case Q3_PLUS_TURNS:
	default:
	{
		if (encoding=="strict")
		{
			switch(p)
			{
			case SHORT_HELIX :
			case START_HELIX :
			case CENTER_HELIX :
			case END_HELIX :
				return HELIX;
			case SHORT_STRAND :
			case START_STRAND :
			case CENTER_STRAND :
			case END_STRAND :
				return STRAND;
			case TURN :
				if (detail!=Q3_STANDARD) return _TURN;
				else return COIL;
			case BEND :
				if (detail==DSSP_NO_EXOTIC) return _BEND;
				else return COIL;
			default :
				return COIL;
			}
		}
		else if (encoding=="loose")
		{
			switch(p)
			{
			case SHORT_HELIX :
			case START_HELIX :
			case CENTER_HELIX :
			case END_HELIX :
			case SHORT_T_HELIX :
			case START_T_HELIX :
			case CENTER_T_HELIX :
			case END_T_HELIX :
			case SHORT_P_HELIX :
			case START_P_HELIX :
			case CENTER_P_HELIX :
			case END_P_HELIX :
				return HELIX;
			case SHORT_STRAND :
			case START_STRAND :
			case CENTER_STRAND :
			case END_STRAND :
			case ISOLATED_STRAND :
				return STRAND;
			case TURN :
				if (detail!=Q3_STANDARD) return _TURN;
				else return COIL;
			case BEND :
				if (detail==DSSP_NO_EXOTIC) return _BEND;
				else return COIL;
			default :
				return COIL;
			}
		}
		else
		{
			cout<<"\nError: unknown encoding: "<<encoding<<endl;
			exit(1);
		}
	}
	}
}




int coil()
{
	if (detail<NO_REGIONS) return RANDOM_COIL;
	else if (detail<DSSP_PLUS_SHORTS) return GENERIC_COIL;
	else if (detail<DSSP_NO_EXOTIC) return DSSP_COIL;
	else return _COIL;
}




Category getQ3(int p, detailLevel dlevel)
{
	if (dlevel<NO_REGIONS)
	{
		if (encoding=="strict")
		{
		switch(p)
		{
			case SHORT_HELIX:
			case START_HELIX:
			case CENTER_HELIX:
			case END_HELIX:
			return HELIX;
			case SHORT_STRAND:
			case START_STRAND:
			case CENTER_STRAND:
			case END_STRAND:
			return STRAND;
			default:
			return COIL;
		}
		}
		else if (encoding=="loose")
		{
		switch(p)
		{
			case SHORT_HELIX:
			case START_HELIX:
			case CENTER_HELIX:
			case END_HELIX:
			case SHORT_T_HELIX:
			case START_T_HELIX:
			case CENTER_T_HELIX:
			case END_T_HELIX:
			case SHORT_P_HELIX:
			case START_P_HELIX:
			case CENTER_P_HELIX:
			case END_P_HELIX:
			return HELIX;
			case SHORT_STRAND:
			case START_STRAND:
			case CENTER_STRAND:
			case END_STRAND:
			case ISOLATED_STRAND:
			return STRAND;
			default:
			return COIL;
		}
	}
	else
	{
		cout<<"\nError: unknown encoding: "<<encoding<<endl;
		exit(1);
	}
	}
	else if (dlevel<DSSP_PLUS_SHORTS)
	{
		if (encoding=="strict")
		{
		switch(p)
		{
			case SHORT_ALPHA_HELIX:
			case ALPHA_HELIX:
			return HELIX;
			case SHORT_BETA_STRAND:
			case BETA_STRAND:
			return STRAND;
			default:
			return COIL;
		}
	}
	else if (encoding=="loose")
	{
		switch(p)
		{
			case SHORT_ALPHA_HELIX:
			case ALPHA_HELIX:
			case SHORT_THIN_HELIX:
			case THIN_HELIX:
			case SHORT_FAT_HELIX:
			case FAT_HELIX:
			return HELIX;
			case SHORT_BETA_STRAND:
			case BETA_STRAND:
			case BETA_BRIDGE:
			return STRAND;
			default:
			return COIL;
		}
	}
	else
	{
		cout<<"\nError: unknown encoding: "<<encoding<<endl;
		exit(1);
	}
	}
	else if (dlevel<DSSP_NO_EXOTIC)
	{
		if (encoding=="strict")
		{
		switch(p)
		{
			case DSSP_HELIX:
			return HELIX;
			case DSSP_STRAND:
			return STRAND;
			default:
			return COIL;
		}
	}
	else if (encoding=="loose")
	{
		switch(p)
		{
			case DSSP_HELIX:
			case DSSP_T_HELIX:
			case DSSP_P_HELIX:
			return HELIX;
			case DSSP_STRAND:
			case DSSP_BRIDGE:
			return STRAND;
			default:
			return COIL;
		}
	}
	else
	{
		cout<<"\nError: unknown encoding: "<<encoding<<endl;
		exit(1);
	}
	}
	else
	{
	switch(p)
	{
	case HELIX:
		return HELIX;
	case STRAND:
		return STRAND;
	default:
		return COIL;
	}
	}
}




char sstrucQ3(char z)
{
	if (encoding=="strict")
	{
		if (z=='H') return 'H';
		else if (z=='e' || z=='E') return 'e';
		else return '.';
	}
	else if (encoding=="loose")
	{
		if (z=='H' || z=='G' || z=='I') return 'H';
		else if (z=='e' || z=='E' || z=='b' || z=='B') return 'e';
		else return '.';
	}
	else
	{
		cout<<"\nError: unknown encoding: "<<encoding<<endl;
		exit(1);
	}
}




char sstruc2char(int p, int definitionLevel)
{
	detailLevel local_definitionLevel = DSSP_STANDARD;
	if (definitionLevel == -1) local_definitionLevel = detail;
	if (p==StructCount) return '?';
	if (local_definitionLevel<NO_REGIONS)
	{
	switch(p)
	{
	case SHORT_HELIX: return 'h';
	case START_HELIX: return 'h';
	case CENTER_HELIX: return 'H';
	case END_HELIX: return 'h';
	case SHORT_T_HELIX: return 'g';
	case START_T_HELIX: return 'g';
	case CENTER_T_HELIX: return 'G';
	case END_T_HELIX: return 'g';
	case SHORT_P_HELIX: return 'i';
	case START_P_HELIX: return 'i';
	case CENTER_P_HELIX: return 'I';
	case END_P_HELIX: return 'i';
	case SHORT_STRAND: return 'e';
	case START_STRAND: return 'e';
	case CENTER_STRAND: return 'E';
	case END_STRAND: return 'e';
	case ISOLATED_STRAND: return 'B';
	case SHORT_COIL: return 'c';
	case TURN: return 'T';
	case BEND: return 'S';
	default: return '.';
	}
	}
	else if (local_definitionLevel<DSSP_PLUS_SHORTS)
	{
	switch(p)
	{
	case ALPHA_HELIX: return 'H';	// long helix residue
	case SHORT_ALPHA_HELIX: return 'h';	// short helix residue
	case BETA_STRAND: return 'E';	// long strand residue
	case SHORT_BETA_STRAND: return 'e';	// short strand residue
	case BETA_BRIDGE: return 'B';	// residue in an isolated strand
	case THIN_HELIX: return 'G';	// long 3-helix residue
	case SHORT_THIN_HELIX: return 'g';	// short 3-helix residue
	case FAT_HELIX: return 'I';	// long 5-helix residue
	case SHORT_FAT_HELIX: return 'i';	// short 5-helix residue
	case TURN_LOOP: return 'T';	// residue on an hydrogen-bonded turn
	case BEND_LOOP: return 'S';	// residue on a bend
	case GENERIC_COIL: return 'c';	// residue in a long coil
	case SHORT_GENERIC_COIL: return '.';	// residue in a short coil
	default: return '.';
	}
	}
	else if (local_definitionLevel<DSSP_NO_EXOTIC)
	{
	switch(p)
	{
	case DSSP_HELIX: return 'H';
	case DSSP_STRAND: return 'E';
	case DSSP_T_HELIX: return 'G';
	case DSSP_P_HELIX: return 'I';
	case DSSP_BRIDGE: return 'B';
	case DSSP_TURN: return 'T';
	case DSSP_BEND: return 'S';
	default: return '.';
	}
	}
	else
	{
	switch(p)
	{
	case HELIX: return 'H';
	case STRAND: return 'E';
	case _TURN: return 'T';
	case _BEND: return 'S';
	default: return '.';
	}
	}
}




int zcode(char c)
{
	switch(c)
	{
	case 'C':
	case 'c': return C_TERM;
	case 'N':
	case 'n': return N_TERM;
	case 'H':
	case 'h': return CENTER;
	default:
		return zoneCount;
	}
}




char zname(int z)
{
	switch(z)
	{
	case 0: return 'C';
	case 1: return 'N';
	case 2: return 'H';
	default: return 'A';
	}
}




int init_pattern(int argc, char *argv[], int (&zone), int &antipattern)
{
	int pattern = -1;
	cout<<"defining secondary structure motifs...."<<endl;
	for (int i = 0; i<argc; i++)
	{
		if (strncmp(argv[i], "-level=", 7)==0)
		{
			int temp;
			sscanf(argv[i], "-level=%d", &temp);
			if (temp<=Q3_STANDARD) detail = (detailLevel)temp;
			else cout<<"Warning: invalid level of detail: "<<temp<<endl<<"setting default...."<<endl;
		}
	}
	StructCount = getStructCount();
	minlength = (int*)malloc(StructCount*sizeof(int));
	for (int i = 0; i<StructCount; i++) minlength[i] = 1;
	for (int i = 0; i<argc; i++)
	{
		sscanf(argv[i], "ss%d", &subdim);
		zone = zcode('A');
	}
	if (subdim<0)
	{
		cout<<"\nError: invalid pattern dimension: "<<subdim<<endl;
		exit(1);
	}
	patternCount = (int)pow((double)StructCount, (double)subdim);
	pairsCount = patternCount+(int)binomialCoefficient(patternCount, 2);
	if ((minlength = (int*)realloc(minlength, patternCount*sizeof(int)))==NULL
		|| (super_reweight = (double*)malloc(patternCount*sizeof(double)))==NULL
		|| (neighboursCount = (int**)malloc(patternCount*sizeof(int*)))==NULL
		|| (auxnneighboursCount = (int**)malloc(patternCount*sizeof(int*)))==NULL
		|| (auxaneighboursCount = (int**)malloc(patternCount*sizeof(int*)))==NULL)
	{
		cerr<<"\nError: memory allocation error (SSP1)."<<endl;
		exit(1);
	}
	for (int i = 0; i<patternCount; i++)
	{
		neighboursCount[i] = NULL;
		auxnneighboursCount[i] = NULL;
		auxaneighboursCount[i] = NULL;
		if ((neighboursCount[i] = (int*)malloc(pairsCount*sizeof(int)))==NULL
			|| (auxnneighboursCount[i] = (int*)malloc(pairsCount*sizeof(int)))==NULL
			|| (auxaneighboursCount[i] = (int*)malloc(pairsCount*sizeof(int)))==NULL)
		{
			cerr<<"\nError: memory allocation error (SSP2)."<<endl;
			exit(1);
		}
	}
	for (int i = 0; i<patternCount; i++)
	{
		minlength[i] = 1;
		super_reweight[i] = 1./(double)patternCount;
		for (int j = 0; j<pairsCount; j++)
		{
			neighboursCount[i][j] = SPARROW_NEIGHBOURS;
			auxnneighboursCount[i][j] = 0;
			auxaneighboursCount[i][j] = 0;
		}
	}
	return pattern;
}




int getMotifsPair(int c0, int c1)
{
	int ID = -1;
	if (c0<0 || c0>=patternCount)
	{
		cout<<"\ninvalid motif: p0="<<c0<<endl;
		exit(1);
	}
	else
	{
		if (c1<0 || c1>=patternCount || c0==c1) ID = c0;
		else if (c0<c1) ID = (c0+1)*patternCount - ((c0+2)*(c0+1))/2 + c1;
		else if (c0>c1) ID = (c1+1)*patternCount - ((c1+2)*(c1+1))/2 + c0;
	}
	return ID;
}




void setShifts()
{
		int D = 0;
		shift = (int*)malloc(subdim*sizeof(int));
		for (int d = 0; d<subdim; d++)
		{
			D = ((2*d-subdim+1)*(gapsize+1));
			shift[d] = D%2<0 ? D%2 : 0 + D/2;
		}
		for (int d = 0; d<subdim; d++)
			if (D%2!=0) shift[d]++;
}




void setPowerBaseUnits()
{
	powFactor = (int*)malloc(subdim*sizeof(int));
	for (int d = 0; d<subdim; d++)
		powFactor[d] = (int)pow((double)StructCount,(double)(subdim-d-1));
}




int buildFilename(char *name, char *ext, scoreMachine *sm, int d, int z, int dir, readMode opt)
{
	FILE *fp;
	char ap[7];
	if (sm->Antipattern==-1) sprintf(ap, "p%04d", sm->Pattern);
	else sprintf(ap, "p%04dvs%04d", sm->Pattern, sm->Antipattern);
	if (dir==1) sprintf(name, "%s/%slev%02dsize%dat%02ul%02u-%cdx_nctg", ext, ap, detail, d, sm->winkey, sm->winlen, zname(z));
	else sprintf(name, "%s/%slev%02dsize%dat%02ul%02u-%csx_nctg", ext, ap, detail, d, sm->winkey, sm->winlen, zname(z));
	if (!sm->power)
	{
	if (strlen(name)+8<111) strcat(name, ".lnr");
	else {cout<<"\ndamn! filename exceeds allowed length."<<endl; exit(1);}
	}
	if (sm->linsave)
	{
	if (strlen(name)+8<111) strcat(name, ".ctg");
	else {cout<<"\ndamn! filename exceeds allowed length."<<endl; exit(1);}
	}
	if (opt==AUXN_MODE)
	{
	if (strlen(name)+strlen(".auxn")<111) strcat(name, ".auxn");
	else {cout<<"\ndamn! filename exceeds allowed length."<<endl; exit(1);}
	}
	else if (opt==AUXA_MODE)
	{
	if (strlen(name)+strlen(".auxa")<111) strcat(name, ".auxa");
	else {cout<<"\ndamn! filename exceeds allowed length."<<endl; exit(1);}
	}
	if (strlen(name)+strlen(ext)<111) strcat(name, ext);
	else {cout<<"\ndamn! filename exceeds allowed length."<<endl; exit(1);}
	if ((fp=fopen(name, "rb"))!=NULL)
	{
	fclose(fp);
	return 0;
	}
	else return 1;
}




void writeMatrix(char *filename, double *Nw, uint size)
{
	uint ic = 0, ij = 0;
	FILE *fp = fopen(filename, "wb");
	for (uint j = 0; j<size; j++)
	{
	ic += j;
	ij = ic+j;
	for (uint i = j; i<size; i++)
	{
	fprintf(fp, "%19.9lf", Nw[ij]);
	fflush(fp);
	ij += i+1;
	}
	fprintf(fp, "\n");
	fflush(fp);
	}
	fclose(fp);
}




void writeMatrix(char *filename, double **Nw, uint size)
{
	FILE *fp = fopen(filename, "wb");
	for (uint i = 0; i<size; i++)
	{
	for (uint j = 0; j<=i; j++)
	{
	fprintf(fp, "%19.12lf", Nw[i][j]);
	fflush(fp);
	}
	fprintf(fp, "\n");
	fflush(fp);
	}
	fclose(fp);
}




void writeMatrix(char *filename, int **Nw, uint size)
{
	FILE *fp = fopen(filename, "wb");
	for (uint i = 0; i<size; i++)
	{
	for (uint j = 0; j<size; j++)
	{
	fprintf(fp, "%9d", Nw[i][j]);
	fflush(fp);
	}
	fprintf(fp, "\n");
	fflush(fp);
	}
	fclose(fp);
}




void writeStats(char *filename, double *x, uint size)
{
	FILE *fp = fopen(filename, "ab");
	for (uint i = 0; i<size; i++)
		fprintf (fp, "%19.16lf\n", x[i]);
	fflush(fp);
	fclose(fp);
}




void storeMatrix(char *matrixfilename, double **Nw, unsigned int size, unsigned int &N)
{
	FILE *fp;
	if (verbosity>0) cout<<"storing stats matrix \""<<matrixfilename<<"\"...."<<endl;
	if ((fp = fopen64(matrixfilename, "wb")) == NULL)
	{
	fprintf(stderr, "Can't open file %s for writing\n", matrixfilename);
	printf("exit h7(%d)\n", 14);
	exit(14);
	}
	for (uint k = 0; k<size; k++)
	{
	for (uint h = k; h<size; h++)
	{
	double entry = Nw[h][k]; // column-wise ordering
	fwrite(&entry, sizeof(double), 1, fp);
	}
	}
	fwrite(&N, sizeof(unsigned int), 1, fp);
	fclose(fp);
}




void storeStats(char *filename, double *x, unsigned int size, unsigned int &N)
{
	FILE *fp;
	if (verbosity>0)	cout<<"storing stats \""<<filename<<"\"...."<<endl;
	if ((fp = fopen(filename, "wb")) == NULL)
	{
	fprintf(stderr, "Can't open file %s for writing\n", filename);
	printf("exit j7(%d)\n", 14);
	exit(14);
	}
	else
	{
	fwrite(x, sizeof(double), size, fp);
	fwrite(&N, sizeof(unsigned int), 1, fp);
	}
	fclose(fp);
}




void storeSolutions(const char *filename, double b, double *w, unsigned int size)
{
	FILE *fp;
	if (verbosity>0)	cout<<"storing solutions \""<<filename<<"\"...."<<endl;
	if ((fp = fopen(filename, "wb")) == NULL)
	{
	fprintf(stderr, "Can't open file %s for writing\n", filename);
	printf("exit j7(%d)\n",14);
	exit(14);
	}
	fwrite(w, sizeof(double), size, fp);
	fwrite(&b, sizeof(double), 1, fp);
	fclose(fp);
}




int getMatrix(char *name, double **(&mtx), uint size, unsigned int *N)
{
	FILE *fs;
	uint cnt = 0;
	if ((fs = fopen(name,"rb"))==NULL)
	{
	cout<<"\nerror opening file \""<<name<<"\"";
	cout<<"\nAborting...."<<endl;
	exit(1);
	}
	for (uint k = 0; k<size; k++)
	{
	for (uint h = k; h<size; h++)
	{
		if (fread(&mtx[h][k], sizeof(double), 1, fs)!=1)
		{
			cout<<"\ndata reading error.";
			return (int)cnt;
		}
		else cnt++;
	}
	}
	if (N!=NULL && fread(N, sizeof(unsigned int), 1, fs)!=1)
		cout<<"\nerror - no sequence count found: "<<endl;
	fclose(fs);
	return (int)cnt;
}




int getStats(char *name, double *(&x), uint size, unsigned int *N)
{
	FILE *fs;
	uint cnt = 0;
	if ((fs = fopen(name,"rb"))==NULL)
	{
	cout<<"\nerror opening file \""<<name<<"\"";
	cout<<"\nAborting...."<<endl;
	exit(1);
	}
	if ((cnt = fread(x, sizeof(double), size, fs))!=size)
	{
	cout<<"\nerror - incomplete data: ";
	cout<<cnt<<"/"<<size<<endl;
	}
	if (N!=NULL && fread(N, sizeof(unsigned int), 1, fs)!=1)
		cout<<"\nerror - no sequence count found: "<<endl;
	fclose(fs);
	return (int)cnt;
}




int getSolutions(const char *name, double &b, double *(&w), uint size)
{
	FILE *fs;
	uint cnt = 0;
	if ((fs = fopen(name, "rb"))==NULL)
	{
		if ( verbosity > 0 )
			cout<<"Warning: file \""<<name<<"\" not found."<<endl;
		return -1;
	}
	if ((cnt = fread(w, sizeof(double), size, fs))!=size || fread(&b, sizeof(double), 1, fs)!=1)
	{
		cout<<"Warning: incomplete data in file \""<<name<<"\": "<<cnt<<"/"<<size<<endl;
		return -1;
	}
	fclose(fs);
	return (int)cnt;
}




int getSVMSolutions(char *name, double *bterm, double *Wsvm, uint size)
{
	char ch = '?';
	char dummy[111];
	double alphay = 0.;
	uint indx = 0, d = 0;
	uint tot = 0, cnt = 0;
	FILE *fp;
	fp = fopen(name,"rb");
	if (fp==NULL)
	{
		perror("aargh!");
		exit(7);
	}
	fgets(dummy, 110, fp); // removing headers and useless stuff
	for (int i = 1; i<10; i++)
	{
	fscanf(fp, "%u", &tot);
	fgets(dummy, 110, fp);
	}
	fscanf(fp, "%lf", &(*bterm));
	*bterm = -(*bterm);
	fgets(dummy, 110, fp);
	cout<<"number of support vectors: "<<tot-1<<endl;
	cout<<"bterm: "<<*bterm<<endl;
	Wsvm = memresVec((long unsigned int)size);
	while (cnt<tot-1)
	{
	fscanf(fp, "%lf", &alphay);
	while (ch!='#')
	{
	if (fscanf(fp, "%u:%u", &indx, &d)>0)
	Wsvm[indx-1] += alphay;
	ch = getc(fp);
	}
	cnt++;
	ch='?';
	}
	fclose(fp);
	return tot-1;
}




int segmentSize(const char *(&s), int k)
{
	int sslength = 0;
	for (int h = k; h<(int)strlen(s) && sstrucQ3(s[h])==sstrucQ3(s[k]); h++) sslength++;
	for (int h = k-1; h>=0 && sstrucQ3(s[h])==sstrucQ3(s[k]); h--) sslength++;
	return sslength;
}




int segmentSize(int *(&s), int protlen, int k)
{
	int sslength = 0;
	for (int h = k; h<protlen && s[h]==s[k]; h++) sslength++;
	for (int h = k-1; h>=0 && s[h]==s[k]; h--) sslength++;
	return sslength;
}




int computeEffectiveNS(int pattern, int *NS, int **neighbours)
{
	int N = 0;
	for (int ID = 0; ID<pairsCount; ID++)
		N += (2*neighbours[pattern][ID]+1)*NS[ID];
	return N;
}




int computeNQuadratic(int totalNS)
{
	int N = 0;
	N = totalNS+(totalNS*(totalNS+1))/2 + 1;
	return N;
}




int countTypes(scoreMachine ***(&sm), int **NS, SampleType *(&ST))
{
	// there'll always be at least one Type:
	// the Type=0 is used for initializations.
	int count = 1;
	bool match = false;
	ST = (SampleType*)malloc(sizeof(SampleType));
	ST[0].init();
	unsigned long memory = (unsigned long)setSize*(2*sizeof(uint)+(NUM_OF_AMINO_ACIDS+3)*(unsigned long)ST[0].length*sizeof(char));
	for (int prnID = 0; prnID<patternCount; prnID++) {
		for (int pairID = 0; pairID<pairsCount; pairID++) {
			for (int s = 0; s<NS[prnID][pairID]; s++)
			{
				match = false;
				for (int i = 0; i<count; i++)
					if (sm[prnID][pairID][s].winlen==ST[i].length && sm[prnID][pairID][s].winkey==ST[i].key)
						match = true;
				if (!match)
				{
					ST = (SampleType*)realloc(ST, sizeof(SampleType)*(count+1));
					if (ST==NULL)
					{
						cout<<"\n>memory reallocation error."<<endl;
						exit(1);
					}
					ST[count].init();
					ST[count].key = sm[prnID][pairID][s].winkey;
					ST[count].length = sm[prnID][pairID][s].winlen;
					memory += (unsigned long)setSize*(2*sizeof(uint)+(NUM_OF_AMINO_ACIDS+3)*(unsigned long)ST[count].length*sizeof(char));
					count++;
				}
			}
		}
	}
	if (verbosity>0)
	{
		cout<<"\ntotal number of sample types: "<<count;
		cout<<"\nmemory required: "<<memory/1024<<" kB."<<endl;
	}
	return count;
}




string matchProfileFormat(string line)
{
	bool gotit = false;
	ostringstream rss;
	size_t chposNA = line.find_first_of(commentchars);
	size_t chposA = line.find_first_not_of(commentchars);
	if (chposA!=string::npos && (chposNA==string::npos || chposA<chposNA))
	{
		uint index = 0;
		char aachar = '#';
		stringstream ss(line);
		while (ss && !gotit)
		{
			string backupbuffer = ss.str();
			// the line is read a character at a time;
			// it will be considered good if the function succeeds in reading a profile from some character on;
			// the backupbuffer serves the scope of storing the substring of characters yet to be read.
// 			cout<<"0. backup buffer is: "<<backupbuffer<<endl;
			if (ss>>aachar)
			{
				index = ss.tellg();
				if (ss.str().size()>index) // inside the string
				{
					int k = 0;
					int entry = 0;
					backupbuffer = ss.str().substr(index);
// 					cout<<"1. backup buffer is: "<<backupbuffer<<endl;
					if (isalpha(aachar)) rss<<aachar<<' ';
					else
					{
						ss.putback(aachar);
						rss<<"! ";
					}
					while (ss>>entry && k<NUM_OF_AMINO_ACIDS)
					{
						rss<<entry<<' ';
						index = ss.tellg();
						if (ss.str().size()>index) // inside the string
						{
							backupbuffer = ss.str().substr(index);
// 							cout<<"3. backup buffer is: "<<backupbuffer<<endl;
						}
						k++;
					}
					if (k<NUM_OF_AMINO_ACIDS)
					{
						rss.str("");
						ss.str(backupbuffer);
						ss.clear();
					}
					else gotit = true;
				}
			}
		}
	}
	return rss.str();
}




double *getProfile(string &rawProfile, int profLength, int dir, int pos)
{
	int m, n;
	int profStart = pos*NUM_OF_AMINO_ACIDS;
	int actualLength = rawProfile.length();
	double *prof = NULL;
	prof = new double[profLength];
	if (verbosity>2) cout<<"profile "<<profStart<<": "<<endl;
	for (m = 0; m<profLength; m++)
	{
		if (dir==1) n = m;
		else n = profLength-1-m;
		char rawEntry = 'E';
		int profIndex = profStart+n;
		if (profIndex>=0 && profIndex<actualLength) rawEntry = rawProfile[profIndex]; // letter code from the DBInterface
		prof[m] = getProfileEntry(rawEntry);
		if (verbosity>2)
		{
			if (m%NUM_OF_AMINO_ACIDS==0) cout<<endl;
			cout<<rawEntry<<flush;
			if (verbosity>3) cout<<'['<<prof[m]<<']'<<flush;
		}
	}
	if (verbosity>2) cout<<endl;
	return prof;
}




double *getCollapsedProfile(string &sequence, int length, int winLength, int profLength, int dir, int pos)
{
	int mp, np;
	int offset = 0;
	double *prof = NULL;
	prof = new double[profLength];
	for (int m = 0; m<profLength; m++) prof[m] = 0.;
	for (mp = 0; mp<winLength; mp++)
	{
		if (dir==1) np = mp;
		else np = winLength-1-mp;
		int profIndex = 0;
		int Index = pos+np;
		if (Index>=0 && Index<length)
		{
			int aaID = aminoacid::astype(sequence[Index]);
			if (dir==-1) aaID = NUM_OF_AMINO_ACIDS-1-aaID;
			profIndex = offset+aaID;
			if (aaID>=0 && aaID<NUM_OF_AMINO_ACIDS && profIndex>=0 && profIndex<profLength);
				prof[profIndex] = 1.;
		}
		else
		{
			for (int aaID = 0; aaID<NUM_OF_AMINO_ACIDS; aaID++)
				prof[offset+aaID] = getProfileEntry('E');
		}
		offset += NUM_OF_AMINO_ACIDS;
	}
	return prof;
}




double getProfileEntry(char rawEntry)
{
	double entry = 0.;
	// the code is converted into an integer and traslated to be used as array index
	int index = TRdecode(rawEntry) + aminoacid::maxValue/2;
	// the array index is used to get float values in (0,1)
	entry = aminoacid::logValue[index];
	return entry;
}




int getNCPartial(int prnID, int *mymotif, int k)
{
	int mrk = prnID;
	for (int d = 0; d<subdim-1 && mrk>=0; d++) // doesn't enter here if subdim==1 of course
	{
		int partmotif = 0;
		int pk = (k-(subdim-1-d)*(gapsize+1));
		if (pk>=0) partmotif = mymotif[pk];
		else partmotif = coil();
		if (mrk/powFactor[d]==partmotif)
		mrk = mrk%powFactor[d];
		else mrk = -1;
	}
	return mrk;
}




int getCNPartial(int prnID, int *mymotif, int h, int protlen)
{
	int mrk = prnID;
	for (int d = 0; d<subdim-1 && mrk>=0; d++) // doesn't enter here if subdim==1 of course
	{
		int partmotif = 0;
		int ph = (h+(subdim-1-d)*(gapsize+1));
		if (ph<protlen) partmotif = mymotif[ph];
		else partmotif = coil();
		if (mrk/powFactor[d]==partmotif)
		mrk = mrk%powFactor[d];
		else mrk = -1;
	}
	return mrk;
}




int getStartingPosition(scoreMachine &sm, int centre, int dir)
{
	if (dir==1) return centre-(int)sm.winkey;
	else return centre-(int)sm.winlen+(int)sm.winkey+1;
}




void initResidue(int k, double **(&score), int *(&motif))
{
	motif[k] = coil();
	for (int p = 0; p<StructCount; p++)
	score[k][p] = -11.;
}




int scoreResidue(int k, int dir, Protein *chain, int protlen, scoreGear *(&SG), double **(&myscore), int *(&mymotif),
	int amode, readMode rmode)
{
	int mrk = NUM_OF_MOTIVES;
	int count = 0;
	int c = 0; // is the centre of the secondary structure pattern
	double maxscore = -11.;
	initResidue(k, myscore, mymotif);
	for (int prn = 0; prn<patternCount; prn++)
	{
		int flag = 1; // tells whether the current pattern "fits in" ---> its score is to be computed (multidimensional patterns)
		if (amode==NO_TRANSITIONS)
		{
			c = k;
			mrk = prn/((patternCount-1)/(StructCount-1));
			flag = prn%((patternCount-1)/(StructCount-1));
		}
		else
		{
			if (dir==1)
			{
				c = k-shift[subdim-1];
				mrk = getNCPartial(prn, mymotif, k);
			}
			else
			{
				c = k+shift[subdim-1];
				mrk = getCNPartial(prn, mymotif, k, protlen);
			}
			if (mrk>=0) flag = 0;
			else flag = 1;
		}
		if (flag==0)
		{
			double score = 0.;
			double *SF = NULL;
			if (verbosity>1) cout<<"pattern "<<prn;
			int zone = (int)CHAIN_CENTER;
			switch(rmode)
			{
				case STD_MODE :
				{
					if (verbosity>1) cout<<" standard score: ";
					SF = getScoreFactors(prn, SG[prn].NSSC, SG[prn].partNS, SG[prn].totalNS, SG[prn].sm,
						neighboursCount, chain->sequence, chain->profile, protlen, zone, c, dir);
					for (int s = 0; s<SG[prn].NSSC; s++)
						score += SF[s]*SG[prn].ssc.first[zone][s];
					break;
				}
				case AUXN_MODE :
				{
					if (SG[prn].auxnNSSC>1)
					{
						if (verbosity>1) cout<<" auxiliary-n score: ";
						SF = getScoreFactors(prn, SG[prn].auxnNSSC, SG[prn].partauxnNS, SG[prn].totalauxnNS, SG[prn].auxn_sm,
											 auxnneighboursCount, chain->sequence, chain->profile, protlen, zone, c, dir);
						for (int s = 0; s<SG[prn].auxnNSSC; s++)
							score += SF[s]*SG[prn].auxn_ssc.first[zone][s];
					}
					else
					{
						if (verbosity>1) cout<<" standard score: ";
						SF = getScoreFactors(prn, SG[prn].NSSC, SG[prn].partNS, SG[prn].totalNS, SG[prn].sm,
											 neighboursCount, chain->sequence, chain->profile, protlen, zone, c, dir);
						for (int s = 0; s<SG[prn].NSSC; s++)
							score += SF[s]*SG[prn].ssc.first[zone][s];
					}
					break;
				}
				case AUXA_MODE :
				{
					if (SG[prn].auxaNSSC>1)
					{
						if (verbosity>1) cout<<" auxiliary-a score: ";
						SF = getScoreFactors(prn, SG[prn].auxaNSSC, SG[prn].partauxaNS, SG[prn].totalauxaNS, SG[prn].auxa_sm,
											 auxaneighboursCount, chain->sequence, chain->profile, protlen, zone, c, dir);
						for (int s = 0; s<SG[prn].auxaNSSC; s++)
							score += SF[s]*SG[prn].auxa_ssc.first[zone][s];
					}
					else
					{
						if (verbosity>1) cout<<" standard score: ";
						SF = getScoreFactors(prn, SG[prn].NSSC, SG[prn].partNS, SG[prn].totalNS, SG[prn].sm,
											 neighboursCount, chain->sequence, chain->profile, protlen, zone, c, dir);
						for (int s = 0; s<SG[prn].NSSC; s++)
							score += SF[s]*SG[prn].ssc.first[zone][s];
					}
					break;
				}
				default :
					cout<<"\nError: undefined scoring mode: "<<rmode<<endl;
					exit(1);
			}
			if (SF!=NULL) delete[] SF;
			if (verbosity>1) cout<<setprecision(6)<<score<<endl;
			if (score>=maxscore)
			{
				maxscore = score;
				mymotif[k] = mrk;
			}
			if (score>=myscore[k][mrk])
				myscore[k][mrk] = score;
		}
	}
// 	normalize_scores(myscore[k], StructCount);
	for (int mrk = 0; mrk<StructCount; mrk++)
		if (myscore[k][mrk]>thresh) count++;
	return count;
}




void normalize_scores(double *(&scores), int n)
{
	double sum = 0.;
	for (int p = 0; p<n; p++)
	{
		if (scores[p]<0)
		{
			double offset = -scores[p];
			for (int q = 0; q<n; q++) scores[q] += offset;
		}
	}
	for (int q = 0; q<n; q++) sum += scores[q];
	for (int q = 0; q<n; q++) scores[q] /= sum;
}




double *getScoreFactors(int pattern, int NSSC, int *pNS, int NS, scoreMachine **(&sm), int **neighbours, string &sequence, string &profile,
	int length, int z, int centre, int dir)
{
	double *fvector = NULL;
	double *rawscore = NULL;
	rawscore = new double[NS];
	int entryIndex = 0;
	for (int p = 0; p<pairsCount; p++)
	{
		for (int sID = 0; sID<pNS[p]; sID++)
		{
			int profLength = (int)sm[p][sID].winlen*NUM_OF_AMINO_ACIDS;
			for (int offset = -neighbours[pattern][p]; offset<=neighbours[pattern][p]; offset++)
			{
				assert(entryIndex<NS);
				if (verbosity>2) cout<<"p"<<p<<"."<<sID<<" offset="<<offset<<endl;
				int winstart = getStartingPosition(sm[p][sID], centre+offset, dir);
				rawscore[entryIndex] = computeRawScore(sm[p][sID], sequence, profile, length, profLength, z, dir, winstart);
				entryIndex++;
			}
		}
	}
	fvector = computeScoreFactors(rawscore, NS, NSSC);
	if (rawscore!=NULL) delete[] rawscore;
	return fvector;
}




double computeRawScore(scoreMachine &sm, string &sequence, string &profile, int length, int profLength, int z, int dir, int start)
{
	double score = 0.;
	double *myProf = NULL;
	int winLength = (int)sm.winlen;
	if (sm.collapse) myProf = getCollapsedProfile(sequence, length, winLength, profLength, dir, start);
	else myProf = getProfile(profile, profLength, dir, start);
	score = sm.getScore(z, myProf, profLength, dir);
	if (verbosity>3)
	{
		cout<<"profile: ";
		for (int q = 0; q<profLength; q++)
		{
			if (q%NUM_OF_AMINO_ACIDS==0) cout<<endl;
			cout<<setprecision(2)<<fixed<<myProf[q]<<" ";
		}
		cout<<setprecision(6)<<fixed<<score<<endl;
	}
	if (myProf!=NULL) delete[] myProf;
	return score;
}




double *computeScoreFactors(double *(&rawscore), int NS, int NSSC)
{
	int SPID = 0;
	double *fvector = new double[NSSC];
	for (int s = 0; s<NSSC; s++) fvector[s] = 0.;
	for (int si = 0; si<NS; si++)
	{
		assert(SPID<NSSC);
		fvector[SPID++] = rawscore[si];
	}
	for (int si = 0; si<NS; si++)
	{
		for (int sj = si; sj<NS; sj++)
		{
			assert(SPID<NSSC);
			fvector[SPID++] = rawscore[si]*rawscore[sj];
		}
	}
	if (SPID==NSSC-1) fvector[SPID] = 1.;
	else if (SPID<NSSC-1) cout<<"\n>entries missing ["<<SPID<<" / "<<NSSC<<"]"<<endl;
	else if (SPID>NSSC-1) cout<<"\n>entries overflow ["<<SPID<<" / "<<NSSC<<"]"<<endl;
	return fvector;
}




void computeSolutionError(uint dim, double **(&dbgmtx), double *(&dbgq), double *wv)
{
	uint imax = 0; // debug mode: index of the most deviating solution element
	double d = 0., deltamax = 0.; // temporary solution value and maximum deviation
	for (uint i = 0; i<dim; i++)
	{
		d = 0.;
		for (uint j = 0; j<dim; j++)
		{
			if (i>j) d += ((double)dbgmtx[i][j])*wv[j];
			else d += ((double)dbgmtx[j][i])*wv[j];
		}
		if (fabs(d-dbgq[i])>deltamax)
		{
			deltamax = fabs(d-dbgq[i]);
			imax = i;
		}
	}
	cout.setf(ios::scientific, ios::floatfield); // floatfield set to scientific
	cout<<"relative linear system modal precision: "<<deltamax<<" / ";
	cout.setf(ios::fixed, ios::floatfield);	// floatfield set to fixed
	cout<<fabs(dbgq[imax])<<endl;
	memrelease(dbgmtx, (long unsigned int)dim);
	memrelease(dbgq, (long unsigned int)dim);
}




void updateTable(uint k, double score, resultsTable &rt, int pattern, int realpattern, ofstream &bs)
{
	char c[2] = ".";
	if (separator=="eu") strcpy(c, ",");
	if (k<scorePlotSize)
	{
		writeScore(bs, score, (const char*)c);
		bs<<'\t'<<flush;
	}
	if (score>thresh) // positive score
	{
		rt.mycnt++;
		if (pattern==realpattern)
		{
			if (k<scorePlotSize)
			{
				writeScore(bs, score, (const char*)c);
				bs<<'\t'<<flush;
			}
			rt.good++;
		}
		else if (k<scorePlotSize)
		{
			bs<<'\t'<<flush;
			writeScore(bs, score, (const char*)c);
		}
	}
	else // negative score
	{
		if (pattern==realpattern)
		{
			if (k<scorePlotSize)
			{
				bs<<'\t'<<flush;
				writeScore(bs, score, (const char*)c);
			}
			rt.miss++;
		}
		else if (k<scorePlotSize)
		{
			writeScore(bs, score, (const char*)c);
			bs<<'\t'<<flush;
		}
	}
	if (k<scorePlotSize)
	{
		bs<<'\t'<<flush;
		if (pattern==realpattern)
		{
			writeScore(bs, score, (const char*)c);
			bs<<'\t'<<flush;
		}
		else
		{
			bs<<'\t'<<flush;
			writeScore(bs, score, (const char*)c);
		}
		bs<<'\t'<<realpattern<<endl;
	}
}




void writeScore(ostream &os, double v, const char *c)
{
	string dataString;
	ostringstream datastream;
	datastream << v;
	long startIndex = datastream.str().find_first_of('.');
	if (startIndex != (long)string::npos) dataString = datastream.str().replace(startIndex, 1, c);
	else dataString = datastream.str();
	os<<dataString<<flush;
}




int** memresInt(long unsigned int dimen)
{
	int **Nw;
	assert((size_t)(dimen*sizeof(int))>0);
	Nw = (int **) malloc(dimen*sizeof(int*));
	if (Nw==NULL)
	{
		cerr<<"\nmatrix memory allocation failure00. Error#"<<errno<<endl;
		perror(" ");
		exit(1);
	}
	for (long unsigned int k = 0; k<dimen; k++)
	{
		Nw[k] = (int *) malloc(dimen*sizeof(int));
		if (Nw[k]==NULL)
		{
			cerr<<"matrix memory allocation failure01. Error#"<<errno<<endl;
			perror(" ");
			exit(1);
		}
		for (long unsigned int h=0; h<dimen; h++) Nw[k][h] = 0;
	}
	return Nw;
}




double** memresFloat(long unsigned int dimen)
{
	double **Nw;
	assert((size_t)(dimen*sizeof(double))>0);
	Nw = (double **) malloc(dimen*sizeof(double*));
	if (Nw==NULL)
	{
		cerr<<"\nmatrix memory allocation failure00. Error#"<<errno<<endl;
		perror(" ");
		exit(1);
	}
	for (long unsigned int k = 0; k<dimen; k++)
	{
		Nw[k] = (double *) malloc((k+1)*sizeof(double));
		if (Nw[k]==NULL)
		{
			cerr<<"matrix memory allocation failure01. Error#"<<errno<<endl;
			perror(" ");
			exit(1);
		}
		for (long unsigned int h = 0; h<k+1; h++) Nw[k][h] = 0.;
	}
	return Nw;
}




double* memresVec(long unsigned int dimen)
{
	double *q;
	assert((size_t)(dimen*sizeof(double))>0);
	q = (double *) malloc(dimen*sizeof(double));
	if (q==NULL)
	{
		fprintf(stderr,"\nvector memory allocation failure. Error#%d\n",errno);
		perror(" ");
		exit(1);
	}
	for (long unsigned int k=0; k<dimen; k++) q[k] = 0.;
	return q;
}




void memrelease(double *(&v), long unsigned int dimen)
{
	if (v!=NULL) free(v);
}




void memrelease(int **(&v), long unsigned int dimen)
{
	for (long unsigned int i = 0; i<dimen; i++)
	if (v[i]!=NULL) free(v[i]);
	if (v!=NULL) free(v);
}




void memrelease(double **(&v), long unsigned int dimen)
{
	for (long unsigned int i = 0; i<dimen; i++)
	if (v[i]!=NULL) free(v[i]);
	if (v!=NULL) free(v);
}




void reInit(double **(&v), long unsigned int dim)
{
	for (long unsigned int k = 0; k<dim; k++)
	for (long unsigned int h = 0; h<k+1; h++) v[k][h] = 0.;
}




void reInit(double *(&v), long unsigned int dim)
{
	for (long unsigned int k = 0; k<dim; k++) v[k] = 0.;
}




void writeVector(const char *filename, int *v, uint size, char format)
{
	FILE *fp;
	fp = fopen(filename, "wb");
	if (fp==NULL) {cerr<<"\nwriting failure."<<endl; exit(1);}
	for (uint i = 0; i<size; i++)
	{
		if (format!='+') fprintf(fp, "%d\n", v[i]);
		else fprintf(fp, "%+d\n", v[i]);
	}
	fclose(fp);
}


uint readVector(const char *filename, int *v, uint size)
{
	FILE *fp;
	uint i = 0;
	fp = fopen(filename, "rb");
	if (fp==NULL) {cerr<<"\nreading failure."<<endl; return i;}
	for (i = 0; i<size; i++)
	if (fscanf(fp, "%d", &v[i])!=1) {cout<<"\nbad input."<<endl; return i;}
	fclose(fp);
	return size;
}


void writeVector(const char *filename, double *v, uint size, char format)
{
	FILE *fp;
	fp = fopen(filename, "wb");
	if (fp==NULL) {cerr<<"\nwriting failure."<<endl; exit(1);}
	for (uint i = 0; i<size; i++)
	{
		if (format!='+') fprintf(fp, "%lf\n", v[i]);
		else fprintf(fp, "%+lf\n", v[i]);
	}
	fclose(fp);
}


uint readVector(const char *filename, double *v, uint size)
{
	FILE *fp;
	uint i = 0;
	fp = fopen(filename, "rb");
	if (fp==NULL) {cerr<<"\nreading failure."<<endl; return i;}
	for (i = 0; i<size; i++)
	if (fscanf(fp, "%lf", &v[i])!=1) {cout<<"\nbad input."<<endl; return i;}
	fclose(fp);
	return size;
}




void copyArray(double **sourceArray, double **(&destinationArray), int dim)
{
	for (int i = 0; i<dim; i++)
		for (int j = i; j<dim; j++)
			destinationArray[j][i] = sourceArray[j][i];
}




void copyArray(double *sourceArray, double *(&destinationArray), int dim)
{
	for (int i = 0; i<dim; i++)
		destinationArray[i] = sourceArray[i];
}




void writeHitCounter(uint **array)
{
	ofstream fp("confidence-hit_table.dat");
	if (fp)
	{
		fp<<"confidenceLevel\t";
		for (int p = 0; p<StructCount; p++)
		{
			fp<<sstruc2char(p)<<" (hits)\t"<<sstruc2char(p)<<" (predicted)\t"<<sstruc2char(p)<<" (total)\t";
			fp<<sstruc2char(p)<<" (pred. ratio)\t"<<sstruc2char(p)<<" (total ratio)\t";
		}
		fp<<endl;
		for (uint confidenceLevel = 0; confidenceLevel<=confidenceMaxLevel; confidenceLevel++)
		{
			fp<<(double)confidenceLevel*confidenceStep<<'\t';
			for (int p = 0; p<StructCount; p++)
			{
				fp<<array[p][confidenceLevel]<<'\t';
				fp<<array[StructCount+p][confidenceLevel]<<'\t';
				fp<<array[2*StructCount+p][confidenceLevel]<<'\t';
				if (array[StructCount+p][confidenceLevel]!=0)
					fp<<(double)array[p][confidenceLevel]/(double)array[StructCount+p][confidenceLevel];
				fp<<'\t';
				if (array[2*StructCount+p][confidenceLevel]!=0)
					fp<<(double)array[p][confidenceLevel]/(double)array[2*StructCount+p][confidenceLevel];
				fp<<'\t';
			}
			fp<<endl;
		}
		fp.close();
	}
}




void writeAminoAcidsStats(resultsTable *(&rt))
{
	ofstream fp("aminoacids-stats.dat");
	uint **total = new uint*[NUM_OF_AMINO_ACIDS+1];
	uint **predicted = new uint*[NUM_OF_AMINO_ACIDS+1];
	for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
	{
		total[aa] = new uint[StructCount];
		predicted[aa] = new uint[StructCount];
	}
	if (fp)
	{
		fp<<"aatype\t";
		for (int p = 0; p<StructCount; p++)
		{
			fp<<sstruc2char(p)<<" (hits)\t"<<sstruc2char(p)<<" (total)\t"<<sstruc2char(p)<<" (predicted)\t";
			fp<<sstruc2char(p)<<" (total ratio)\t"<<sstruc2char(p)<<" (pred. ratio)\t";
		}
		fp<<endl;
		for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
		{
			fp<<aminoacid::ascode(aa)<<'\t';
			for (int p = 0; p<StructCount; p++)
			{
				total[aa][p] = 0;
				predicted[aa][p] = 0;
				for (int q = 0; q<StructCount; q++)
				{
					total[aa][p] += rt[p].confusionVector[aa][q];
					predicted[aa][p] += rt[q].confusionVector[aa][p];
				}
				fp<<rt[p].confusionVector[aa][p]<<'\t';
				fp<<total[aa][p]<<'\t'<<predicted[aa][p]<<'\t';
				if (total[aa][p]!=0)
					fp<<(double)rt[p].confusionVector[aa][p]/(double)total[aa][p];
				fp<<'\t';
				if (predicted[aa][p]!=0)
					fp<<(double)rt[p].confusionVector[aa][p]/(double)predicted[aa][p];
				fp<<'\t';
			}
			fp<<endl;
		}
		fp.close();
	}
	for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
	{
		delete[] total[aa];
		delete[] predicted[aa];
	}
	delete[] total;
	delete[] predicted;
}




string encodeProfile(string &profline)
{
	int k = 0;
	string encodedprofline = "";
	istringstream pss(profline);
	while (pss && k<NUM_OF_AMINO_ACIDS)
	{
		short rawEntry = -1;
		pss>>rawEntry;
		encodedprofline += TRencode(rawEntry);
		k++;
	}
	return encodedprofline;
}




Protein *readChain(inputMode imode, string sequenceFile, string profileFile)
{
	Protein *tmp = new Protein;
	if (imode==FASTA_MODE)
	{
		ifstream seqfs(sequenceFile.c_str());
		if (seqfs.is_open())
		{
			string line;
			while (getline(seqfs, line))
			{
				size_t chposNA = line.find_first_of(commentchars);
				size_t chposA = line.find_first_not_of(commentchars);
				if (chposA!=string::npos && (chposNA==string::npos || chposA<chposNA))
				{
					for (uint k = 0; k<line.size(); k++)
					{
						if (line[k]!=' ' and line[k]!='\t') // TODO are there any other white spaces?
						{
							if (isalpha(line[k]))
							{
								char res = line[k];
								if (res=='U' or res=='u') res = 'x';
								tmp->sequence += res;
							}
							else tmp->sequence += seqbreakchar;
							tmp->positions += CHAIN_CENTER; // dummy value here
						}
					}
				}
			}
			seqfs.close();
		}
		else
		{
			cout<<"a problem occurred while opening the sequence file "<<sequenceFile<<endl;
			exit(1);
		}
	}
	ifstream profilefs(profileFile.c_str());
	if (profilefs.is_open())
	{
		string line;
		string profline;
		unsigned int res = 0;
		while (profilefs)
		{
			if (imode==FASTA_MODE)
			{
// 				in FASTA_MODE a sequence is available as well:
// 				the program slides along the sequence until the first valid amino acid entry is found.
// 				for every non valid amino acid entry a default line is added to the profile
				while (res<tmp->sequence.size() && !isalpha(tmp->sequence[res]))
				{
					tmp->profile += DEFAULT_PROFILE;
					res++;
				}
			}
			if (getline(profilefs, line))
			{
				profline = matchProfileFormat(line);
				if (imode==FASTA_MODE && res<tmp->sequence.size())
				{
					if (profline.size()>2)
					{
						if (profline[0] == toupper(tmp->sequence[res])) // sequence and profile match
						{
							string buffer = profline.substr(2);
							tmp->profile += encodeProfile(buffer);
							res++;
						}
						else
						{
							ostringstream errorMessage;
							errorMessage << "sequence and profile do not match in residue " << res << "!" << endl;
							errorMessage << "expected: " << tmp->sequence[res] << "; received: " << profline[0] << endl;
							errorMessage << "line: " << line;
							std::cerr << errorMessage.str().c_str() << std::endl;
							exit ( 1 );
						}
					}
				}
				else if (profline.size()>2)
				{
					tmp->sequence += profline[0];
					tmp->positions += CHAIN_CENTER;
					string buffer = profline.substr(2);
					tmp->profile += encodeProfile(buffer);
				}
			}
		}
	}
	else
	{
		cout<<"a problem occurred while opening the profile file "<<profileFile<<endl;
		exit(1);
	}
	return tmp;
}




int factorial(int number)
{
	int temp;
	if (number <= 1) return 1;
	temp = number * factorial(number - 1);
	return temp;
}




int binomialCoefficient(int n, int k)
{
	if (n<k)
	{
	cerr<<"\nerror: illegal quantities in binomial coefficient."<<endl;
	exit(1);
	}
	else if (n==k) return 1;
	else return factorial(n) / ( factorial(k) * factorial(n-k) );
}




double binomialProbability(int n, int k, double p)
{
	double b;
	if (p>1. || p<0.)
	{
	cerr<<"\nerror: illegal quantities in binomial probability."<<endl;
	exit(1);
	}
	b = pow(p,(double)k) * pow((1.-p),(double)(n-k));
	b *= binomialCoefficient(n, k);
	return b;
}




double determinant(double **a, int n)
{
	int i, j, j0, j1;
	double det = 0.;
	double **m = NULL;

	if (n<1) { /* Error */ }
	else if (n==1) det = a[0][0];
	else if (n==2) det = a[0][0]*a[1][1]-a[1][0]*a[0][1];
	else
	{
	det = 0;
	for (j0 = 0; j0<n; j0++)
	{
	m = (double**) malloc((n-1)*sizeof(double*));
	for (i = 0; i<n-1; i++) m[i] = (double*) malloc((n-1)*sizeof(double));
	for (i = 1; i<n; i++)
	{
	j1 = 0;
	for (j = 0; j<n; j++)
	{
	if (j==j0) continue;
	m[i-1][j1] = a[i][j];
	j1++;
	}
	}
	det += pow(-1.0,1.0+j0+1.0) * a[0][j0] * determinant(m,n-1);
	for (i = 0; i<n-1; i++) free(m[i]);
	free(m);
	}
	}
	return(det);
}




double computeCorrelationCoefficient(resultsTable *(&rC), int catalogSize)
{
	double genMCC = 0.;
	double **confusionMatrix = NULL;
	confusionMatrix = new double*[catalogSize];
	for (int k = 0; k<catalogSize; k++)
	{
		confusionMatrix[k] = NULL;
		confusionMatrix[k] = new double[catalogSize];
		for (int h = 0; h<catalogSize; h++)
		{
			confusionMatrix[k][h] = 0.;
			for (int aa = 0; aa<NUM_OF_AMINO_ACIDS+1; aa++)
				confusionMatrix[k][h] += (double)rC[h].confusionVector[aa][k];
		}
	}
	genMCC = computeCorrelationCoefficient(confusionMatrix, catalogSize);
	if (verbosity>0) printConfusionMatrix(confusionMatrix, catalogSize);
	for (int k = 0; k<catalogSize; k++)
		if (confusionMatrix[k]!=NULL) delete[] confusionMatrix[k];
	if (confusionMatrix!=NULL) delete[] confusionMatrix;
	return genMCC;
}




double computeCorrelationCoefficient(double **(&confusionMatrix), int catalogSize)
{
	double numerator = 0.;
	double denominator = 0.;
	double sum_over_k_A = 0.;
	double sum_over_k_B = 0.;
	for (int k = 0; k<catalogSize; k++)
	{
		double sum_over_i_A = 0.;
		double sum_over_i_B = 0.;
		double sum_over_ij_A = 0.;
		double sum_over_ij_B = 0.;
		for (int i = 0; i<catalogSize; i++)
		{
			for (int j = 0; j<catalogSize; j++)
			{
				numerator += confusionMatrix[k][k]*confusionMatrix[i][j] - confusionMatrix[k][i]*confusionMatrix[j][k];
				if (j!=k)
				{
					sum_over_ij_A += confusionMatrix[j][i];
					sum_over_ij_B += confusionMatrix[i][j];
				}
			}
			sum_over_i_A += confusionMatrix[k][i];
			sum_over_i_B += confusionMatrix[i][k];
		}
		sum_over_k_A += sum_over_i_A*sum_over_ij_A;
		sum_over_k_B += sum_over_i_B*sum_over_ij_B;
	}
	denominator = sqrt(sum_over_k_A*sum_over_k_B+TINY);
	return numerator/denominator;
}




void printConfusionMatrix(double **(&confusionMatrix), int catalogSize)
{
	ofstream fs("confusionMatrix.csv", ios::app);
	fs<<"# columns mean >predicted as<, lines span the actual secondary structure motifs."<<endl;
	fs<<"# [for the Q3 case the order is helix, strand, coil]"<<endl;
	for (int k = 0; k<catalogSize; k++)
	{
		for (int h = 0; h<catalogSize; h++)
			fs<<confusionMatrix[h][k]<<'\t';
		fs<<endl;
	}
	fs.close();
}




void waitChild(int pid)
{
	int wpid, status;
	if (!isotropy && pid>0)
	{
	wpid = waitpid(pid, &status, 0);
	if (wpid == -1)
	perror("waitpid");
	if (WIFEXITED(status)) {
	printf("child exited, status=%d\n", WEXITSTATUS(status));

	} else if (WIFSIGNALED(status)) {
	printf("child killed (signal %d)\n", WTERMSIG(status));
	exit(1);
	} else if (WIFSTOPPED(status)) {
	printf("child stopped (signal %d)\n", WSTOPSIG(status));
	exit(1);
	}
	}
}




uint computeTotalResidues(Protein *(&pept))
{
	uint total = pept->sequence.length()-2*(uint)indent;
	return total;
}




int comparefloat (const void *a, const void *b)
{
        int response = 0;
        if ( *(double*)a - *(double*)b > 0. ) response = 1;
        else if ( *(double*)a - *(double*)b < 0. ) response = -1;
        return response;
}




uint fastaTransform(double input)
{
	uint output = 0;
	if (input>=1.) output = 9;
	else if (input>=0.) output = (uint)(input*10.);
	return output;
}

