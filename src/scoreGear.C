#include <fstream>
#include <iostream>
#include "scoreGear.h"
#include "toolz.h"

using namespace std;

scoreMachine::scoreMachine()
{
	this->Xdx = NULL;
	this->Xsx = NULL;
	this->Wdx = NULL;
	this->Wsx = NULL;
	this->bdx = NULL;
	this->bsx = NULL;
	this->reweight = NULL;
	this->power = true;
	this->linear = true;
	this->linsave = false;
	this->collapse = false;
	this->dim = 0;
	this->Ntot = 0;
	this->Mtot = 0;
	this->winlen = SPARROW_WINDOW_LENGTH;
	this->winkey = MAXLENGTH;
	this->isSet = false;
	this->isInit = false;
	this->isAvailable = false;
}



scoreMachine::~scoreMachine()
{
	if (this->isInit)
		delete[] this->reweight;
}



void scoreMachine::init(int first, int second)
{
// 	default values
	this->Pattern = first;
	this->Antipattern = second;
	this->reweight = new double[patternCount];
	for (int prn = 0; prn<patternCount; prn++) this->reweight[prn] = 0.5;
	if (first>=0 && first<patternCount) this->reweight[first] = 1.0;
	this->isSet = false;
	this->isInit = true;
	this->isAvailable = true;
}



void scoreMachine::setArraysSpace()
{
	assert(zoneCount>-1);
	uint *N = new uint[this->winlen];
	uint *M = new uint[this->winlen];
	for (uint qq = 0; qq<this->winlen; qq++)
	{
		N[qq] = (uint)(this->linear?NUM_OF_AMINO_ACIDS:0);
		M[qq] = (uint)(this->power?NUM_OF_AMINO_ACIDS:0);
		this->Ntot += N[qq];
		this->Mtot += M[qq];
	}
	this->dim = (uint)(this->Ntot+(this->Mtot*(this->Mtot+1))/2);
	if (this->winkey>=this->winlen) this->winkey = (uint)(this->winlen/2);
	delete[] N;
	delete[] M;
	Xdx = (double **) malloc((zoneCount+1)*sizeof(double*));
	Xsx = (double **) malloc((zoneCount+1)*sizeof(double*));
	Wdx = (double **) malloc((zoneCount+1)*sizeof(double*));
	Wsx = (double **) malloc((zoneCount+1)*sizeof(double*));
	bdx = (double *) malloc((zoneCount+1)*sizeof(double));
	bsx = (double *) malloc((zoneCount+1)*sizeof(double));
	if (Xdx==NULL or Xsx==NULL or Wdx==NULL or Wsx==NULL or bdx==NULL or bsx==NULL)
	{
		cout<<"\nError: scoring function memory allocation problem (SF)."<<endl;
		exit(1);
	}
	for (int z = 0; z<zoneCount+1; z++)
	{
		Xdx[z] = (double *) malloc(dim*sizeof(double));
		Xsx[z] = (double *) malloc(dim*sizeof(double));
		Wdx[z] = (double *) malloc(dim*sizeof(double));
		Wsx[z] = (double *) malloc(dim*sizeof(double));
		bdx[z] = 0.;
		bsx[z] = 0.;
		if (Xdx[z]==NULL or Xsx[z]==NULL or Wdx[z]==NULL or Wsx[z]==NULL)
		{
		cout<<"\nError: scoring function memory allocation problem (SFz"<<z<<")."<<endl;
		exit(1);
		}
		for (uint d = 0; d<dim; d++)
		{
		Xdx[z][d] = 0.;
		Xsx[z][d] = 0.;
		Wdx[z][d] = 0.;
		Wsx[z][d] = 0.;
		}
	}
	this->isSet = true;
}



void scoreMachine::releaseArraysSpace()
{
	if (this->isSet)
	{
		for (int z = 0; z<zoneCount+1; z++)
		{
			free(Xdx[z]);
			free(Xsx[z]);
			free(Wdx[z]);
			free(Wsx[z]);
		}
		free(Xdx);
	free(Xsx);
	free(Wdx);
	free(Wsx);
	free(bdx);
	free(bsx);
	}
}



void scoreMachine::getComponents(readMode rmode)
{
	char namedx[111];
	char namesx[111];
	char ext[5] = ".prm";
	for (int z = 0; z<zoneCount+1; z++)
	{
		if (buildFilename(namedx, ext, this, subdim, z, 1, rmode)!=0)
			buildFilename(namedx, ext, this, subdim, zoneCount, 1, rmode);
		if (buildFilename(namesx, ext, this, subdim, z, -1, rmode)!=0)
			buildFilename(namesx, ext, this, subdim, zoneCount, -1, rmode);
		if (verbosity>1) cout<<"reading file \""<<namedx<<"\"...."<<endl;
		if (isotropy) strcpy(namesx, namedx);
		if (getSolutions(namedx, this->bdx[z], this->Wdx[z], this->dim)<0
			|| getSolutions(namesx, this->bsx[z], this->Wsx[z], this->dim)<0)
		{
			cout<<"\nError: couldn't retrieve scoring function components."<<endl;
			cout<<"Aborted."<<endl;
			exit(1);
		}
	}
}



double scoreMachine::computeConstTerm(int pattern, int z, int dir)
{
	double *w = dir==1 ? this->Wdx[z] : this->Wsx[z];
	double *x = dir==1 ? this->Xdx[z] : this->Xsx[z];
	double sum = this->reweight[pattern];
	for (uint k = 0; k<this->dim; k++) sum -= w[k]*x[k];
	this->bdx[z] = sum;
	return sum;
}



double scoreMachine::getScore(int zone, double *prof, int profLength, int dir)
{
	double sum = 0.;
	int jm, jpmp, jluke; // sequence vector a.a./position linear and quadratic component}
	if (!this->isAvailable) return sum;
	for (jm = 0; jm<profLength; jm++)
	{
	if (linear)
	{
	if (dir==1) sum += Wdx[zone][jm]*prof[jm];
	else sum += Wsx[zone][jm]*prof[jm];
	}
	if (power)
	{
	jluke = (int)Ntot+jm;
	for (jpmp = 0; jpmp<profLength; jpmp++)
	{
	jluke += jpmp;
	if (jpmp>=jm)
	{
	if (dir==1) sum += Wdx[zone][jluke]*prof[jm]*prof[jpmp];
	else sum += Wsx[zone][jluke]*prof[jm]*prof[jpmp];
	}
	}
	}
	}
	if (dir==1) return sum+bdx[zone];
	else return sum+bsx[zone];
}



bool scoreMachine::isFeasible()
{
	unsigned int patternLength = (unsigned int)((subdim-1)*(gapsize+1)+1);
	if (this->winlen<patternLength)
	{
		cout<<"inadequate window length detected."<<endl;
		this->isAvailable = false;
		return false;
	}
	else if (this->winlen%2!=patternLength%2)
	{
		cout<<"window/pattern parity mismatch detected."<<endl;
		this->isAvailable = false;
		return false;
	}
	else return true;
}




scoreGear::scoreGear()
{
	this->sm = NULL;
	this->auxn_sm = NULL;
	this->auxa_sm = NULL;
	this->partNS = NULL;
	this->partauxnNS = NULL;
	this->partauxaNS = NULL;
	this->totalNS = 0;
	this->totalauxnNS = 0;
	this->totalauxaNS = 0;
	this->NSSC = 0;
	this->auxnNSSC = 0;
	this->auxaNSSC = 0;
	this->ssc.first = NULL;
	this->ssc.second = NULL;
	this->auxn_ssc.first = NULL;
	this->auxn_ssc.second = NULL;
	this->auxa_ssc.first = NULL;
	this->auxa_ssc.second = NULL;
	this->Pattern = patternCount;
	this->aux = false;
}



scoreGear::~scoreGear() {}



void scoreGear::setup(bool use_aux)
{
	this->aux = use_aux;
	this->setScoreTools();
	this->getSuperCoefficients(STD_MODE, this->NSSC);
	if (this->aux)
	{
		this->getSuperCoefficients(AUXN_MODE, this->auxnNSSC);
		this->getSuperCoefficients(AUXA_MODE, this->auxaNSSC);
	}
}



void scoreGear::setScoreTools()
{
	this->partNS = new int[pairsCount];
	this->partauxnNS = new int[pairsCount];
	this->partauxaNS = new int[pairsCount];
	for (int pID = 0; pID<pairsCount; pID++)
	{
		this->partNS[pID] = SPARROW_SCORE_TYPES;
		this->partauxnNS[pID] = 0;
		this->partauxaNS[pID] = 0;
// 		cout<<"pattern "<<this->Pattern<<" pair "<<pID<<" scorers count: "<<this->partNS[pID]<<endl;
	}
	this->totalNS = computeEffectiveNS(this->Pattern, this->partNS, neighboursCount);
	this->totalauxnNS = computeEffectiveNS(this->Pattern, this->partauxnNS, auxnneighboursCount);
	this->totalauxaNS = computeEffectiveNS(this->Pattern, this->partauxaNS, auxaneighboursCount);
	this->NSSC = computeNQuadratic(this->totalNS);
// 	cout<<"total number of super-parameters: "<<this->NSSC<<endl;
	if (this->aux)
	{
		this->auxnNSSC = computeNQuadratic(this->totalauxnNS);
		this->auxaNSSC = computeNQuadratic(this->totalauxaNS);
// 		cout<<"total number of n-auxiliary super-parameters: "<<this->auxnNSSC<<endl;
// 		cout<<"total number of a-auxiliary super-parameters: "<<this->auxaNSSC<<endl;
	}
	this->sm = new scoreMachine*[pairsCount];
	this->auxn_sm = new scoreMachine*[pairsCount];
	this->auxa_sm = new scoreMachine*[pairsCount];
	for (int p0 = 0; p0<patternCount; p0++)
	{
		int pID = 0;
		for (int p1 = p0+1; p1<patternCount; p1++)
		{
			pID = getMotifsPair(p0, p1);
			this->sm[pID] = this->getScorers(STD_MODE, this->partNS[pID], p0, p1);
			this->auxn_sm[pID] = this->getScorers(AUXN_MODE, this->partauxnNS[pID], p0, p1);
			this->auxa_sm[pID] = this->getScorers(AUXA_MODE, this->partauxaNS[pID], p0, p1);
		}
		pID = getMotifsPair(p0);
		this->sm[pID] = this->getScorers(STD_MODE, this->partNS[pID], p0);
		this->auxn_sm[pID] = this->getScorers(AUXN_MODE, this->partauxnNS[pID], p0);
		this->auxa_sm[pID] = this->getScorers(AUXA_MODE, this->partauxaNS[pID], p0);
	}
}



scoreMachine *scoreGear::getScorers(readMode rmode, int NS, int p0, int p1)
{
	int pID = 0;
	scoreMachine *tmpsm = NULL;
	pID = getMotifsPair(p0, p1);
	if (NS>0)
	{
		tmpsm = new scoreMachine[NS];
		for (int s = 0; s<NS; s++)
		{
			tmpsm[s].init(p0, p1);
			tmpsm[s].setArraysSpace();
			if (tmpsm[s].isFeasible())
			{
				tmpsm[s].getComponents(rmode);
				if (verbosity>4)
				{
					stringstream ss;
					ss<<"sc"<<this->Pattern<<"."<<pID<<"."<<s<<".dat";
					writeVector(ss.str().c_str(), tmpsm[s].Wdx[zoneCount], tmpsm[s].dim);
				}
			}
		}
	}
	return tmpsm;
}



void scoreGear::getSuperCoefficients(readMode rmode, int N)
{
	string tag = "";
	if (N>1)
	{
		pair <double**, double**> rssc;
		switch(rmode)
		{
			case STD_MODE :
				break;
			case AUXN_MODE :
				tag = ".auxn";
				break;
			case AUXA_MODE :
				tag = ".auxa";
				break;
			default :
				cout<<"\nError: undefined read mode: "<<rmode<<endl;
				exit(1);
		}
		rssc.first = new double*[zoneCount+1];
		rssc.second = new double*[zoneCount+1];
		for (int z = 0; z<=zoneCount; z++)
		{
			rssc.first[z] = NULL;
			rssc.second[z] = NULL;
			rssc.first[z] = new double[N];
			rssc.second[z] = new double[N];
			for (int i = 0; i<N; i++)
			{
				rssc.first[z][i] = 0.;
				rssc.second[z][i] = 0.;
			}
			this->getSuperCoefficients(z, tag, N, rssc, rmode);
		}
		switch(rmode)
		{
			case STD_MODE :
				this->ssc = rssc;
				break;
			case AUXN_MODE :
				this->auxn_ssc = rssc;
				break;
			case AUXA_MODE :
				this->auxa_ssc = rssc;
				break;
			default :
				cout<<"\nError: undefined read mode: "<<rmode<<endl;
				exit(1);
		}
	}
}



void scoreGear::getSuperCoefficients(int z, string &tag, int N, pair <double**, double**> &rssc, readMode rmode)
{
	char namedx[111];
	char namesx[111];
	sprintf(namedx, ".prm/p%04dlev%02dsize%d-%cdx%s.rmz", this->Pattern, detail, subdim, zname(z), tag.c_str());
	sprintf(namesx, ".prm/p%04dlev%02dsize%d-%csx%s.rmz", this->Pattern, detail, subdim, zname(z), tag.c_str());
	if (getSolutions(namedx, rssc.first[z][N-1], rssc.first[z], N-1)<0)
	{
		sprintf(namedx, ".prm/p%04dlev%02dsize%d-%cdx%s.rmz", this->Pattern, detail, subdim, zname(zoneCount), tag.c_str());
		if (getSolutions(namedx, rssc.first[z][N-1], rssc.first[z], N-1)<0)
		{
			cout<<"\nFile "<<namedx<<" error: impossible to procede."<<endl;
			exit(1);
		}
	}
	if (verbosity>4) cout<<"pattern="<<this->Pattern<<" zone="<<z<<" NC-super-parameters:"<<endl;
	for (int i = 0; i<N; i++)
	{
		rssc.second[z][i] = rssc.first[z][i];
		if (verbosity>4) cout<<rssc.first[z][i]<<endl;
	}
	if (!isotropy)
	{
		if (getSolutions(namesx, rssc.second[z][N-1], rssc.second[z], N-1)<0)
		{
			sprintf(namesx, ".prm/p%04dlev%02dsize%d-%csx%s.rmz", this->Pattern, detail, subdim, zname(zoneCount), tag.c_str());
			if (getSolutions(namesx, rssc.second[z][N-1], rssc.second[z], N-1)<0)
			{
				cout<<"\nFile "<<namesx<<" error: impossible to procede."<<endl;
				exit(1);
			}
		}
		if (verbosity>4)
		{
			cout<<"pattern="<<this->Pattern<<" zone="<<z<<" CN-super-parameters:"<<endl;
			for (int i = 0; i<N; i++)
				cout<<rssc.second[z][i]<<endl;
		}
	}
}



void scoreGear::eraseScoreTools()
{
	for (int pairID = 0; pairID<pairsCount; pairID++)
	{
		for (int s = 0; s<this->partNS[pairID]; s++) this->sm[pairID][s].releaseArraysSpace();
		for (int s = 0; s<this->partauxnNS[pairID]; s++) this->auxn_sm[pairID][s].releaseArraysSpace();
		for (int s = 0; s<this->partauxaNS[pairID]; s++) this->auxa_sm[pairID][s].releaseArraysSpace();
		if (this->sm[pairID]!=NULL) delete[] this->sm[pairID];
		if (this->auxn_sm[pairID]!=NULL) delete[] this->auxn_sm[pairID];
		if (this->auxa_sm[pairID]!=NULL) delete[] this->auxa_sm[pairID];
	}
	if (this->sm!=NULL) delete[] this->sm;
	if (this->auxn_sm!=NULL) delete[] this->auxn_sm;
	if (this->auxa_sm!=NULL) delete[] this->auxa_sm;
	if (this->partNS!=NULL) delete[] this->partNS;
	if (this->partauxnNS!=NULL) delete[] this->partauxnNS;
	if (this->partauxaNS!=NULL) delete[] this->partauxaNS;
}



void scoreGear::eraseSuperCoefficients()
{
	for (int z = 0; z<=zoneCount; z++)
	{
		if (this->ssc.first[z]!=NULL) delete[] this->ssc.first[z];
		if (this->ssc.second[z]!=NULL) delete[] this->ssc.second[z];
		if (this->aux)
		{
			if (this->auxn_ssc.first[z]!=NULL) delete[] this->auxn_ssc.first[z];
			if (this->auxn_ssc.second[z]!=NULL) delete[] this->auxn_ssc.second[z];
			if (this->auxa_ssc.first[z]!=NULL) delete[] this->auxa_ssc.first[z];
			if (this->auxa_ssc.second[z]!=NULL) delete[] this->auxa_ssc.second[z];
		}
	}
	if (this->ssc.first!=NULL) delete[] this->ssc.first;
	if (this->ssc.second!=NULL) delete[] this->ssc.second;
	if (this->aux)
	{
		if (this->auxn_ssc.first!=NULL) delete[] this->auxn_ssc.first;
		if (this->auxn_ssc.second!=NULL) delete[] this->auxn_ssc.second;
		if (this->auxa_ssc.first!=NULL) delete[] this->auxa_ssc.first;
		if (this->auxa_ssc.second!=NULL) delete[] this->auxa_ssc.second;
	}
}

