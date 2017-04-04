
#include <fstream>
#include <iostream>
#include "myClassifier.h"

using namespace std;

myClassifier::myClassifier(int NC, int ns, int SC, string myprompt)
{
	this->isInit = false;
	this->isUsable = false;
	this->classesCount = (uint)NC;
	this->featsCount = (uint)SC;
	this->cclass = new uint[ns];
	this->lambda = TINY;
	this->sampleID = new uint[this->classesCount];
	this->partCount = new uint[this->classesCount];
	this->partWeight = new double*[this->classesCount];
	this->featsVector = new double**[this->classesCount];
	this->featsWeightVector = new double*[this->classesCount];
	this->statsVector = new double*[this->classesCount];
	this->statsMatrix = new double**[this->classesCount];
	this->results = new resultsTable[this->classesCount];
	for (uint c = 0; c<this->classesCount; c++)
	{
		this->sampleID[c] = 0;
		this->partCount[c] = 0;
		this->partWeight[c] = new double[this->classesCount];
		for (uint d = 0; d<this->classesCount; d++) partWeight[c][d] = 1./(double)classesCount;
		this->featsVector[c] = NULL;
		this->featsWeightVector[c] = NULL;
		this->featsWeightVector[c] = new double[this->featsCount+1];
		for (uint nf = 0; nf<this->featsCount+1; nf++) this->featsWeightVector[c][nf] = 0.;
		this->statsVector[c] = memresVec((long unsigned int)this->featsCount);
		this->statsMatrix[c] = memresFloat((long unsigned int)this->featsCount);
		this->results[c].init((int)classesCount);
	}
	this->prompt = myprompt;
	this->solutionfilename = myprompt;
	this->cfgfilename = myprompt+".cfg";
}

myClassifier::~myClassifier()
{
	delete[] this->cclass;
	delete[] this->sampleID;
	for (uint c = 0; c<this->classesCount; c++)
	{
		delete[] this->partWeight[c];
		for (uint nf = 0; nf<this->featsCount; nf++)
			if (this->isInit && featsVector[c][nf]!=NULL) delete[] featsVector[c][nf];
		if (this->featsVector[c]!=NULL) delete[] this->featsVector[c];
		if (this->featsWeightVector[c]!=NULL) delete[] this->featsWeightVector[c];
		memrelease(this->statsVector[c], (long unsigned int)this->featsCount);
		memrelease(this->statsMatrix[c], (long unsigned int)this->featsCount);
	}
	delete[] this->featsVector;
	delete[] this->featsWeightVector;
	delete[] this->statsMatrix;
	delete[] this->statsVector;
	delete[] this->results;
	delete[] this->partWeight;
	delete[] this->partCount;
}


void myClassifier::init(int nsamples)
{
	uint total = 0;
	for (uint c = 0; c<this->classesCount; c++)
	{
		total += this->partCount[c];
		this->featsVector[c] = NULL;
		this->featsVector[c] = new double*[this->featsCount];
		for (uint nf = 0; nf<this->featsCount; nf++)
		{
			this->featsVector[c][nf] = NULL;
			if (this->partCount[c]>0) this->featsVector[c][nf] = new double[this->partCount[c]];
			for (uint i = 0; i<this->partCount[c]; i++) this->featsVector[c][nf][i] = 0.;
		}
	}
	if (total!=(uint)nsamples)
	{
		cout<<"\n>samples count mismatch: aborting...."<<endl;
		exit(1);
	}
	this->readParameters();
	this->isInit = true;
}


void myClassifier::readParameters()
{
	uint ct = 0;
	char ch = '*';
	string opt = "";
	stringstream ss;
	ifstream fp(this->cfgfilename.c_str(), ios::in);
	if (!fp) cout<<this->prompt<<": no configuration file found. setting default parameters...."<<endl;
	else // the input file exists, so we read the parameters from it
	{
	while (fp)
	{
	fp>>ch;
	if (ch=='#') getline(fp, opt);
	else
	{
	fp.putback(ch);
	fp>>opt;
	if (opt=="class")
	{
		fp>>ct;
		if (ct>=this->classesCount)
		{
			cout<<"\n>invalid class."<<endl;
			exit(1);
			}
		}
	else if (opt=="RWFs")
	{
	uint c = 0;
	double sum = 0.;
	getline(fp, opt);
	ss<<opt;
	while (c<this->classesCount && ss>>this->partWeight[ct][c])
	{
	if (this->partWeight[ct][c]<0.)
		this->partWeight[ct][c] = 0.;
	sum += this->partWeight[ct][c];
	c++;
	}
	sum = 1./(sum+TINY);
	for (c = 0; c<this->classesCount; c++)
	this->partWeight[ct][c] *= sum;
	}
	else if (opt=="lambda") fp>>this->lambda;
	else cout<<"unknown parameter: "<<opt<<endl;
	}
	}
	fp.close();
	}
}


void myClassifier::buildScoringFunctions()
{
	double *normFactor = new double[this->classesCount];
	double *dbgq = NULL; // debug right hand side vector
	double **dbgmtx = NULL; // debug matrix
	for (uint c = 0; c<this->classesCount; c++)
		myClassifier::learnStats(this->featsVector[c], this->featsCount, this->partCount[c], this->statsMatrix[c], this->statsVector[c]);
	for (uint ct = 0; ct<this->classesCount; ct++)
	{
		double *x = memresVec((long unsigned int)this->featsCount);
		double **a = memresFloat((long unsigned int)this->featsCount);
		if (verbosity>0)
		{
			dbgq = memresVec((long unsigned int)this->featsCount);
			dbgmtx = memresFloat((long unsigned int)this->featsCount);
		}
		for (uint c = 0; c<this->classesCount; c++) normFactor[c] = partWeight[ct][c]/((double)partCount[c]+TINY);
		this->featsWeightVector[ct][this->featsCount] = this->partWeight[ct][ct];
		for (uint nf = 0; nf<this->featsCount; nf++)
		{
			for (uint c = 0; c<this->classesCount; c++) x[nf] += normFactor[c]*this->statsVector[c][nf];
			this->featsWeightVector[ct][nf] = this->partWeight[ct][ct]*(this->statsVector[ct][nf]/((double)partCount[ct]+TINY)-x[nf]);
			if (verbosity>0) dbgq[nf] = this->featsWeightVector[ct][nf];
		}
		for (uint nf = 0; nf<this->featsCount; nf++)
		{
			for (uint mf = nf; mf<this->featsCount; mf++)
			{
				for (uint c = 0; c<this->classesCount; c++)
					a[mf][nf] += normFactor[c]*this->statsMatrix[c][mf][nf];
				a[mf][nf] -= x[mf]*x[nf];
				if (mf==nf) a[mf][nf] += this->lambda;
				if (verbosity>0) dbgmtx[mf][nf] = a[mf][nf];
			}
		}
		if (verbosity>0) cout<<prompt<<"> solving linear equations system...."<<endl;
		LLt_factorize(a, (long int)this->featsCount);
		LLt_backsubst(a, (long int)this->featsCount, this->featsWeightVector[ct]);
		for (uint nf = 0; nf<this->featsCount; nf++) this->featsWeightVector[ct][this->featsCount] -= featsWeightVector[ct][nf]*x[nf];
		if (verbosity>0) computeSolutionError(this->featsCount, dbgmtx, dbgq, this->featsWeightVector[ct]);
		memrelease(x, (long unsigned int)this->featsCount);
		memrelease(a, (long unsigned int)this->featsCount);
	}
	delete[] normFactor;
}


void myClassifier::addFeatureInstance(int cc, uint featID, double feat)
{
	if (this->sampleID[cc]<this->partCount[cc])
		this->featsVector[cc][featID][this->sampleID[cc]] = feat;
	else
	{
		cout<<"\n>sample count override."<<endl;
		exit(1);
	}
}


void myClassifier::readSolutions()
{
	char buf[9];
	for (uint ct = 0; ct<this->classesCount; ct++)
	{
		sprintf(buf, "%u.prm", ct);
		string tmpname = this->solutionfilename+buf;
		getSolutions(tmpname.c_str(), this->featsWeightVector[ct][this->featsCount], this->featsWeightVector[ct], this->featsCount);
	}
	this->isUsable = true;
}


void myClassifier::writeSolutions()
{
	char buf[9];
	for (uint ct = 0; ct<this->classesCount; ct++)
	{
		sprintf(buf, "%u.prm", ct);
		string tmpname = this->solutionfilename+buf;
		storeSolutions(tmpname.c_str(), this->featsWeightVector[ct][this->featsCount], this->featsWeightVector[ct], this->featsCount);
	}
}


void myClassifier::incrementSampleID(int cc)
{
	if (this->sampleID[cc]<this->partCount[cc]) sampleID[cc]++;
}


void myClassifier::updateTable(int myClass, int realClass)
{
	this->results[myClass].mycnt++;
	this->results[realClass].cnt++;
	if (myClass==realClass)
		this->results[realClass].good++;
	else this->results[realClass].miss++;
}


int myClassifier::getClass(double *targetVector)
{
	int winnerClass = classesCount;
	double maxscore = -11.;
	for (uint ct = 0; ct<this->classesCount; ct++)
	{
		double score = this->featsWeightVector[ct][this->featsCount];
		for (uint nf = 0; nf<this->featsCount; nf++) score += targetVector[nf]*this->featsWeightVector[ct][nf];
		if (this->isUsable && verbosity>0) cout<<"class "<<ct<<" score: "<<score<<endl;
		if (score>=maxscore)
		{
			winnerClass = ct;
			maxscore = score;
		}
	}
	return winnerClass;
}


double myClassifier::getClassScore(int cclass, double *targetVector)
{
	double score = this->featsWeightVector[cclass][this->featsCount];
	for (uint nf = 0; nf<this->featsCount; nf++)
		score += targetVector[nf]*this->featsWeightVector[cclass][nf];
	return score;
}


int myClassifier::learnStats(double **featsArray, uint NF, uint size, double **(&Nw), double *(&Qw), uint Nlin, bool lin, bool pow)
{
	uint jm, jpmp, jsms, jtmt; // sequence vector a.a./position linear and quadratic component
	uint jluke = 0, jleia = 0; // special indexes
	double *tmpQuadEntry = new double[size]; // stores temporary products of two featsArray values
	double *tmpCubeEntry = new double[size]; // stores temporary products of three featsArray values
	for (jm = 0; jm<NF; jm++)	// first residue loop
	{
	if (lin)
	{
	for (int i = 0; i<(int)size; i++)
	Qw[jm] += featsArray[jm][i];
	}
	jluke = Nlin+jm;
	for (jpmp = 0; jpmp<NF; jpmp++)	// second residue loop
	{
	jluke += jpmp;
	if (jpmp>=jm)
	{
	for (int i = 0; i<(int)size; i++)
	{
	tmpQuadEntry[i] = featsArray[jm][i]*featsArray[jpmp][i];
	if (lin) Nw[jpmp][jm] += tmpQuadEntry[i];
	if (pow) Qw[jluke] += tmpQuadEntry[i];
	}
	}
	if (pow)
	{
	for (jsms = 0; jsms<NF; jsms++)	// third residue loop
	{
	jleia = Nlin+jsms;
	if (jpmp>=jm)
	{
	for (int i = 0; i<(int)size; i++)
	{
	tmpCubeEntry[i] = tmpQuadEntry[i]*featsArray[jsms][i];
	if (lin) Nw[jluke][jsms] += tmpCubeEntry[i];
	}
	}
	for (jtmt = 0; jtmt<NF; jtmt++)	// fourth residue loop
	{
	jleia += jtmt;
	if (jtmt>=jsms && jpmp>=jm && jluke>=jleia)
	for (int i = 0; i<(int)size; i++)
	Nw[jluke][jleia] += featsArray[jtmt][i]*tmpCubeEntry[i];
	}
	}
	}
	}
	}
	delete[] tmpQuadEntry;
	delete[] tmpCubeEntry;
	return 0;
}


