#include <iostream>
#include <fstream>
#include "aatype.h"

using namespace std;

const string aminoacid::onelettercodes = "ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv";

double* aminoacid::logValue = NULL;

aminoacid::aminoacid(){}
aminoacid::~aminoacid(){}

void aminoacid::init()
{
	logValue = new double[maxValue];
	for (int i = 0; i<(int)maxValue; i++)
	logValue[i] = 1./(1.+exp((double)((int)maxValue/2-i)));
}

void aminoacid::kill()
{
	delete[] logValue;
}


int aminoacid::astype(char z)
{
	int j;
	switch (z)
	{
	case 'A': {j=0;	break;}
	case 'R': {j=1;	break;}
	case 'N': {j=2;	break;}
	case 'D': {j=3;	break;}
	case 'C': {j=4;	break;}
	case 'Q': {j=5;	break;}
	case 'E': {j=6;	break;}
	case 'G': {j=7;	break;}
	case 'H': {j=8;	break;}
	case 'I': {j=9;	break;}
	case 'L': {j=10; break;}
	case 'K': {j=11; break;}
	case 'M': {j=12; break;}
	case 'F': {j=13; break;}
	case 'P': {j=14; break;}
	case 'S': {j=15; break;}
	case 'T': {j=16; break;}
	case 'W': {j=17; break;}
	case 'Y': {j=18; break;}
	case 'V': {j=19; break;}
	default : {j=20;}
	}
	return j;
}


char aminoacid::ascode(int a)
{
	char c;
	switch (a)
	{
	case 0:	{c='A'; break;}
	case 1:	{c='R'; break;}
	case 2:	{c='N'; break;}
	case 3:	{c='D'; break;}
	case 4:	{c='C'; break;}
	case 5:	{c='Q'; break;}
	case 6:	{c='E'; break;}
	case 7:	{c='G'; break;}
	case 8:	{c='H'; break;}
	case 9:	{c='I'; break;}
	case 10: {c='L'; break;}
	case 11: {c='K'; break;}
	case 12: {c='M'; break;}
	case 13: {c='F'; break;}
	case 14: {c='P'; break;}
	case 15: {c='S'; break;}
	case 16: {c='T'; break;}
	case 17: {c='W'; break;}
	case 18: {c='Y'; break;}
	case 19: {c='V'; break;}
	default : {c='X';}
	}
	return c;
}



int *aminoacid::astype(char *z, int length)
{
	int *ncode = NULL;
	ncode = new int[length];
	for (int k = 0; k<length; k++)
		ncode[k] = astype(z[k]);
	return ncode;
}


char *aminoacid::ascode(int *a, int length)
{
	char *lcode = NULL;
	lcode = new char[length];
	for (int k = 0; k<length; k++)
		lcode[k] = ascode(a[k]);
	return lcode;
}

