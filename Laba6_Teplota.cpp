#include "pch.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#define PRES double
#define NXB 15
#define NX NXB*3+1
#define NYB 12
#define NY NYB*3+1
#define REP 3000
#define EPSL 1.e-5
#define LL 1.7
#define TEM1 5.0
#define TEM2 15.0
#define HX 0.2
#define HY 0.3

using namespace std;

void maxpvr(PRES *t1, PRES *del, PRES *maxdel)
{
	PRES d = fabs(*del) / fabs(*t1);
	if (d > *maxdel) 
		*maxdel = d;
}

int main(int argc, char **argv)
{
	ofstream foutT("c:/work/dT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	int i, i1, i2, i3,
		j, j1, j2, j3,
		k = 0, rp = 0;
	bool v = 0;
	PRES T1 = TEM1, T2 = TEM2, h = HX, r = HY, tx, t0, t1, del, maxdel = 0.0;
	PRES T[NY][NX];
	PRES lam = LL;
	PRES eps = EPSL;
	int prz = 1;
	int nT = 0;
	PRES alf_1 = -h / r;
	PRES alf_2 = -r / h;
	PRES alf_3 = alf_2 * 0.5;
	PRES alf_4 = alf_1 * 0.5;
	PRES gam_1 = -2 * (alf_1 + alf_2);
	PRES gam_2 = -1.5 * (alf_1 + alf_2);
	PRES gam_3 = -(alf_1 + alf_2);
	PRES gam_4 = -(alf_3 + alf_4);

	i1 = NXB;
	i2 = 2 * NXB;
	i3 = 3 * NXB;
	j1 = NYB; 
	j2 = 2 * NYB; 
	j3 = 3 * NYB;
	rp = REP;

	for (j = 0; j <= j3; j++) 
	{
		for (i = 0; i <= i3; i++) 
		{ 
			T[j][i] = 0.0; 
		}
	}
	//[A,B]
	for (i = 0; i <= i1; i++) 
		T[j3][i] = T1;
	//[I,J]
	for (i = i2; i <= i3; i++) 
		T[0][i] = T2;

	while (k < rp && prz == 1) 
	{
		k++;
		for (j = 0; j <= j3; j++) 
		{
			for (i = 0; i <= i3; i++) 
			{
				v = 0;
				t0 = T[j][i];
				tx = 0;
				//D
				if (i == i3 && j == j3)
					tx = -(alf_4*T[j - 1][i] + alf_3 * T[j][i - 1]) / gam_4;
				//E
				else if (i == 0 && j == j2)
					tx = -(alf_3*T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_4;
				//F
				else if (i == i1 && j == j2)
					tx = -(alf_4 * T[j-1][i] + alf_3 * T[j][i-1] + alf_2 * T[j][i+1] + alf_1 * T[j+1][i]) / gam_2;
				//G
				else if (i == i1 && j == j1)
					tx = -(alf_3*T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_4;
				//H
				else if (i == i2 && j == j1)
					tx = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_2;

				//горизонтали
				//BD
				else if (i1 < i && i < i3 && j == j3)
					tx = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j-1][i]) / gam_3;

				//EF
				else if (0 < i && i < i1 && j == j2)
					tx = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j+1][i]) / gam_3;

				//GH
				else if (i1 < i && i < i2 && j == j1)
					tx = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j+1][i]) / gam_3;

				//вертикали
				//AE
				else if (i == 0 && j2 < j && j < j3)
					tx = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_3;

				//FG
				else if (i == i1 && j1 < j && j < j2)
					tx = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_3;

				//HI
				else if (i == i2 && 0 < j && j < j1)
					tx = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) / gam_3;

				//DJ
				else if (i == i3 && 0 < j && j < j3)
					tx = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_4 * T[j + 1][i]) / gam_3;

				//области
				//ABEF
				else if (0 < i && i <= i1 && j2 < j && j < j3)
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;

				//BCHG
				else if (i1 < i && i <= i2 && j1 < j && j < j3)
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;

				//CDJI
				else if (i2 < i && i < i3 && 0 < j && j < j3)
					tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) / gam_1;

				if (tx != 0)
				{
					del = lam * (tx - t0);
					t1 = t0 + del;
					T[j][i] = t1;
					maxpvr(&t1, &del, &maxdel);
				}
			}
		}
		nT++; 
		PRES w = maxdel;
		foutT.write((char*)&w, sizeof w);
		if (maxdel < eps) 
			prz = 0; 
		maxdel = 0.0;
	}
	foutT.close();
	ofstream fouT("nT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fouT.write((char*)&nT, sizeof nT);
	fouT.close();
	ofstream fout("Pole.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	for (j = 0; j < NY; j++) 
	{
		for (i = 0; i < NX; i++)
		{ 
			PRES w = T[j][i];
			fout.write((char*)&w, sizeof w);
		}
	}
	fout.close();
	int n_x = NX; int n_y = NY;
	ofstream fou("Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
	fou.write((char*)&n_x, sizeof n_x);
	fou.write((char*)&n_y, sizeof n_y);
	fou.close();
}