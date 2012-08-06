/*
����2�����M�`��
custom-imp2.c
2004.6/26��	�Ƃ肠����:custom-imp.c
2004.10/4	�ꕔ����(�����t�@�C���Ȃǂ��܂�:custom-imp2.c)

�����@(�A��@)���g�p
�M�������l��
�����o�͌^

Nd=40,Ns=100��,1�P�ʎ��Ԃ̌v�Z��10�����x������
*/

#include <stdio.h>
#include "joken2.h"

// �o�̓t�@�C��(joken.h�̒��Œ�`����ꍇ�̓R�����g�A�E�g���Ă�������)
#define FILENAME "dat/imp-all2.dat"

double A[Nd*Ns][Nd*Ns];
double D[Nd*Ns];

#include "gauss.h"


void keisu(void) {
// �s��v�Z�pA,D�쐬
	int i, j;
	double Fd, Fs;
	
	for (i=0; i<Nd*Ns; i++) {
		for (j=0; j<Nd*Ns; j++) {
			A[i][j] = 0;
		}
	}
	
	for (i=0; i<Nd; i++) {
		for (j=0; j<Ns; j++) {
			coefficient(i,j);
			Fd = a*ht/(hd*hd);
			Fs = a*ht/(hs*hs);
			
			D[Nd*i+j] = Told[i][j] + ht*b*S[i][j]*(hd*hs);
			
			A[Nd*i+j][Nd*i+j] = 1 + 2*(Fd + Fs);
			switch (i) {
				case 0:
					A[Nd*i+j][Nd*(i+1)+j] = -2*Fd;
					break;
				case Nd-1:
					A[Nd*i+j][Nd*(i-1)+j] = -2*Fd;
					break;
				default:
					A[Nd*i+j][Nd*(i-1)+j] = -Fd;
					A[Nd*i+j][Nd*(i+1)+j] = -Fd;
			}
			switch (j) {
				case 0:
					A[Nd*i+j][Nd*i+(j+1)] = -2*Fs;
					break;
				case Ns-1:
					A[Nd*i+j][Nd*i+(j-1)] = -2*Fs;
					break;
				default:
					A[Nd*i+j][Nd*i+(j-1)] = -Fs;
					A[Nd*i+j][Nd*i+(j+1)] = -Fs;
			}
		}
	}
}

void replace(void) {
	int i, j;
	for (i=0; i<Nd; i++) {
		for (j=0; j<Ns; j++) {
			Tnew[i][j] = D[Nd*i+j];	// Tnew
		}
	}
}

void output(int t) {
	int d, s;
	double t1, d1, s1;
	
	printf ("%d ",t);
	FILE *fp;
	if (t==0) {
		fp = fopen(FILENAME, "w");
	}
	else {
		fp = fopen(FILENAME, "a");
	}
//	fprintf(fp, "IGOR\nwaves/d/O wave_time,wave_depth,wave_surface,wave_thermo\nbegin\n");
	for (d=0; d<Nd; d++) {
		for (s=0; s<Ns; s++) {
				t1 = t*ht;
				d1 = d*hd;
				s1 = s*hs;
			fprintf(fp, "%f,%f,%f,%f\n",t1,d1,s1,Tnew[d][s]);
		}
	}
	//	fprintf(fp, "end\n");
	fclose(fp);
	printf("-OK\n");
}

main(void) {
	int t,d,s;

// �����l�ݒ�
	generation(0);
	for (d=0; d<Nd; d++) {
		for (s=0; s<Ns; s++) {
			Tnew[d][s] = shoki(d,s) + ht*b*S[d][s]*(hd*hs);
		}
	}
	output(0);
	
	for (t=1; t<=Nt; t++) {
		for (d=0; d<Nd; d++) {
			for (s=0; s<Ns; s++) {
				Told[d][s] = Tnew[d][s];
			}
		}
		generation(t);
		keisu();
		gauss(Nd*Ns);
		replace();
		output(t);
	}
	
	return 0;
}

