/*
非定常2次元熱伝導
custom-crank2.c
2004.6/26頃	とりあえず:custom-crank.c
2004.10/4	一部改正(条件ファイルなども含め):custom-crank2.c

差分法(クランク-ニコルソン法)を使用
熱発生も考慮
逐次出力型
*/

#include <stdio.h>
#include "joken2.h"

// 出力ファイル(joken.hの中で定義する場合はコメントアウトしてください)
#define FILENAME "dat/crank-all2.dat"

double A[Nd*Ns][Nd*Ns];
double D[Nd*Ns];

#include "gauss.h"


void keisu(void) {
// 行列計算用A,D作成
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
			
			D[Nd*i+j] = 2*(1 - Fd - Fs)*Told[i][j] + 2*ht*b*S[i][j]*(hd*hs);
			switch (i) {
				case 0:
					D[Nd*i+j] = D[Nd*i+j] + 2*Fd*Told[i+1][j];
					break;
				case Nd-1:
					D[Nd*i+j] = D[Nd*i+j] + 2*Fd*Told[i-1][j];
					break;
				default:
					D[Nd*i+j] = D[Nd*i+j] + Fd*(Told[i-1][j] + Told[i+1][j]);
			}
			switch (j) {
				case 0:
					D[Nd*i+j] = D[Nd*i+j] + 2*Fs*Told[i][j+1];
					break;
				case Ns-1:
					D[Nd*i+j] = D[Nd*i+j] + 2*Fs*Told[i][j-1];
					break;
				default:
					D[Nd*i+j] = D[Nd*i+j] + Fs*(Told[i][j-1] + Told[i][j+1]);
			}
			
			A[Nd*i+j][Nd*i+j] = 2*(1 + Fd + Fs);
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
			Tnew[i][j] = D[Nd*i+j];
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

// 初期値設定
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

