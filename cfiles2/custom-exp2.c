/*
非定常2次元熱伝導
custom-exp2.c
2004.6/26頃	とりあえず:custom-exp.c
2004.10/4	一部改正(条件ファイルなども含め):custom-exp2.c

差分法(陽解法)を使用
熱発生も考慮
逐次出力型
*/

#include <stdio.h>
#include "joken2.h"

// 出力ファイル(joken.hの中で定義する場合はコメントアウトしてください)
#define FILENAME "dat/exp-all2.dat"


void thermo(void) {
	int i, j;
	double Fd, Fs;
	double T1, T2;
	
	for (i=0; i<Nd; i++) {
		for (j=0; j<Ns; j++) {
			coefficient(i,j);
			Fd = a*ht/(hd*hd);
			Fs = a*ht/(hs*hs);
			
			switch (i) {
				case 0:
					T1 = 2*Told[i+1][j];
					break;
				case Nd-1:
					T1 = 2*Told[i-1][j];
					break;
				default:
					T1 = Told[i-1][j] + Told[i+1][j];
			}
			switch (j) {
				case 0:
					T2 = 2*Told[i][j+1];
					break;
				case Ns-1:
					T2 = 2*Told[i][j-1];
					break;
				default:
					T2 = Told[i][j-1] + Told[i][j+1];
			}
			Tnew[i][j] = Fd*T1 + Fs*T2 + (1-2*(Fd+Fs))*Told[i][j] + ht*b*S[i][j]*(hd*hs);
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

int shusoku(void) {
	coefficient(0,0);
	if (ht <= (hd*hd)/(2*a) && ht <= (hs*hs)/(2*a)) {
		printf("Convergence-condition : OK!\n");
		return 1;
	}
	else {
		printf("Convergence-condition : Bad & STOP Program\n");
		return 0;
	}
}

main(void) {
	int t,d,s;
	int i;

// 収束条件判定
	i = shusoku();
	if (i==0) {
		return 0;
	}
	
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
		thermo();
		output(t);
	}
	
	return 0;
}

