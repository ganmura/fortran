// ガウスの消去法

void gauss(int N) {
	double f;
	int i,j,k;
	
	//前進過程
	for (i=0; i<N-1; i++) {
		f = A[i][i];
		for (j=i; j<N; j++) {
			A[i][j] = A[i][j]/f;
		}
		D[i] = D[i]/f;
		for (k=i+1; k<N; k++) {
			if (A[k][i]!=0) {
				f = A[k][i];
				for (j=i; j<N; j++) {
					A[k][j] = A[k][j] - f*A[i][j];
				}
				D[k] = D[k] - f*D[i];
			}
		}
	}
	f = A[N-1][N-1];
	A[N-1][N-1] = A[N-1][N-1]/f;
	D[N-1] = D[N-1]/f;
	
	//後退過程
	for (i=N-1; i>0; i--) {
		for (k=i-1; k>=0; k--) {
			if (A[k][i]!=0) {
				f = A[k][i];
				for (j=i; j>=0; j--) {
					A[k][j] = A[k][j] - f*A[i][j];
				}
				D[k] = D[k] - f*D[i];
			}
		}
	}
}
