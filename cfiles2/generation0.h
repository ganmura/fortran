// 発生熱量
// S : 単位体積･単位時間あたりの熱発生率[W/mm^2]

void generation(int t) {
	int i,j;
	
	// Sリセット
	for (i=0; i<Nd; i++) {
		for (j=0; j<Ns; j++) {
			S[i][j] = 0;
		}
	}
	
	// 熱発生率指定

}
