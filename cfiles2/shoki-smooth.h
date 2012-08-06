// 初期値設定用マクロ

double shoki(int d, int s) {
// t=0における初期値
	double func;
	func = 300;
	
	func = func + d;
	func = func + s;
	
	return func;
}
