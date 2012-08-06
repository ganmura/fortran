// 出力ファイル(プログラム本体で定義する場合はコメントアウトしてください)
//#define FILENAME "dat/all2.dat"

// 時間格子点数
#define Nt 50
// 空間格子点数
#define Nd 20
#define Ns 20

// 時間格子点間隔
#define ht 0.00001
// 空間格子点間隔
#define hd 0.1
#define hs 0.1


// ********** 物質の設定 **********
// 密度[g/mm3]
#define rho 0.00233
// 比熱[J/g K]
#define cp 0.864
// 熱伝導率[W/mm K]
#define lambda 0.168

// 熱拡散率[mm2/s] > 上の3定数から決定
//#define a

double Tnew[Nd][Ns];
double Told[Nd][Ns];
double S[Nd][Ns];


// **********マクロのための設定**********

// generation.hのための設定
	// 投入熱量[W]=[J/s]
	#define p 1500
	
	//投入効率[0-1]
	#define pe 1

// **********マクロの設定ここまで**********

// 初期温度の設定
#include "shoki.h"
// 熱量投入の設定
#include "generation2.h"

double a;
double b;
void coefficient(int d, int s) {
// 指定の深さ,表面位置における係数を指定
	b = 1/ (rho * cp);
	a = b * lambda;
}
