// �o�̓t�@�C��(�v���O�����{�̂Œ�`����ꍇ�̓R�����g�A�E�g���Ă�������)
//#define FILENAME "dat/all2.dat"

// ���Ԋi�q�_��
#define Nt 50
// ��Ԋi�q�_��
#define Nd 20
#define Ns 20

// ���Ԋi�q�_�Ԋu
#define ht 0.00001
// ��Ԋi�q�_�Ԋu
#define hd 0.1
#define hs 0.1


// ********** �����̐ݒ� **********
// ���x[g/mm3]
#define rho 0.00233
// ��M[J/g K]
#define cp 0.864
// �M�`����[W/mm K]
#define lambda 0.168

// �M�g�U��[mm2/s] > ���3�萔���猈��
//#define a

double Tnew[Nd][Ns];
double Told[Nd][Ns];
double S[Nd][Ns];


// **********�}�N���̂��߂̐ݒ�**********

// generation.h�̂��߂̐ݒ�
	// �����M��[W]=[J/s]
	#define p 1500
	
	//��������[0-1]
	#define pe 1

// **********�}�N���̐ݒ肱���܂�**********

// �������x�̐ݒ�
#include "shoki.h"
// �M�ʓ����̐ݒ�
#include "generation2.h"

double a;
double b;
void coefficient(int d, int s) {
// �w��̐[��,�\�ʈʒu�ɂ�����W�����w��
	b = 1/ (rho * cp);
	a = b * lambda;
}
