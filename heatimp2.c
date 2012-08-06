/* 
************************************************************************
*  heat.c:   差分法を用いた熱伝導方程式の解法 （陰解法）               *
************************************************************************
*/

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define size 11			/* 格子数 */

main()
{ 
   int i,j,max,jj;
   double pi, h, r, k, x, dx, dt, u[size][2];     
   double uu[size], a1[size], a2[size], a3[size], a4[size], s[size], p[size];
   
   FILE *output;			/* heat.datにデータを保存 */
   output = fopen("heat.dat","w");
   
   pi=10;

   max = 900;
   k=1.0;
   h=2.0/pi;
   dx=pi/(size-1);
   r  = k * dt / dx*dx;
   dt=3.0*dx*dx/k;

   for(j=0; j<size; j++){
        x=dx*j;
        if(x < pi/2.0){
          u[j][0]=h*x;
	      }else{
	        u[j][0]=h*(pi-x);
	      }
      }

   for(i=0; i<2; i++){ /* 境界条件*/
      u[0][i]      = 0.;
      u[size-1][i] = 0.;
   } 
   for(j=0 ; j<size; j++){
       fprintf(output, "%f\n", u[j][0]);
   }
   fprintf(output, "\n");			/* gnuplotのための空行 */
   
   for(i=1; i<=max; i++){ /* max時間ステップにわたってループ */
     for( j=0; j<=size-1; j++){
       a1[j]= r;
       a2[j]=-2.*r-1.;
       a3[j]=r;
       if(j==0){
	        a4[j]=(2*r-1)*u[j][0]+r*u[j+1][0];
       }else{
	        if(j==size-1){
	          a4[j]=r*u[j-1][0]+(2*r-1)*u[j][0];
	        }else{
	          a4[j]=r*u[j-1][0]+(2*r-1)*u[j][0]+r*u[j+1][0];
	        }
        }
      }
      //
      p[0]=a2[0];
      s[0]=a4[0];
      for( j=1; j<=(size-1); j++){
        p[j]=a2[j]-a1[j]*a3[j-1]/p[j-1];
        s[j]=a4[j]-a1[j]*s[j-1]/p[j-1];
      }
      //
      uu[size-1]=s[size-1]/p[size-1];
      for( j = 1; j <= (size-1 - 1); j++ ){
       jj=size-1-j;
       uu[jj] = (s[jj]-a3[jj]*uu[jj+1])/p[jj];
      }
      for( j = 1; j <= (size-1 - 1); j++ ){
        u[j][1] = uu[j];
      }

      for(j=0 ; j<size; j++){
	      fprintf(output, "%f\n", u[j][1]);
      }
      fprintf(output, "\n");			/* gnuplotのための空行 */
      for(j=0; j<size; j++) u[j][0]=u[j][1];    /* 新しい時間を古い時間に */
   }
   printf("data stored in heat.dat\n");
   fclose(output);
}      

