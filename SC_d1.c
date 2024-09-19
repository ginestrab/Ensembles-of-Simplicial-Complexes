/**************************************************************************************************
 * If you use this code, please cite  
 * O.T. Courtney and G. Bianconi 
 * "Generalized network structures: the configuration model and the canonical ensemble of
 * simplicial complexes"
 * Phys. Rev. E 93, 062311 (2016)
***************************************************************************************************
 * Code that  generates random simplicial complexes with scale-free generalized degree 
 * distribution.
 *
 * The option to use a Poisson distributed generalized degree distribution has also been included.  
 * The necessary code may be found in comments at the relevant points.
 *
 * This code uses:
 * N  Number of nodes in the simplicial complex
 * m  The minimum of the scale-free distribution
 * gamma2  Exponent of the scale-free distribution
 * lambda  Expected value of the Poisson distribution (commented-out)
 * Avoid  Whether or not 'back-tracking' is allowed when illegal matchings are proposed
 * (Avoid==1 allowed, Avoid==0 not allowed)
 * NX  Maximum number of 'back-tracks' before matching process restarts from an unmatched network
 *************************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 10000
#define m 1
#define gamma2 2.3
/* #define lambda 10 */
#define Avoid 1
#define NX 80

int *kgi,*kg,***tri;

/*************************************************************************************************/
/* Randomly select an unmatched stub. Choose takes as its input a random number between 0 and the 
total number of stubs and gives as its output the index of the node of the selected stub */
int Choose(double x){
	int i1,i;
	for (i=0;i<N;i++){
		x-=kgi[i];
		if (x<0){
			i1=i;
			break;
		}
	}
	return(i1);
}
/*************************************************************************************************/

int main(int argc, char** argv){
	int i,j,nrun,j2,i1,i2,i3,naus,*knng,*pkg,*k,**l,*pk,*knn,n,**a,*Ck;
	double xaus, x;
	char filec[60];

	FILE *fp,*gp;

	gp=fopen("edge_list.txt","w");
	srand48(time(NULL));
	kgi=(int*)calloc(N,sizeof(int));
	kg=(int*)calloc(N,sizeof(int));
	k=(int*)calloc(N,sizeof(int));
	a=(int**)calloc(N,sizeof(int*));
	knng=(int*)calloc(N,sizeof(int));
	pkg=(int*)calloc(N,sizeof(int));
	knn=(int*)calloc(N,sizeof(int));
	pk=(int*)calloc(N,sizeof(int));
	Ck=(int*)calloc(N,sizeof(int));

	for(i=0;i<N;i++){
		a[i]=(int*)calloc(N,sizeof(int));
	}

	xaus=4;  

	while(xaus>2){
	/***********************************************************************************************/
	/* Initialization */
		for(i=0;i<N;i++){
		/* Nodes are assigned desired generalized degree according to a scale-free distribution */
			kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
			/* kgi[i]= poisson(lambda); */
			while(kgi[i]>(N-1)){
			/* Desired generalized degrees are re-drawn if they exceed the maximum possible generalized degree of a node (natural cut-off) */
				kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
				/* kgi[i]= poisson(lambda); */
			}
			kg[i]=0;  /* Generalized degree of node i intially set to 0 */
			k[i]=0;  /* Degree of node i intially set to 0 */
			for(j=0;j<N;j++){
				a[i][j]=0;
			}
		}
		xaus=0;
		for(i=0;i<N;i++){
			xaus+=kgi[i];
		}
		naus=0; /* Back-track counter initially set to zero */
	/***********************************************************************************************/
	/* Stubs matched */
		while((xaus>3)&&(naus<1+Avoid*NX)){
			/* Randomly select two nodes proportional to the number of unmatched stubs they have remaining. */
			x=xaus*drand48();
			i1=Choose(x);
			kg[i1]++;
			kgi[i1]--;
			xaus--;

			x=xaus*drand48();
			i2=Choose(x);
			kg[i2]++;
			kgi[i2]--;
			xaus--;

			/* Check proposed matching is legal */
			if((i1!=i2)&&(a[i1][i2]==0)){
				/* Proposed matching legal. Create link */
				a[i1][i2]=1;
				a[i2][i1]=1;
			}
			else{
				/* Proposed matching illegal. Back-track and increment back-track counter by one */
				naus++;
				if(Avoid==1){
					kg[i1]--;
					kgi[i1]++;
					kg[i2]--;
					kgi[i2]++;
				}
			}
		}
	}
/*************************************************************************************************/
/* Degrees calculated */
	for (i=0;i<N;i++){
		for(j=i+1;j<N;j++){
			if(a[i][j]>0){
				k[i]++;
				k[j]++;
			}
		}
	}
/*************************************************************************************************/
/* Print list of edges to file */
	for (i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(a[i][j]==1){
				fprintf(gp,'%d %d\n',i,j)
			}
		}

	}
/*************************************************************************************************/
	fclose(gp);

	return 0;
}


