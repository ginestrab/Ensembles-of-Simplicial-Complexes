/**************************************************************************************************
 * If you use this code, please cite G. Bianconi and O.T. Courtney
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
#define gamma2 2.8
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
/* Check if a tetrahedron exists. Takes four nodes as an input and outputs 1 if there already exists a triangle incident to them and 0 otherwise. */
int Check(i1,i2,i3,i4){
	int in,c=0;
	if ((i1==i2)||(i1==i3)||(i1==i4)||(i2==i3)||(i2==i4)||(i3==i4)){
		c=1;
	}
	if (c==0){
        for(in=0;in<kg[i1]-1;in++){
            if((tri[i1][0][in]==i2)&&(tri[i1][1][in]==i3)&&(tri[i1][2][in]=i4)){
                c=1;
                break;
            }
            if((tri[i1][0][in]==i2)&&(tri[i1][1][in]==i4)&&(tri[i1][2][in]=i3)){
                c=1;
                break;
            }
            if((tri[i1][0][in]==i3)&&(tri[i1][1][in]==i2)&&(tri[i1][2][in]=i4)){
                c=1;
                break;
            }
            if((tri[i1][0][in]==i3)&&(tri[i1][1][in]==i4)&&(tri[i1][2][in]=i2)){
                c=1;
                break;
            }
            if((tri[i1][0][in]==i4)&&(tri[i1][1][in]==i2)&&(tri[i1][2][in]=i3)){
                c=1;
                break;
            }
            if((tri[i1][0][in]==i4)&&(tri[i1][1][in]==i3)&&(tri[i1][2][in]=i2)){
                c=1;
                break;
            }
        }
	}
	return(c);
}

/*************************************************************************************************/
/* Create tetrahedron. Takes four nodes as an input and creates a tetrahedron incident to them. */
void Tetrahedron(int i1, int i2,int i3, int i4){
	int iaus;

	tri[i1][0]=(int*)realloc(tri[i1][0],kg[i1]*sizeof(int));
	tri[i1][1]=(int*)realloc(tri[i1][1],kg[i1]*sizeof(int));
	tri[i1][2]=(int*)realloc(tri[i1][2],kg[i1]*sizeof(int));

	tri[i2][0]=(int*)realloc(tri[i2][0],kg[i2]*sizeof(int));
	tri[i2][1]=(int*)realloc(tri[i2][1],kg[i2]*sizeof(int));
	tri[i2][2]=(int*)realloc(tri[i2][2],kg[i2]*sizeof(int));

	tri[i3][0]=(int*)realloc(tri[i3][0],kg[i3]*sizeof(int));
	tri[i3][1]=(int*)realloc(tri[i3][1],kg[i3]*sizeof(int));
	tri[i3][2]=(int*)realloc(tri[i3][2],kg[i3]*sizeof(int));

	tri[i4][0]=(int*)realloc(tri[i4][0],kg[i4]*sizeof(int));
	tri[i4][1]=(int*)realloc(tri[i4][1],kg[i4]*sizeof(int));
	tri[i4][2]=(int*)realloc(tri[i4][2],kg[i4]*sizeof(int));

	tri[i1][0][kg[i1]-1]=i2;
	tri[i1][1][kg[i1]-1]=i3;
	tri[i1][2][kg[i1]-1]=i4;

	tri[i2][0][kg[i2]-1]=i1;
	tri[i2][1][kg[i2]-1]=i3;
	tri[i2][2][kg[i2]-1]=i4;

	tri[i3][0][kg[i3]-1]=i1;
	tri[i3][1][kg[i3]-1]=i2;
	tri[i3][2][kg[i3]-1]=i4;

	tri[i4][0][kg[i4]-1]=i1;
	tri[i4][1][kg[i4]-1]=i2;
	tri[i4][2][kg[i4]-1]=i3;
}
/*************************************************************************************************/

int main(int argc, char** argv){
	int i,j,j2,i1,i2,i3,i4,naus,*knng,*pkg,*k,**l,*pk,*knn,n,**a,*Ck;
	double xaus, x;
	char filec[60];

	FILE *gp;

	srand48(time(NULL));

	gp=fopen("edge_list.txt","w");

	kgi=(int*)calloc(N,sizeof(int));
	kg=(int*)calloc(N,sizeof(int));
	k=(int*)calloc(N,sizeof(int));
	a=(int**)calloc(N,sizeof(int*));
	knng=(int*)calloc(N,sizeof(int));
	pkg=(int*)calloc(N,sizeof(int));
	knn=(int*)calloc(N,sizeof(int));
	pk=(int*)calloc(N,sizeof(int));
	Ck=(int*)calloc(N,sizeof(int));
	tri=(int***)calloc(N,sizeof(int**));
	for(i=0;i<N;i++){
            a[i]=(int*)calloc(N,sizeof(int));
            tri[i]=(int**)calloc(3,sizeof(int*));
            tri[i][0]=NULL;
            tri[i][1]=NULL;
            tri[i][2]=NULL;
    }

	xaus=4;
	while(xaus>3){
	/***********************************************************************************************/
	/* Initialization */
        for(i=0;i<N;i++){
        /* Nodes are assigned desired generalized degree according to a scale-free distribution */
            kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
            /* kgi[i]= poisson(lambda); */
            while(kgi[i]>pow(N,3.)/6){
            /* Desired generalized degrees are re-drawn if they exceed the maximum possible generalized degree of a node (natural cut-off) */
                kgi[i]=(int)(m*pow(drand48(),-1./(gamma2-1.)));
                /* kgi[i]= poisson(lambda); */
            }
            kg[i]=0;  /* Generalized degree of node iintially set to 0 */
            k[i]=0;   /* Degree of node i intially set to 0 */
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
        while((xaus>4)&&(naus<1+Avoid*NX)){
        /* Randomly select four nodes proportional to the number of unmatched stubs they have remaining. */
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

            x=xaus*drand48();
            i3=Choose(x);
            kg[i3]++;
            kgi[i3]--;
            xaus--;

            x=xaus*drand48();
            i4=Choose(x);
            kg[i4]++;
            kgi[i4]--;
            xaus--;

            /* Check proposed matching is legal */
            if((Check(i1,i2,i3,i4)==0)){
            /* Proposed matching legal. Create tetrahedron and links */
                Tetrahedron(i1,i2,i3,i4);

                a[i1][i2]=1;
                a[i2][i1]=1;

                a[i1][i3]=1;
                a[i3][i1]=1;

                a[i1][i4]=1;
                a[i4][i1]=1;

                a[i2][i3]=1;
                a[i3][i2]=1;

                a[i2][i4]=1;
                a[i4][i2]=1;

                a[i3][i4]=1;
                a[i4][i3]=1;
            }
            else{
            /* Proposed matching illegal. Back-track and increment back-track counter by one */
                naus++;
                if(Avoid==1){
                    kg[i1]--;
                    kgi[i1]++;
                    kg[i2]--;
                    kgi[i2]++;
                    kg[i3]--;
                    kgi[i3]++;
                    kg[i4]--;
                    kgi[i4]++;
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


