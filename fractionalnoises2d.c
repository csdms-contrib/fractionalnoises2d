#ifdef HAVE_MALLOC_H
# include<malloc.h>
#endif
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1
#define PI 3.141592653589793

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(idum)
int *idum;
{
        static int iset=0;
        static float gset;
        float fac,r,v1,v2;
        float ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP


main()
{    float fact,*precyear,betax,sum,mean,sddev,sddevcomp;
     int lengthx,lengthy;
     int *nn,i,j,idum,index;
     char outputfile[30];
     FILE *fp;

     printf("\nName of output file: ");
     scanf("%s",outputfile);
     printf("\nLength of desired noise: (factor of two, please) ");
     scanf("%d",&lengthx);
     lengthx=lengthx*2;
     nn=ivector(1,2);
     nn[1]=lengthx;
     nn[2]=lengthx;
     lengthy=lengthx;
     printf("\nMean: ");
     scanf("%f",&mean);
     printf("\nStandard Dev: ");
     scanf("%f",&sddev);
     printf("\nBeta: ");
     scanf("%f",&betax);
     printf("\nRandom number seed: (any negative integer) ");
     scanf("%d",&idum);
     fp=fopen(outputfile,"w");
     precyear=vector(1,2*lengthx*lengthy);
     for (i=1;i<=lengthx;i++)
      for (j=1;j<=lengthy;j++)
       {precyear[2*(i-1)*lengthy+2*j-1]=gasdev(&idum);
        precyear[2*(i-1)*lengthy+2*j]=0.0;}
     fourn(precyear,nn,2,1);
     for (j=1;j<=lengthy;j++)
      {precyear[2*j-1]=0.0;
       precyear[2*j]=0.0;
       precyear[2*(j-1)*lengthy+1]=0.0;
       precyear[2*(j-1)*lengthy+2]=0.0;}
     for (i=1;i<=lengthx/2;i++)
      for (j=1;j<=lengthy/2;j++)
       {fact=pow(sqrt(i*i+j*j)/(2*lengthx),(betax+1)/2.0);
        index=2*i*lengthy+2*j+1;
        precyear[index]=precyear[index]/(lengthx*lengthy*fact);
        precyear[index+1]=precyear[index+1]/(lengthx*lengthy*fact);
        index=2*i*lengthy+2*lengthy-2*j+1;
        precyear[index]=precyear[index]/(lengthx*lengthy*fact);
        precyear[index+1]=precyear[index+1]/(lengthx*lengthy*fact); 
        index=2*lengthx*lengthy-2*(i-1)*lengthy-2*(lengthy-j)+1;
        precyear[index]=precyear[index]/(lengthx*lengthy*fact);
        precyear[index+1]=precyear[index+1]/(lengthx*lengthy*fact);
        index=2*lengthx*lengthy-2*lengthy-2*(i-1)*lengthy+2*(lengthy-j)+1;
        precyear[index]=precyear[index]/(lengthx*lengthy*fact);
        precyear[index+1]=precyear[index+1]/(lengthx*lengthy*fact);}
     fourn(precyear,nn,2,-1);
     sum=0.0;
     for (i=1;i<=lengthx/2;i++)
      for (j=1;j<=lengthy/2;j++)
       sum+=precyear[2*(i-1)*lengthy+2*j-1];
     sum=sum/(lengthx*lengthy);
     sddevcomp=0.0;
     for (i=1;i<=lengthx/2;i++)
      for (j=1;j<=lengthy/2;j++)
       sddevcomp+=(precyear[2*(i-1)*lengthy+2*j-1]-sum)*
        (precyear[2*(i-1)*lengthy+2*j-1]-sum);
     sddevcomp=sqrt(sddevcomp/(lengthx*lengthy));
     for (i=1;i<=lengthx/2;i++)
      for (j=1;j<=lengthy/2;j++)
       {index=2*(i-1)*lengthy+2*j-1;
        precyear[index]=(precyear[index]-mean)/sddevcomp;
        precyear[index]=(precyear[index]*sddev)+mean;         
        fprintf(fp,"%f\n",precyear[index]);}
     fclose(fp);
} 
