// This code perfoms the steepest descent-like iterative search algorithm within the space of
// rotation symmetric S-boxes (RSSBs), 2-RSSBs, 5-RSSBs, and the concatenation of RSSBs in dimension 10.

// Selçuk Kavut
// Dept. of Computer Engineering
// Balýkesir University, Turkey

// The code is compiled by Microsoft Visual C++ 2010.

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
//The files "mersenne.cpp" and "randomc.h" are used to generate random S-box.
#include "mersenne.cpp" 
#include "randomc.h"  

//n: Number of variables
#define n 10
#define N 1024

int cek=0,NL;
long double FW21[N+1],FW22[N+1];
//NIT: Number of iterations.
//NR : Number of runs.
int NIT,NR,SZ[N],h,gn;
//PR: Contains the permutations corresponding to the symmetric S-box.
signed char PR[4][N]={{1,2,3,4,5,6,7,8,9,0},	//permutation for RSSBs
					  {0,2,3,4,5,6,7,8,9,1},	//permutation for the concatenation of RSSBs
					  {2,3,4,5,6,7,8,9,0,1},	//permutation for 2-RSSBs
					  {5,6,7,8,9,0,1,2,3,4}		//permutation for 5-RSSBs
};
signed char p[n];
signed char BS[N];

int *R   =(int *)malloc(N*N*sizeof(int));
int *ROT =(int *)malloc(N*N*sizeof(int));
int *DBP =(int *)malloc(N*n*sizeof(int));
int *FWec=(int *)malloc(N*(N-1)*sizeof(int));
int *Sbc =(int *)malloc(N*n*sizeof(int));
int *Bc  =(int *)malloc(N*n*sizeof(int));
int *G   =(int *)malloc(N*n*sizeof(int));
int *Rb  =(int *)malloc(N*n*sizeof(int));
signed char *A  =(signed char *)malloc(N*N*sizeof(signed char));
signed char *HD =(signed char *)malloc(N*N*2*sizeof(signed char));

//The file "H10.txt" contains the Walsh-Hadamard matrix.
FILE *outh=fopen("H10.txt", "r");
//The file "ANF10.txt" is used to compute algebraic degree.
FILE *outa=fopen("ANF10.txt", "r");
//Result are written to the file "S10_P.txt".
FILE *outr=fopen("S10_P.txt", "w");

int _tmain(int argc, _TCHAR* argv[])
{
	long double fact(int k);
	void repres();
	void orbit(int k);
	int muk(const void *x,const void *y);
	void init ();
	void RSBOX(int *S);
	void dec2bin(int d, int *b);
	void ACOR(int *FW, int *AC); 
	int anf(signed char *T);
	void walsh(signed char *T, int *FW);
	long double sumsse(int *FW);
	long double sumssea(int *FWxt);
	int findmaxor(int *DUA);
	int findmingl(int *tt);
	int findmaxgl(int *tt);
	int findmaxac(int *tt);
	int findmaxwh(int *tt);
	int rot2d(int d);
	int rotSch(int *S);
	void S2T(int *S, unsigned int *ts);

	time_t rawtime;
	struct tm * timeinfo;

	signed char x,T[N];
	int i,j,k,tmp,wxor,KUP,cnt,K,it,jt,kt,t0,k1,t1,jx,Cnt,indx,CHK,NLS,ACS,ADS,SS[N];
	unsigned int ts[N/2];
	int S[N],FW[N],AC[N],NLs[N-1],ACs[N-1],ADs[N-1];
	long double I,Maxi,LI;
	int r=0,Sc[N],MO,CS,PS,ST[N];

	int Cx=0,Rx[N],wxa=N,NLA=-1,ACA=N;
	int *INS=(int *)malloc(N*(N-1)*sizeof(int));
	int *SZS=(int *)malloc(N*(N-1)*sizeof(int));

	int *STF=(int *)malloc(N*2*sizeof(int));
	int *FWxt=(int *)malloc(N*(N-1)*sizeof(int));
	int *DUA=(int *)malloc(N*N*sizeof(int));
	int *FWex=(int *)malloc(N*(N-1)*sizeof(int));
	int *CF=(int *)malloc((N-1)*n*sizeof(int));
	int *DB=(int *)malloc(N*n*sizeof(int));
	int *Sb=(int *)malloc(N*n*sizeof(int));
	signed char *SR=(signed char *)malloc(N*(N-1)*sizeof(signed char));
	signed char *HDM=(signed char *)malloc(N*N*sizeof(signed char));
	signed char *HD2=(signed char *)malloc(N*N*2*sizeof(signed char));
	unsigned int *Q3=(unsigned int *)malloc((N/2)*1000*sizeof(unsigned int));
	int *MEM=(int *)malloc(N*(N-1)*5*sizeof(int));
	long double *COST=(long double *)malloc(N*N*sizeof(long double));

	printf("\n# of runs: ");
	scanf("%d",&NR);
	printf("\n# of iterations: (choose less than 1000): ");
	scanf("%d",&NIT);
	printf("\nSelect permutation (1: RSSBs, 2: Concatenations, 3: 2-RSSBs, 4: 5-RSSBs): ");
	scanf("%d",&PS);
	if (PS>4 || PS<1 || NIT>1000 || NR<=0)
		return 0;
	printf("\n# of runs = %d \n#of iterations = %d",NR,NIT);
	printf("\n(Permutation )p = ");
	PS=PS-1;
	for (i=0;i<n;i++)
		printf("%d ",PR[PS][i]);

	for (i=0;i<n;i++)
		p[i]=PR[PS][i];
	fprintf(outr,"\n# of runs = %d \n#of iterations = %d",NR,NIT);
	fprintf(outr,"\n(Permutation )p = ");
	for (i=0;i<n;i++)
		fprintf(outr,"%d ",PR[PS][i]);
	
	//repres(): Generates ortbits and their sizes.
	repres();

	printf("\n(Total # of orbits) gn=%d",gn);
	fprintf(outr,"\n(Total # of orbits) gn=%d",gn);

	for (i=0;i<N;i++)
		ST[i]=0;
	for (i=0;i<gn;i++)
		ST[SZ[i]]=ST[SZ[i]]+1;
	j=0;
	for (i=0;i<gn;i++)
		if (ST[i]!=0)
		{
			*(STF+j*2)=i;
			*(STF+j*2+1)=ST[i];
			j=j+1;
		}
	K=0;
	for (i=0;i<j;i++)
	{
		printf("\nOrbit size:%d...number of them:%d",*(STF+i*2),*(STF+i*2+1));
		fprintf(outr,"\nOrbit size:%d...number of them:%d",*(STF+i*2),*(STF+i*2+1));
		K=K+(*(STF+i*2))*(*(STF+i*2+1));
	}

	if (K!=N)
	{
		printf("\nError K=%d",K);
		return 0;
	}

	LI=0;
	for (i=0;i<j;i++)
	{
		I=(long double) (pow((long double) (*(STF+i*2)),(long double) (*(STF+i*2+1))));
		LI=LI+ log(I)/log(2.0)+fact(*(STF+i*2+1)); 
	}
	printf("\nTotal # of these symmetric (and bijective) S-boxes=2^%f",LI);
	fprintf(outr,"\nTotal # of these symmetric (and bijective) S-boxes=2^%f",LI);

	// ROT: Contains all the permutations of every input of an S-box.
	for (i=0;i<N;i++)
	{
		k=i;
		for (j=0;j<N;j++)
		{
			*(ROT+i*N+j)=k;
			k=rot2d(k);
		}
	}
	
	//HDM: Walsh-Hadamard matrix.
	i=0;j=0;
	while (!feof(outh))
	{
		fscanf(outh,"%d ",&x);

		*(HDM+i*N+j)=x;
		j=j+1;
		if (j==N)
		{
			i=i+1;
			j=0;
		}
	}
	fclose(outh);
	
	//A: Used to compute the algebraic degree.
	i=0;j=0;
	while (!feof(outa))
	{
		fscanf(outa,"%d ",&x);

		*(A+i*N+j)=x;
		j=j+1;
		if (j==N)
		{
			i=i+1;
			j=0;
		}
	}
	fclose(outa);
	
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
		{
			*(HD+i*N*2+j*2)  =  *(HDM+i*N+j);
			*(HD+i*N*2+j*2+1)=-(*(HDM+i*N+j));
		}
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
		{
			*(HD2+i*N*2+j*2)  = 2*(*(HDM+i*N+j));
			*(HD2+i*N*2+j*2+1)=-2*(*(HDM+i*N+j));
		}

	// FW22: Pre-computed values used to efficiently compute the cost function.
	I=0;
	for (i=0;i<N+1;i++)
	{
		FW21[i]=1;
		for (j=0;j<2;j++)
			FW21[i]=I*FW21[i];
		I=I+1;
	}
	for (i=0;i<N+1;i++)
		FW22[i]=(FW21[i]-N)*(FW21[i]-N);

	//DB: Binary values of all the inputs of an S-box.
	for (i=0;i<N;i++)
		for (j=0;j<n;j++)
			*(DB+i*n+j)=(i&(1<<(n-1-j)))>>(n-1-j);

	//DBP: Polar form of DB.
	for (i=0;i<N;i++)
		for (j=0;j<n;j++)
			if (*(DB+i*n+j)==1)
				*(DBP+i*n+j)=-1;
			else
				*(DBP+i*n+j)=1;

	//BS: Hamming weights of all the inputs of an S-box.
	for (i=0;i<N;i++)
	{
		BS[i]=0;
		for (j=0;j<n;j++)
			BS[i]=BS[i]+*(DB+i*n+j);
	}

	//CF: Binary values of the non-zero orbit representatives.
	for (i=1;i<gn;i++)
		for (j=0;j<n;j++)
			*(CF+(i-1)*n+j)=(*(R+i*N)&(1<<(n-1-j)))>>(n-1-j);

	//CS: Number of all possible pairs of orbits with the same orbit size.
	CS=0;
	for (i=0;i<gn-1;i++)
		for (j=i+1;j<gn;j++)
			if (SZ[i]==SZ[j])
			{
				*(INS+CS*2)=i;
				*(INS+CS*2+1)=j;
				SZS[CS]=SZ[i];
				CS=CS+1;
			}

	//Cx: Number of orbits with orbit size > 1.
	for (i=0;i<gn;i++)
		if (SZ[i]>1)
		{
			Rx[Cx]=i;
			Cx=Cx+1;
		}

	for (KUP=1;KUP<=NR;KUP++)
	{
	
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	printf ("\n(KUP=%d) Current local time and date: %s", KUP,asctime(timeinfo));
	fprintf (outr,"\n(KUP=%d) Current local time and date: %s", KUP,asctime(timeinfo));

	fclose(outr);
	outr=fopen("S10_P.txt", "a");

	//RSBOX: Generates an S-box randomly.
	RSBOX(S);

	//Sx: Inverse function
	//Uncommenting the following lines, the initial S-box can be chosen as the inverse function. In this case, the permutation should be selected accordingly (that is, PS can be 1, 3, or 4).
	//int Sx[]={0,970,917,525,811,321,27,392,599,300,642,417,54,635,784,609,175,87,600,404,261,80,834,1001,108,863,247,6,545,938,195,709,350,627,174,876,177,209,808,155,522,254,160,138,645,352,979,840,216,896,703,480,494,568,12,379,67,852,853,937,390,178,395,419,700,453,231,56,348,533,729,646,354,117,418,258,593,322,310,699,21,344,508,619,320,206,276,17,267,542,704,496,935,167,657,569,432,411,769,726,383,772,960,554,988,901,113,355,24,993,758,361,134,106,681,190,683,73,851,1007,780,748,356,867,790,663,838,532,377,549,906,1003,462,799,112,618,696,654,43,835,435,752,269,605,708,849,234,919,836,507,516,486,163,191,644,39,620,932,375,203,42,550,688,152,1016,990,215,93,640,987,412,425,552,909,34,16,534,36,61,498,385,405,992,221,847,819,334,286,291,1000,115,153,864,281,822,30,515,742,429,907,766,505,521,159,897,975,85,740,953,37,779,255,226,270,710,166,48,312,963,776,493,183,722,702,268,868,212,337,339,930,380,66,343,295,146,251,679,365,991,695,537,954,473,587,712,391,711,26,557,601,303,235,653,898,41,211,754,336,75,926,789,20,983,746,924,308,575,88,224,142,213,556,369,389,285,592,86,563,647,809,870,193,481,504,538,274,187,933,393,967,675,188,468,894,815,233,649,759,1014,483,9,357,972,250,326,697,382,1011,265,579,78,813,217,744,841,347,750,865,406,736,84,5,77,547,353,908,304,580,1009,735,957,345,430,912,186,472,257,227,951,228,824,682,850,232,81,331,795,315,68,928,32,623,45,324,72,107,122,301,996,871,770,111,810,561,961,237,442,478,671,272,615,610,668,643,572,158,582,128,977,55,230,446,306,100,705,180,562,768,621,273,60,245,7,288,461,62,858,428,791,492,509,488,1010,900,19,181,318,803,771,915,927,97,170,415,457,413,883,11,74,63,535,444,510,427,452,171,540,423,397,198,332,902,96,775,624,140,903,948,529,1006,986,745,366,887,421,889,381,658,536,638,713,727,424,65,674,1022,678,414,837,845,760,394,132,641,686,626,590,882,292,526,502,905,335,242,730,958,959,632,367,749,51,282,885,299,946,512,151,757,401,839,782,974,399,220,52,594,91,920,179,1002,606,761,470,731,283,201,773,149,82,400,422,739,485,774,672,196,150,720,829,816,555,202,40,1012,943,3,469,866,825,438,616,589,127,69,176,420,448,240,284,701,426,980,89,721,738,28,778,323,570,129,161,861,172,821,103,520,271,248,595,796,717,363,386,277,962,689,1008,692,53,95,548,1015,374,945,843,266,786,1013,911,309,327,929,376,814,936,971,765,243,607,531,466,613,275,76,495,558,1005,724,966,8,18,249,714,622,921,143,500,588,652,15,371,965,764,591,999,370,530,639,135,83,156,388,603,351,434,680,465,33,659,637,694,859,477,805,707,13,812,629,449,617,168,463,10,373,154,44,71,278,706,296,793,916,608,252,137,978,995,94,447,628,891,753,690,125,860,1017,801,918,372,685,944,368,514,785,454,290,879,684,456,236,625,114,341,116,677,669,464,823,162,565,662,947,567,792,630,239,136,305,833,79,64,539,223,50,90,384,648,634,144,31,214,246,244,450,602,913,969,560,719,718,517,543,222,725,597,723,99,451,899,70,474,503,884,955,956,329,319,875,544,511,207,934,197,832,313,441,263,964,121,479,316,886,141,661,256,890,931,487,110,297,460,501,892,877,612,586,200,881,387,98,360,408,101,506,513,433,219,806,546,210,120,862,490,872,14,673,576,942,922,260,124,398,693,650,856,346,559,1019,984,133,1018,666,976,407,997,633,777,818,38,279,362,4,636,311,583,294,519,994,807,185,831,553,194,687,340,528,830,941,914,518,826,820,743,698,22,139,148,458,126,489,47,314,888,574,1020,459,854,184,904,145,342,118,57,58,846,923,794,985,396,631,664,551,781,25,192,317,527,123,225,968,280,359,783,874,873,737,35,763,989,676,949,767,467,416,732,482,751,443,842,445,755,660,762,950,293,952,49,204,253,728,403,105,431,436,848,471,130,199,325,173,1021,578,333,715,828,409,651,2,667,147,497,604,788,855,264,982,259,410,349,581,229,756,157,287,741,92,584,59,29,973,1004,827,787,524,670,573,484,691,437,880,893,338,895,208,241,733,734,330,475,476,102,364,564,218,747,611,598,289,869,716,1,585,302,939,491,205,802,378,655,46,541,998,925,262,798,857,440,169,104,878,165,238,182,109,817,656,358,804,981,614,189,23,499,131,940,596,439,119,566,328,402,307,523,577,298,571,164,665,800,797,844,910,455,1023};
	//for (i=0;i<N;i++)
	//	S[i]=Sx[i];

	fprintf(outr,"\nSbox=[");
	for (i=0;i<N;i++)
		fprintf(outr,"%d,",S[i]);
	fprintf(outr,"];");
	
	//Sb: Binary version of S.
	for (i=0;i<N;i++)
		for (j=0;j<n;j++)
			*(Sb+i*n+j)=*(DB+S[i]*n+j);

	//SR: Component functions obtained the by linear combinations corresponding to the orbit representatives.
	for (i=0;i<gn-1;i++)
		for (j=0;j<N;j++)
			*(SR+i*N+j)=0;

	for (i=0;i<gn-1;i++)
	{
		for (j=0;j<N;j++)
		{
			for (k=0;k<n;k++)
				*(SR+i*N+j)=(*(SR+i*N+j))^((*(Sb+j*n+k))*(*(CF+i*n+k)));
			T[j]=*(SR+i*N+j);
		}
		//walsh: Computes Walsh-Hadamard spectrum.
		walsh(T,FW);
		for (j=0;j<N;j++)
			*(FWex+j*(N-1)+i)=FW[j];
		//ACOR: Computes autocorrelation spectrum.
		ACOR(FW,AC);
		//NLs: Nonlinearities of the component functions.
		NLs[i]=N/2-findmaxwh(FW)/2;
		//ACS: Absolute indicators of the component functions.
		ACs[i]=findmaxac(AC);
		//ADs: Algebraic degrees of the component functions.
		ADs[i]=anf(T);
	}

	for (i=0;i<N;i++)
		for (j=0;j<gn-1;j++)
			*(FWxt+i*(N-1)+j)=*(FWex+i*(N-1)+j);

	//DUA: Difference distribution table.
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			*(DUA+j*N+i)=0;

	for (i=1;i<N;i++)
		for (j=0;j<N;j++)
		{
			tmp=S[j^i]^S[j];
			*(DUA+tmp*N+i)=*(DUA+tmp*N+i)+1;
		}

	//wxor: Differential uniformity.
	wxor=findmaxor(DUA);
	//NLS: Nonlinearity.
	NLS=findmingl(NLs);
	//ACS: Absolute indicator.
	ACS=findmaxgl(ACs);
	//ADS: Algebraic degree.
	ADS=findmaxgl(ADs);
	
	printf("\n\nStart KUP=%d NL=%d AC=%d AD=%d XOR=%d COST=%f",KUP,findmingl(NLs),findmaxgl(ACs),findmaxgl(ADs),wxor,sumssea(FWxt));
	fprintf(outr,"\n\nStart KUP=%d NL=%d AC=%d AD=%d XOR=%d COST=%f",KUP,findmingl(NLs),findmaxgl(ACs),findmaxgl(ADs),wxor,sumssea(FWxt));
	fclose(outr);
	outr=fopen("S10_P.txt", "a");

	//Q3: Stores all the iteration outputs.
	S2T(S,ts);
	cnt=0;
	for (i=0;i<N/2;i++)
	{
		*(Q3+cnt)=ts[i];
		cnt=cnt+1;
	}

	//Iterations are performed by the following loop.
	for (K=1;K<=NIT;K++)
	{
		Cnt=0;
		Maxi=1.6e+100;
		//The following loop generates the neighours obtained by swapping two orbits having the same size.
		//Then the corresponding cost is computed.
		for (i=0;i<CS;i++) 
		{
			MO=SZS[i];
			for (it=0;it<gn-1;it++) 
				for (jt=0;jt<MO;jt++)
				{

					t0=SR[it*N+*(R+*(INS+i*2)*N  +jt) ];
					t1=SR[it*N+*(R+*(INS+i*2+1)*N+jt) ];
					if (t0 != t1)
					{
						k1=*(R+*(INS+i*2)*N+jt);
						for (jx=0;jx<N;jx++)
							FWxt[jx*(N-1)+it]=FWxt[jx*(N-1)+it]+*(HD2+k1*N*2+jx*2+t1);
					}

					t0=SR[it*N+*(R+*(INS+i*2+1)*N+jt) ];
					t1=SR[it*N+*(R+*(INS+i*2)*N  +jt) ];
					if (t0 != t1)
					{
						k1=*(R+*(INS+i*2+1)*N+jt);
						for (jx=0;jx<N;jx++)
							FWxt[jx*(N-1)+it]=FWxt[jx*(N-1)+it]+*(HD2+k1*N*2+jx*2+t1);
					}
				}
				*(MEM+Cnt*5)  =i;
				*(MEM+Cnt*5+1)=0;
				*(MEM+Cnt*5+2)=0;
				*(MEM+Cnt*5+3)=0;
				*(MEM+Cnt*5+4)=1;


				COST[Cnt]=sumssea(FWxt);
				
				if (COST[Cnt]<=Maxi)
				{
					Maxi=COST[Cnt];
					indx=Cnt;
				}
				Cnt=Cnt+1;

				for (kt=0;kt<N;kt++)
					for (jt=0;jt<gn-1;jt++)
						FWxt[kt*(N-1)+jt]=FWex[kt*(N-1)+jt];

		}

		//The following loop generates the neighours obtained by applying the permutations to every orbit having size > 1.
		//Then the corresponding cost is computed.
		for (i=0;i<Cx;i++) 
		{
			MO=SZ[Rx[i]];
			for (k=1;k<MO;k++)
			{
				for (it=0;it<gn-1;it++) 
					for (jt=0;jt<MO;jt++)
					{
						t0=SR[it*N+*(R+Rx[i]*N+jt)];
						t1=SR[it*N+*(R+Rx[i]*N+((jt+k)%MO))];
						if (t0 != t1)
						{
							k1=*(R+Rx[i]*N+jt);
							for (jx=0;jx<N;jx++)
								FWxt[jx*(N-1)+it]=FWxt[jx*(N-1)+it]+*(HD2+k1*N*2+jx*2+t1);
						}

					}
					*(MEM+Cnt*5)  =i;
					*(MEM+Cnt*5+1)=k;
					*(MEM+Cnt*5+2)=0;
					*(MEM+Cnt*5+3)=0;
					*(MEM+Cnt*5+4)=2;

					COST[Cnt]=sumssea(FWxt);
				
					if (COST[Cnt]<=Maxi)
					{
						Maxi=COST[Cnt];
						indx=Cnt;
					}
					Cnt=Cnt+1;

					for (kt=0;kt<N;kt++)
						for (jt=0;jt<gn-1;jt++)
							FWxt[kt*(N-1)+jt]=FWex[kt*(N-1)+jt];
			}
		}

		
		for (it=0;it<N;it++) Sc[it]=S[it];

		//The S-box having the lowest cost in the neighbourhood is assigned to Sc.
		if (*(MEM+indx*5+4)==1) 
		{	i=*(MEM+indx*5);	MO=SZS[i];
			for (jt=0;jt<MO;jt++)
			{
				Sc[ *(R+*(INS+i*2)*N+jt)  ]=S[ *(R+*(INS+i*2+1)*N+jt)];
				Sc[ *(R+*(INS+i*2+1)*N+jt)]=S[ *(R+*(INS+i*2)*N+jt)  ];
			}}
		else if (*(MEM+indx*5+4)==2) 
		{	i=*(MEM+indx*5);	k=*(MEM+indx*5+1);	MO=SZ[Rx[i]];
			for (jt=0;jt<MO;jt++)
				Sc[ *(R+Rx[i]*N+jt)]=S[ *(R+Rx[i]*N+((jt+k)%MO)) ];}

		//It is checked whether Sc is generated previously.
		S2T(Sc,ts);
		for (i=0;i<K;i++)
		{
			CHK=0;
			for (j=0;j<N/2;j++)
				if (*(Q3+i*N/2+j)==ts[j])
					CHK=CHK+1;
				else
					break;
			if (CHK==N/2)
				break;
		}

		//If Sc is one of the previous iteration outputs, the following loop runs and a new S-box having the smallest possible cost is chosen. 
		while (CHK==N/2)
		{
			Maxi=1.69e+100;

			for (i=0;i<N;i++)	Sc[i]=S[i];

			COST[indx]=1.7e+100;
			for (i=0;i<Cnt;i++)
				if (COST[i]<=Maxi)
				{
					Maxi=COST[i];
					indx=i;
				}	

		if (*(MEM+indx*5+4)==1)
		{	i=*(MEM+indx*5);	MO=SZS[i];
			for (jt=0;jt<MO;jt++)
			{
				Sc[ *(R+*(INS+i*2)*N+jt)  ]=S[ *(R+*(INS+i*2+1)*N+jt)];
				Sc[ *(R+*(INS+i*2+1)*N+jt)]=S[ *(R+*(INS+i*2)*N+jt)  ];
			}}
		else if (*(MEM+indx*5+4)==2) 
		{	i=*(MEM+indx*5);	k=*(MEM+indx*5+1);	MO=SZ[Rx[i]];
			for (jt=0;jt<MO;jt++)
				Sc[ *(R+Rx[i]*N+jt)]=S[ *(R+Rx[i]*N+((jt+k)%MO)) ];}

			S2T(Sc,ts);
			for (i=0;i<K;i++)
			{
				CHK=0;
				for (j=0;j<N/2;j++)
					if (*(Q3+i*N/2+j)==ts[j])
						CHK=CHK+1;
					else
						break;
				if (CHK==N/2)
					break;
			}
		}

		for (i=0;i<N;i++)	S[i]=Sc[i];
		
		//The iteration output is stored in Q3.
		S2T(S,ts);
		for (i=0;i<N/2;i++)
		{
			*(Q3+cnt)=ts[i];
			cnt=cnt+1;
		}

		for (i=0;i<N;i++)
			for (j=0;j<n;j++)
				Sb[i*n+j]=DB[S[i]*n+j];

		for (i=0;i<gn-1;i++)
			for (j=0;j<N;j++)
				SR[i*N+j]=0;

		for (i=0;i<gn-1;i++)
		{
			for (j=0;j<N;j++)
				for (k=0;k<n;k++)
				{
					SR[i*N+j]=SR[i*N+j]^(Sb[j*n+k]*CF[i*n+k]);
					T[j]=SR[i*N+j];
				}
			walsh(T,FW);
			for (j=0;j<N;j++)
				FWex[j*(N-1)+i]=FW[j];
			ACOR(FW,AC);
			NLs[i]=N/2-findmaxwh(FW)/2;
			ACs[i]=findmaxac(AC);
			ADs[i]=anf(T);
		}

		for (i=0;i<N;i++)
			for (j=0;j<gn-1;j++)
				FWxt[i*(N-1)+j]=FWex[i*(N-1)+j];

		for (i=0;i<N;i++)
			for (j=0;j<N;j++)
				DUA[j*N+i]=0;

		//rotSch: Checks whether S is symmetric S-box.
		if (rotSch(S)==1) 
		{
			printf("error.Srot");
			return 0;
		}

		for (i=0;i<N;i++)
			SS[i]=S[i];

		//In the following, it is checked whether S is bijective.
		qsort(SS,N,sizeof(int),muk);
		for (i=0;i<N;i++)
		{
			if (SS[i]!=i)
			{
				printf("\nerror.sort");
				return 0;
			}
		}

		for (i=1;i<N;i++)
			for (j=0;j<N;j++)
			{
				tmp=S[j^i]^S[j];
				DUA[tmp*N+i]=DUA[tmp*N+i]+1;
			}

		wxor=findmaxor(DUA);
		NLS=findmingl(NLs);
		ACS=findmaxgl(ACs);
		ADS=findmaxgl(ADs);

		if ((K/100)*100==K)
			printf("\n\n KUP=%d K=%d NL=%d AC=%d d=%d X=%d r=%d C=%f (%d, %d)",KUP,K,NLS,ACS,ADS,wxor,r,sumssea(FWxt),Cnt,*(MEM+indx*5+4));

		//Best results are saved.
		if (NLS>=NLA || wxor<=wxa || ACS<=ACA)
		{
			r=r+1;
			if (NLS>NLA)	NLA=NLS;
			if (wxor<wxa)	wxa=wxor;
			if (ACS<ACA)	ACA=ACS;
			fprintf(outr, "\n\n (New) KUP=%d K=%d NL=%d AC=%d d=%d X=%d r=%d COST=%f (%d, %d) ; NLA =%d ACA=%d XA=%d",KUP,K,NLS,ACS,ADS,wxor,r,sumssea(FWxt),Cnt,*(MEM+indx*5+4),NLA,ACA,wxa);
			printf("\n\n (New) KUP=%d K=%d NL=%d AC=%d d=%d X=%d r=%d COST=%f (%d, %d) ; NLA =%d ACA=%d XA=%d",KUP,K,NLS,ACS,ADS,wxor,r,sumssea(FWxt),Cnt,*(MEM+indx*5+4),NLA,ACA,wxa);
			fprintf(outr,"\nSbox=[");
			for (i=0;i<N;i++)
				fprintf(outr,"%d,",S[i]);
			fprintf(outr,"];");
			fclose(outr);
			outr=fopen("S10_P.txt", "a");
		}

	}

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	printf ("\n(KUP=%d) Current local time and date: %s", KUP,asctime(timeinfo));
	fprintf (outr,"\n(KUP=%d) Current local time and date: %s", KUP,asctime(timeinfo));
	fclose(outr);
	outr=fopen("S10_P.txt", "a");

	}

	fclose(outr);
	return 0;
}

void S2T(int *S, unsigned int *ts)
{
	int i,j,k=N-1;
	for (i=0;i<N/2;i++)
	{
		ts[i]=0;
		for (j=0;j<2;j++)
		{
			ts[i]=ts[i]^(S[k]<<(12*j));
			k=k-1;
		}
	}
}

int rot2d(int d)
{
	int e,i;
	e=0;
	for (i=0;i<n;i++)
		e=e^(((d&(1<<(n-1-p[i])))>>(n-1-p[i]))<<(n-1-i));
	return e;
}

int rotSch(int *S)
{
	int i,j,k=0;

	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
		{
			if (S[*(ROT+i*N+j)]!=*(ROT+S[i]*N+j)) 
			{
				k=1;
				break;
			}
		}
			
	return k;
}


void walsh(signed char *T, int *FW)
{
	int i,j;
	for (i=0;i<N;i++)
	{
		FW[i]=0;
		for (j=0;j<N;j++)
			FW[i]=FW[i]+*(HD+j*N*2+i*2+T[j]);
	}
}

void ACOR(int *FW, int *AC) 
{ 
	int i, j, k, a; 
	for (i=0; i<N; i++)
		AC[i]=FW[i]*FW[i]; 
	for (i=1; i<N; i=i<<1) 
		for (j=0; j<N; j=j+(i<<1)) 
			for (k=j; k<j+i; k++)
			{
				a=AC[k]-AC[k+i]; 
				AC[k]=AC[k]+AC[k+i]; 
				AC[k+i]=a; 
			} 
			for (i=0; i<N; i++)
				AC[i]>>=n; 
} 

int anf(signed char *T)
{
	int i,j,AF[N],Maxi=0;
	for (i=0;i<N;i++)
	{
		AF[i]=0;
		for (j=0;j<N;j++)
			AF[i]=AF[i]+*(A+i*N+j)*T[j];
		AF[i]=(AF[i]%2)*BS[i];
		if (AF[i]>Maxi)
			Maxi=AF[i];
	}
	return Maxi;
}

int muk(const void *x,const void *y)
{
	int f;
	f=((*(int *)x)-(*(int *)y));
	return f;
}

void dec2bin(int d, int *b)
{
	int i;
	for (i=0;i<n;i++)
		b[i]=(d&(1<<i))>>i;

}


long double sumsse(int *FW)
{
	int i;
	long double sum=0;
	for (i=0;i<N;i++)
		if (FW[i]<0)
			sum=sum+FW22[-FW[i]];
		else
			sum=sum+FW22[FW[i]];
	return sum;
}


long double sumssea(int *FWxt)
{
	int i,j;
	long double sum=0,sumx;
	
	for (j=0;j<gn-1;j++)
	{
		sumx=0;
		for (i=0;i<N;i++)
			if (*(FWxt+i*(N-1)+j)<0)
				sumx=sumx+FW22[-(*(FWxt+i*(N-1)+j))];
			else
				sumx=sumx+FW22[  *(FWxt+i*(N-1)+j) ];
		sum=sum+sumx;
	}

	return sum;
}

int findmaxac(int *tt)
{
	int i,D,Maxi=-1;
	for (i=1;i<N;i++)
	{
		D=tt[i];
		if (tt[i]<0)
			D=-tt[i];
		if (D>Maxi)
			Maxi=D;
	}
	return Maxi;
}

int findmaxwh(int *tt)
{
	int i,D,Maxi=-1;
	for (i=0;i<N;i++)
	{
		D=tt[i];
		if (tt[i]<0)
			D=-tt[i];
		if (D>Maxi)
			Maxi=D;
	}
	return Maxi;
}

int findmaxor(int *DUA)
{
	int i,j,D,Maxi=-1;
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
		{
			D=*(DUA+i*N+j);
			if (D>Maxi)
				Maxi=D;
		}
	return Maxi;
}

int findmaxgl(int *tt)
{
	int i,D,Maxi=-1;
	for (i=0;i<gn-1;i++)
	{
		D=tt[i];
		if (tt[i]<0)
			D=-tt[i];
		if (D>Maxi)
			Maxi=D;
	}
	return Maxi;
}

int findmingl(int *tt)
{
	int i,D,Mini=5000;
	for (i=0;i<gn-1;i++)
	{
		D=tt[i];
		if (tt[i]<0)
			D=-tt[i];
		if (D<Mini)
			Mini=D;
	}
	return Mini;
}

void RSBOX(int *S)
{
	int ir,i,j,k,chk;
	long int seed;
	int PI[N],SS[N];

	for (i=0;i<gn-1;i++)
	{
		chk=0;
		while (chk==0)
		{
			chk=1;
			cek=cek+100000;
			seed=(long int) time(0)+cek;
			TRandomMersenne rg(seed);
			ir=rg.BRandom();
			j=ir&(N-1);

			if (j>=gn || j==0 || SZ[j]!=SZ[i+1])	chk=0;
			for (k=0;k<i;k++)
				if (PI[k]==j)
					chk=0;
			if (chk==1)	PI[i]=j;
		}
	}

	for (i=0;i<gn-1;i++)
	{
		k=SZ[i+1];
		while (k>=SZ[i+1])
		{
			cek=cek+100000;
			seed=(long int) time(0)+cek;
			TRandomMersenne rg(seed);
			ir=rg.BRandom();
			k=ir&(0xF);
			if (k<SZ[i+1])
				for (j=0;j<SZ[i+1];j++)
					S[*(R+(i+1)*N+j)]=*(R+PI[i]*N+((j+k)%(SZ[i+1])));
		}

	}

	S[0]=0;

	if (rotSch(S)==1)
	{
		printf("error.Srot");
		for (i=0;i<N;i++)
			printf("\n%d %d",i,S[i]);
	}

	for (i=0;i<N;i++)
		SS[i]=S[i];

	qsort(SS,N,sizeof(int),muk);
	for (i=0;i<N;i++)
	{
		if (SS[i]!=i)
		{
			printf("\nerror.sort");
			break;
		}
	}

}

void orbit(int k)
{
	int i,j,l,chk=1,b[n],a[n];
	h=0;
	for (i=0;i<n;i++)
	{
		*(G+h*n+i)=*(Rb+k*n+i);
		b[i]=*(Rb+k*n+i);
	}
	while (chk!=0)
	{
		h=h+1;
		for (i=0;i<n;i++)
			a[i]=b[p[i]];
		for (i=0;i<n;i++)
			b[i]=a[i];
		for (i=0;i<h;i++)
		{
			l=0;
			for (j=0;j<n;j++)
				if (*(G+i*n+j)==a[j])
					l=l+1;
			if (l==n)
			{
				chk=0;
				break;
			}
		}
		if (chk!=0)
			for (i=0;i<n;i++)
				*(G+h*n+i)=a[i];
	}
}


void repres()
{
	int i,j,k=0;
	for (i=0;i<N;i++)
	{
		for (j=0;j<n;j++)
			*(Bc+i*n+n-1-j)=(i&(1<<j))>>j;
	}
	for (i=0;i<N;i++)
		if (*(Bc+i*n)!=-1)
		{
			for (j=0;j<n;j++)
				*(Rb+k*n+j)=*(Bc+i*n+j);
			orbit(k);
			SZ[k]=h;
			for (j=0;j<h;j++)
				*(Bc+(512*(*(G+j*n)) + 256*(*(G+j*n+1)) + 128*(*(G+j*n+2)) + 64*(*(G+j*n+3)) + 32*(*(G+j*n+4)) + 16*(*(G+j*n+5)) + 8*(*(G+j*n+6)) + 4*(*(G+j*n+7)) + 2*(*(G+j*n+8)) + 1*(*(G+j*n+9)))*n)=-1;
			for (j=0;j<h;j++)
				*(R+k*N+j)=512*(*(G+j*n)) + 256*(*(G+j*n+1)) + 128*(*(G+j*n+2)) + 64*(*(G+j*n+3)) + 32*(*(G+j*n+4)) + 16*(*(G+j*n+5)) + 8*(*(G+j*n+6)) + 4*(*(G+j*n+7)) + 2*(*(G+j*n+8)) + 1*(*(G+j*n+9));
			k=k+1;
			gn=k;
		}
}

long double fact(int k)
{
	long double E=0,LI=0;
	int i;
	for (i=1;i<=k;i++)
	{
		E=E+1;
		LI=LI + log(E)/log(2.0);
	}
	return LI;
}




