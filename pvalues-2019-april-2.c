/* P-VALUES FOR CATEGORY AND CLASSIFICATION QSAR */
/* PROBABILITIES FOR EACH CLASS ARE GIVEN */

/* AUTHOR: A. GOLBRAIKH  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MESS(a) printf("\nError opening file %s.\n",a);
#define MAXLENGTH 256

FILE *of;
char ofname[MAXLENGTH];

void compilation_date(void)
{
	char buff[13];
	char m[4];
	int d,y;
	sprintf(buff,"%s\n",__DATE__);
	sscanf(buff,"%s%d%d",m,&d,&y);
	printf("\nCOMPILATION DATE: %s %d, %d\n",m,d,y);
}

void open_file(void)
{
	if((of=fopen(ofname,"w"))==NULL) { MESS(ofname); exit(EXIT_FAILURE); }
}

long number_of_compounds(int dim, long *n)
{
	int i;
	long NTot=0;
	for(i=0;i<=dim;i++) NTot+=n[i];
	return NTot;
}

double Total_Number_of_Combinations(int dim, long n)
{
	return n*log(dim);
}

long find_max(int dim,long *n)
{
	int i;
	long q=n[0];
	for(i=1;i<=dim;i++)
		if(q<n[i]) q=n[i];
	return q;
}

long find_characteristic(double p)
{
	long i=0;
	while(i<p) i++;
	return i;
}

long calculate_maxerrors(int dim, long *maxerror, long *n, int option)
{
	int i;
	long maxtotalerror=0;
	for(i=0;i<=dim;i++)
	{
		int t;
		if(option==2) 
		{
			t = (i>dim-i) ? i : dim-i;
			maxerror[i]=t*n[i];
		}
		else maxerror[i]=n[i];
		maxtotalerror+=maxerror[i];
	}
	return maxtotalerror;
}

double logfactorial(long k)
{
	double q=0;
	long i;
	for(i=1;i<=k;i++) q+=log(i);
	return q;
}

void multinomial(int dim, int dim1, int NClass, int n, double comb,
				 double q, double logfactn, long error, 
				 int option, double *distriberrors, double *combprobs, long div, double *prob)
{
	int i;
	if(dim==0)
	{	
		double b;
//		printf("%d ",n);
		b=logfactorial(n);
		comb+=b;
		q+=b-n*log(prob[dim1]);
		if(option==2) error+=abs(NClass-dim1)*n;
		else if(NClass+dim-dim1) error+=n;
//	      printf("%le %ld\n",exp(logfactn-q-div),error);
		distriberrors[error]+=exp(logfactn-q-div);
		combprobs[error]+=exp(logfactn-comb-div);
	}
	else
	{
		for(i=0;i<=n;i++)
		{
//			printf("%d ",i);
			double b;
			if(i>0)
			{
				b=logfactorial(i-1);	
				comb-=b;
				q-=b-(i-1)*log(prob[dim1-dim]);
				if(option==2) error-=abs(NClass+dim-dim1)*(i-1);
				else if(NClass+dim-dim1) error-=(i-1);
			}
			b=logfactorial(i);
			comb+=b;
			q+=b-i*log(prob[dim1-dim]);
			if(option==2) error+=abs(NClass+dim-dim1)*i;
			else if(NClass+dim-dim1) error+=i;
			multinomial(dim-1,dim1,NClass,n-i,comb,q,logfactn,error,option,distriberrors,combprobs,div,prob);
		}
	}
}

void totalerrordistribution(double *distribtotalerror, double **distriberrors,
							long maxtotalerror, long *maxerror, double totaldistrib, 
							long totalerror,int dim, int dim1, long j, long n)
{
	long i;
//	long j;
//	for(j=0; j<=maxtotalerror;j++)
//	for(j=0; j<=1;j++)
	{
		if(dim==0)
		{
			if(j<=maxerror[dim]) totaldistrib*=distriberrors[dim][j];
			else totaldistrib=0;
//			printf("j=%ld dim=%d distriberrors[%d][%ld]=%le totaldistrib=%le\n",
//				j,dim,dim,j,distriberrors[dim][j],totaldistrib);
			distribtotalerror[n]+=totaldistrib;
//			printf("\tdistribtotalerror[%ld]=%le\n",n,distribtotalerror[n]);
		}
		else
		{
			for(i=0;i<=j;i++)
			{
				if(dim==dim1 && i==0) totaldistrib=1;
				if(i<=maxerror[dim])
				{
					if(i>0) totaldistrib/=distriberrors[dim][i-1]; 
					totaldistrib*=distriberrors[dim][i];
//					printf("i=%ld j=%ld dim=%d distriberrors[%d][%ld]=%le totaldistrib=%le\n",
//						i,j,dim,dim,i,distriberrors[dim][i],totaldistrib);
					totalerrordistribution(distribtotalerror,distriberrors,
							maxtotalerror,maxerror,totaldistrib,totalerror,dim-1,dim1,j-i,n);
				}
				else
				{
					totaldistrib=0;
//					totalerrordistribution(distribtotalerror,distriberrors,
//						maxtotalerror,maxerror,totaldistrib,totalerror,dim-1,dim1,j,n); 
					break;

				}
			}
		}
	}
}

void calculate_normdistriberrors(long p,double *logs, double *norms)
{
	long i;
	double s=0;
	for(i=0;i<=p;i++) s+=logs[i];
	for(i=0;i<=p;i++) norms[i]=logs[i]/s;
}

void calculate_cumulative(long p,double *norms, double *cums)
{
	long i;
	cums[0]=norms[0];
	for(i=0;i<p;i++) cums[i+1]=cums[i]+norms[i+1];
}

void main(int argc, char *argv[])
{
	long *n; /* numbers of compounds in each class */
	long maxn; /* maximum of n */
	long NTot; /* total number of compounds */
	long *maxerror; /* max error for each class */
	long error; /* current error for multinomial */
	int option; /* 1 - classification, 2 - category */
	double **logdistriberrors; /* error distribution for classes (no of cases) */
	double **normdistriberrors; /* error distribution for classes (density) */
	double **cumdistriberrors; /* cumulative error distribution (density) */
	double **cumlogdistriberrors; /* cumulative error distribution (no of cases) */
	double **combprobs; /* probablity for error distribution (no of cases) */
	int dim; /* number of classes */
	double logfactn; /* log factorial for classes */
	double logfactNTot; /* log factorial total (max of classes) */
	long div; /* Characteristic of logfactNTot */
	double comb; /* number of combinations for split */
	double q; /* probability of split */
	long maxtotalerror=0; /* maximum total error */
	double *logdistribtotalerror; /* total error distribution (no of cases) */
	double *cumlogdistribtotalerror; /* cumulative error distribution (no of cases) */
	double *normdistribtotalerror; /* total error distribution (density) */
	double *cumdistribtotalerror; /* cumulative error distribution (density) */
	double **prob; /* probability array */
	int i;
	long j;

	printf("CALCULATION OF P-values for classification and category QSAR\n");
	compilation_date();
	if(argc<10 || argc!=atoi(argv[3])*(atoi(argv[3])+1)+4)
	{
		printf("\nUSAGE: pvalues9 outputfile option number_of_classes n1 n2 ... nn p11 p12 ... p21 p22 ... ... pnn\n\n");
		printf("\toutputfile is the name of the output file\n");
		printf("\toption is 1 for classification and 2 for category\n");
		printf("\tn1,n2,...,nn are the counts of compounds in each class\n\n");
		printf("\tpij are the probablilities of compounds of class i\n");
		printf("\tdistributed between n classes j=1,2,...,n\n\n");
		printf("\tIt should be: pi1+pi2+...+pin=1\n\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		strcpy(ofname,argv[1]);
		option=atoi(argv[2]);
		open_file();

			if(option==1) printf("\nCLASSIFICATION QSAR\n");
			else if(option==2) printf("\nCATEGORY QSAR\n");
			else 
			{
				printf("\nWRONG OPTION\n");
				exit(EXIT_FAILURE);
			}

			if(option==1) fprintf(of,"\nCLASSIFICATION QSAR\n");
			else if(option==2) fprintf(of,"\nCATEGORY QSAR\n");
			else 
			{
				fprintf(of,"\nWRONG OPTION: %d\n",option);
				exit(EXIT_FAILURE);
			}

		dim=atoi(argv[3])-1;
 
             if(option==1)
             {
			printf("\nNumber of Classes: %d\n",dim+1);
			fprintf(of,"\nNumber of Classes: %d\n",dim+1);
             }

  		  if(option==2)
             {
			printf("\nNumber of Categories: %d\n",dim+1);
			fprintf(of,"\nNumber of Categories: %d\n",dim+1);
             }

		n=malloc((dim+1)*sizeof(long));
		prob=(double**) malloc((dim+1)*sizeof(double *));
		for(i=0;i<=dim;i++) prob[i]=malloc((dim+1)*sizeof(double));
		for(i=0;i<=dim;i++) n[i]=atoi(argv[i+4]);
		for(i=0;i<=dim;i++)
		for(j=0;j<=dim;j++) prob[i][j]=atof(argv[dim+5+(dim+1)*i+j]);
	}

		printf("Probabilities:\n\n");
		for(i=0;i<=dim;i++)
		{
			printf("\tClass %d: ",i);
			for(j=0;j<=dim;j++) printf(" %le",prob[i][j]);
			printf("\n");
		}

		fprintf(of,"Probabilities:\n\n");
		for(i=0;i<=dim;i++)
		{
			fprintf(of,"Class %d:",i);
			for(j=0;j<=dim;j++) fprintf(of," %le",prob[i][j]);
			fprintf(of,"\n");
		}

	NTot=number_of_compounds(dim,n);

        if(option==1)
        { 
		printf("\nTotal Number of Compounds in All Classes: %ld\n",NTot);
		fprintf(of,"\nTotal Number of Compounds in All Classes: %ld\n",NTot);
        }

        if(option==2)
        { 
		printf("\nTotal Number of Compounds in All Categories: %ld\n",NTot);
		fprintf(of,"\nTotal Number of Compounds in All Categories: %ld\n",NTot);
        }

	maxerror=malloc((dim+1)*sizeof(long));
	maxtotalerror=calculate_maxerrors(dim,maxerror,n,option);

		for(i=0;i<=dim;i++)
		{ 
			if(option==1) 
                { printf("\nNumber of Compounds in Class %d: %ld; ",i+1,n[i]);
            	  printf("Maximum error: %ld",maxerror[i]);
                }
                if(option==2)
                { 
			  fprintf(of,"\nNumber of Compounds in Class %d: %ld; ",i+1,n[i]);
            	  fprintf(of,"Maximum error: %ld",maxerror[i]);
                }
		}
		printf("\n\nMaximum Total Error: %ld\n\n",maxtotalerror);
		fprintf(of,"\n\nMaximum Total Error: %ld\n\n",maxtotalerror);
    
	logdistriberrors=(double **)malloc(sizeof(double*)*(dim+1));
	for(i=0;i<=dim;i++) logdistriberrors[i]=(double *) calloc(maxerror[i]+1,sizeof(double));
	cumlogdistriberrors=(double **)malloc(sizeof(double*)*(dim+1));
	for(i=0;i<=dim;i++) cumlogdistriberrors[i]=(double *) calloc(maxerror[i]+1,sizeof(double));
	combprobs=(double **)malloc(sizeof(double*)*(dim+1));
	for(i=0;i<=dim;i++) combprobs[i]=(double *) calloc(maxerror[i]+1,sizeof(double));
	normdistriberrors=(double **)malloc(sizeof(double*)*(dim+1));
	for(i=0;i<=dim;i++) normdistriberrors[i]=(double *) calloc(maxerror[i]+1,sizeof(double));
	cumdistriberrors=(double **)malloc(sizeof(double*)*(dim+1));
	for(i=0;i<=dim;i++) cumdistriberrors[i]=(double *) calloc(maxerror[i]+1,sizeof(double));
	logdistribtotalerror=calloc(maxtotalerror+1,sizeof(double));
	normdistribtotalerror=calloc(maxtotalerror+1,sizeof(double));
	cumdistribtotalerror=calloc(maxtotalerror+1,sizeof(double));
	cumlogdistribtotalerror=calloc(maxtotalerror+1,sizeof(double));
//	logfactNTot=Total_Number_of_Combinations(dim+1,NTot);
//	div=find_characteristic(logfactNTot);
//	printf("div=%ld\n",div);
//	maxn=find_max(dim,n);
//	logfactNTot=Total_Number_of_Combinations(dim+1,maxn);
//	div=find_characteristic(logfactNTot);
	div=0;
	for(i=0;i<=dim;i++)
	{
		error=0; 
		q=0;
		comb=0;
		logfactn=logfactorial(n[i]);
//		printf("%d %le\n",n[i],exp(logfactn));
		multinomial(dim,dim,i,n[i],comb,q,logfactn,error,option,logdistriberrors[i],combprobs[i],div,prob[i]);
		calculate_normdistriberrors(maxerror[i],logdistriberrors[i],normdistriberrors[i]);
		calculate_cumulative(maxerror[i],normdistriberrors[i],cumdistriberrors[i]);
		calculate_cumulative(maxerror[i],logdistriberrors[i],cumlogdistriberrors[i]);
	}

		if(option==1)
           {
		  printf("\n     Distribution of errors for classes\n");
		  fprintf(of,"\n     Distribution of errors for classes\n");
		}
		if(option==2)
           {
		  printf("\n     Distribution of errors for categories\n");
		  fprintf(of,"\n     Distribution of errors for categories\n");
		}
			for(i=0;i<=dim;i++)
			{
				if(option==1) printf("\n       Class %d\n\n",i+1);
				if(option==2) printf("\n	    Category %d\n\n",i+1);
				printf("--------------------------------------------------------------------------------------------\n");
				printf("  Error     Combinations      Probability      Cumulative        Density        Cumulative\n");
				printf("--------------------------------------------------------------------------------------------\n");
				for(j=0;j<=maxerror[i];j++) 
					printf("  %5ld      %8.3le        %8.3le        %8.3le      %8.3le       %8.3le\n",
					       j,combprobs[i][j],logdistriberrors[i][j],cumlogdistriberrors[i][j],
						     normdistriberrors[i][j],cumdistriberrors[i][j]);
			}

			for(i=0;i<=dim;i++)
			{
				fprintf(of,"\n       Class %d\n\n",i+1);
				fprintf(of,"--------------------------------------------------------------------------------------------\n");
				fprintf(of,"  Error     Combinations      Probability      Cumulative        Density        Cumulative\n");
				fprintf(of,"--------------------------------------------------------------------------------------------\n");
				for(j=0;j<=maxerror[i];j++) 
					fprintf(of,"  %5ld      %8.3le        %8.3le        %8.3le      %8.3le       %8.3le\n",
					       j,combprobs[i][j],logdistriberrors[i][j],cumlogdistriberrors[i][j],
						     normdistriberrors[i][j],cumdistriberrors[i][j]);
			}

//		printf("\nmaxtotalerror=%ld\n\n",maxtotalerror);
	
	for(j=0;j<=maxtotalerror;j++)
	{
		totalerrordistribution(logdistribtotalerror,combprobs, //logdistriberrors,
			maxtotalerror,maxerror,1,0,dim,dim,j,j);
	}
	calculate_normdistriberrors(maxtotalerror,logdistribtotalerror,normdistribtotalerror);
	calculate_cumulative(maxtotalerror,normdistribtotalerror,cumdistribtotalerror);
	calculate_cumulative(maxtotalerror,logdistribtotalerror,cumlogdistribtotalerror);
		
		printf("\n     Total distribution of errors\n\n");
		printf("---------------------------------------------------------------------------\n");
		printf("  Error     Combinations      Cumulative       Density        Cumulative\n");
		printf("---------------------------------------------------------------------------\n");
		for(j=0;j<=maxtotalerror;j++) 
			printf("  %5ld      %8.3le       %8.3le      %8.3le       %8.3le\n",j,
			       logdistribtotalerror[j],cumlogdistribtotalerror[j],
				   normdistribtotalerror[j],cumdistribtotalerror[j]);

		fprintf(of,"\n     Total distribution of errors\n\n");
		fprintf(of,"---------------------------------------------------------------------------\n");
		fprintf(of,"  Error     Combinations      Cumulative       Density        Cumulative\n");
		fprintf(of,"---------------------------------------------------------------------------\n");
		for(j=0;j<=maxtotalerror;j++) 
			fprintf(of,"  %5ld      %8.3le       %8.3le      %8.3le       %8.3le\n",j,
			       logdistribtotalerror[j],cumlogdistribtotalerror[j],
				   normdistribtotalerror[j],cumdistribtotalerror[j]);
		
	for(i=0;i<=dim;i++) free(logdistriberrors[i]);
	free(logdistriberrors); 
	for(i=0;i<=dim;i++) free(cumlogdistriberrors[i]);
	free(cumlogdistriberrors); 
	for(i=0;i<=dim;i++) free(combprobs[i]);
	free(combprobs); 
	for(i=0;i<=dim;i++) free(normdistriberrors[i]);
	free(normdistriberrors); 
	for(i=0;i<=dim;i++) free(cumdistriberrors[i]);
	free(cumdistriberrors); 
	free(maxerror);
	free(logdistribtotalerror);
	free(normdistribtotalerror);
	free(cumdistribtotalerror);
	free(cumlogdistribtotalerror);
	for(i=0;i<=dim;i++) free(prob[i]);
	free(prob);
	free(n);
}
