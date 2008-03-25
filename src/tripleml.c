#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

double caseIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1);
double caseIIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1);
double p0(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p1(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p2(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p4(double gamma0, double gamma1, double theta0, double theta1, int tree);
double pxxx(double gamma0, double gamma1, double theta0, double theta1);
double pxxy(double gamma0, double gamma1, double theta0, double theta1);
double pyxx(double gamma0, double gamma1, double theta0, double theta1);
double pxyx(double gamma0, double gamma1, double theta0, double theta1);
double pxyz(double gamma0, double gamma1, double theta0, double theta1);
void tripleProb(double *gamma0, double *gamma1, double *theta0, double *theta1, double *p0,double *p1,double *p2,double *p3,double *p4);
double findMaxtrixElement(int row, int col, int nrow, double *a);
double rndu (void);
int nextNucleotide (int nucleotide, double branch);
void simnucleotides (int *nucleotide, int *length, double *branch, int *next);


int main(void)
{
	int i, nucleotide[1000], length=1000, next[1000];
	double branch=1.0;

	for(i=0; i<1000; i++)
	{
		nucleotide[i] = i % 4 + 1;
		printf("n %d\n",nucleotide[i]);
	}

	simnucleotides(nucleotide, &length, &branch, next);

	for(i=0; i<1000; i++)
	{
		printf("n %d %d\n",nucleotide[i], next[i]);
	}
	return(1);
}

/*nucleotide: 0 A 1 G 2 C 3 T*/
void simnucleotides (int *nucleotide, int *length, double *branch, int *next)
{
	int i;

	for(i=0; i< (*length); i++)
	{
		next[i] = nextNucleotide(nucleotide[i], (*branch));
	}

}



int nextNucleotide (int nucleotide, double branch)
{
	double prob[4], ran = rndu();
	int i;


	for(i=0; i<4; i++)
	{
		if(i == (nucleotide-1))
			prob[i] = 0.25*(1+3*exp(-4*branch/3));
	 	else
			prob[i] = 0.25*(1-exp(-4*branch/3));
	}

	if(ran < prob[0])
		return(1);
	else if(ran < prob[0] + prob[1])
		return(2);
	else if(ran < prob[0] + prob[1] + prob[2])
		return(3);
	else 
		return(4);
		
}	

double findMaxtrixElement(int row, int col, int nrow, double *a)
{
	return (a[nrow*(col-1) + row-1]);
}

double rndu (void)
{
	return((rand() / ((double)RAND_MAX + 1)));
}

void tripleProb(double *gamma0, double *gamma1, double *theta0, double *theta1, double *p0,double *p1,double *p2,double *p3,double *p4)
{ 
	
	*p0 = pxxx(*gamma0,*gamma1,*theta0,*theta1);
	*p1 = pxxy(*gamma0,*gamma1,*theta0,*theta1);
	*p2 = pyxx(*gamma0,*gamma1,*theta0,*theta1);
	*p3 = *p2;
	*p4 = pxyz(*gamma0,*gamma1,*theta0,*theta1);

}



double caseIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1)
{
	double p;
    	p = 4*(1-exp((a1-2/theta1)*gamma0))/(theta0*theta1*(-a0+2/theta0)*(-a1+2/theta1));
    	return(p);
}

double caseIIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1)
{
	double p;
    	p = 12.0/(theta0*theta0*(-a0+2/theta0)*(-a1+6/theta0));
    	return(p);
} 




double p0(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = phy + 3*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8.0*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 6*exp(-(8.0*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 1 + 3*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8.0*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 6*exp(-(8.0*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
    	return (p);
}


double p1(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 3*phy + 9*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 6*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 3 + 9*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 6*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
printf("P %f %f %f %f %d %f %f %f\n",gamma0,  gamma1, theta0,  theta1, tree, a0,  a1, p);
    	return (p);

}

double p2(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 3*phy - 3*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 3 - 3*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;

    	return (p);

}



double p4(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 6*phy - 6*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 12*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 12*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 6 - 6*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 12*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 12*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
   
    	return (p);

}

double pxxx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
   	p = p0(gamma0, gamma1, theta0, theta1, 0) + 3 * p0(gamma0, gamma1, theta0, theta1, 1);
   	return (p);
}


double pxxy(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;   
	p = p1(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1); 
   	return(p);
}


double pyxx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;   
	p = p2(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1);  
   	return(p);
}


double pxyx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
   	p = p2(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1); 
   	return(p);
}


double pxyz(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
	p = p4(gamma0, gamma1, theta0, theta1, 0) + 3*p4(gamma0, gamma1, theta0, theta1, 1);
   	return(p);
}








