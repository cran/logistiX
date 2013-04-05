#include <R.h>
#include <stdio.h>
//#include <stdafx.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <malloc.h>
#define PREC 0.00000001
#define LPREC 18.42
#define MLPREC 0.0001

/*********************************************************************
 *                                                                   *
 *                                                                   *
 *                 Exact logistic regression                         *
 *                 SE for simulation:                                *
 *                                                                   *
 *                 supports handling of multiple input data sets     *
 *                 of the same size and structure                    *
 *                                                                   *
 *                 also supports STRATA!!!                           *
 *                                                                   *
 *                                                                   *
 *                                                                   *
 *                 Last change 07/26/2000                            *
 *                 includes test output 05/08/2001                   *
 *																	 *
 *																	 *
 *																	 *
 *                 Adaption for R 2012-09, last change 11/05/2012    *
 *                 Georg Heinze & Tobias Ladner                      *
 *                                                                   *
 *                                                                   *
 *********************************************************************/

struct t {
		double *value;
		long double counter;
  long double score;
  long double probability;
  long double LR;
		struct t *next;
		 };

int listlength(struct t* head)
{
	int ll=0;
	struct t* point;
	point=head->next;
	while(point != point->next)
	{
		ll++;
		point=point->next;
	}
	return(ll);
}

int minimum(int a, int b)
{
 if (a<b)
  return(a);
 else
  return(b);
 }

int maximum(int a, int b)
 {
 if (a>b)
  return(a);
 else
  return(b);
 }

long double coco(int a, int b)
 {
 long double combi;
 long double coeff=1;
 if(b==0)
  return(1.0);
 if(b==a)
  return(1.0);
 for(combi=(long double)a; combi>(long double)maximum(b, a-b); combi--)
  coeff*=combi;
 for(combi=(long double)minimum(b, a-b); combi>=1.0; combi--)
  coeff=coeff/combi;
 return(coeff);
 }

double absolute(double zahl)
 {
 if(zahl>0.0)
  return(zahl);
 else
  return(-zahl);
 }

double Data(double *array, int obs, int var, int num_var)
 {
 if(var==0)
  return(1.0);
 else
  return(array[obs*num_var+var-1]);
 }

double y(double *array, int obs, int num_var)
 {
 return(array[obs*num_var+num_var-1]);
 }

int read_comment(FILE *com)
 {
 char c;
 int i=0;
 fscanf(com, " %c", &c);
 do
  {
  fscanf(com, "%c", &c);
  i++;
  }
  while(!(i>90) && !(c=='*'));
 if(i>90)
  return(0);
 else
  return(1);
 }

void LXcalc_prob(struct t *head, double theta, int var)
 {
 /*
   calc_prob is designed so that an overflow is avoided
 */

 long double sum, summand, d1, d2, d3;
 struct t* out;
 struct t* in;

 out=head->next;
 while(out->next !=out)
  {
  sum=0.0;
  in=head->next;
  while(in->next != in)
   {
   d1=log(in->counter);
   d2=log(out->counter);
   d3=theta*(in->value[var]-out->value[var]);
   if((d1-d2+d3) > LPREC)
    summand=1.0/PREC;
   else
    summand=exp(d1-d2+d3);
   sum+=summand;
   in=in->next;
   }
  out->probability=1.0/sum;
  out=out->next;
  }
 }

void LXcalc_scores(struct t* head, double beta, int var, double * t_obs, double *score_obs, double *z)
{
 struct t* point;
 double sum_ttp=0;
 double mean=0;
 double variance=0;

 LXcalc_prob(head, beta, var);
 
 point=head->next;
 while(point->next != point)
 {
  mean+=point->value[var]*point->probability;
  sum_ttp+=point->value[var]*point->value[var]*point->probability;
  point=point->next;
 }
 variance=sum_ttp-mean*mean;
 
 point=head->next;
 while(point->next != point)
 {
  point->score=(point->value[var]-mean)*(point->value[var]-mean)/variance;
  if(point->value[var]==t_obs[var])
  {
	  score_obs[0]=point->score;
	  z[0]=(point->value[var]-mean)/sqrt(variance);
  }
  point=point->next;
 }
}

double likelihood(double *t_obs, int var, double beta, struct t* currhead)
 {
/* double sum=0.0; */
 double likelihood=0.0;
 struct t *point;
 LXcalc_prob(currhead, beta, var);

 point=currhead->next;
 while(point->next != point)
  {
  if(point->value[var]==t_obs[var])
   likelihood=point->probability;
  point=point->next;
  }

 return(likelihood);
 }

double LogXact_estimates(double *t_obs, int var, struct t* currhead, 
                         int method)
 {

 /*

 Calculates Median Unbiased (method=1)
        and Maximum Likelihood (method=2)
 estimates of Variable var given all others.
 method=3: If ML-estimates do not exist, a MU-estimate will be calculated.

 */

 
 double Min, Max, lower, upper, currbeta, beta,
        lowbound, upbound, like;
 struct t *point;
 Min=currhead->next->value[var];
 Max=Min;

 point=currhead->next;

 while(point->next != point)
  {
  if(point->value[var] > Max)
   Max=point->value[var];
  if(point->value[var] < Min)
   Min=point->value[var];
  point=point->next;
  }

 if(Min==Max)
  return(9999.0);

 if(method==2 && Max==t_obs[var])
  return(999.0);
 if(method==2 && Min==t_obs[var])
  return(-999.0);
 if(method==3 && (Max==t_obs[var] || Min==t_obs[var]))
  method=1;
 if(method==3)
  method=2;

 currbeta=0.0;
 lowbound=-10.0;
 upbound=10.0;
 beta=upbound+1.0;

 if(method==1)         /* Median Unbiased Estimator */
  do
   {
   lower=0.0;
   upper=0.0;
   LXcalc_prob(currhead, currbeta, var);
   point=currhead->next;
   while(point->next != point)
    {
    if(point->value[var]<=t_obs[var])
     lower+=point->probability;
    if(point->value[var]>=t_obs[var])
     upper+=point->probability;
    point=point->next;
    }

   if(t_obs[var]==Max)
	 lower=0.5;
	if(t_obs[var]==Min)
	 upper=0.5;


	if(upper-lower>PREC)
	 {
	 upbound=currbeta;
	 currbeta=(lowbound+currbeta)/2;
	 }
	else
	 if(lower-upper>PREC)
	  {
	  lowbound=currbeta;
	  currbeta=(upbound+currbeta)/2;
	  }
	 else
	  beta=currbeta;
	} while(beta !=currbeta);
 else if(method==2)        // ML-estimator
  do
	{
	like=likelihood(t_obs,var,currbeta,currhead);
	if(likelihood(t_obs,var,currbeta-MLPREC,currhead)>like)
	 {
	 upbound=currbeta;
	 currbeta=(lowbound+currbeta)/2;
	 }
	else
	if(likelihood(t_obs,var,currbeta+MLPREC,currhead)>like)
	  {
	  lowbound=currbeta;
	  currbeta=(upbound+currbeta)/2;
	  }
	else
	  beta=currbeta;
   } while(beta !=currbeta);

 LXcalc_prob(currhead,0.0,var);
 return(beta);
 }

void LXcalc_LR(double *t_obs, int var, struct t* currhead, int IP, double beta0)
{
 /* calculates the LR statistic for each t(var) */

 struct t* point;
 double *t_loop;
 int j;
 double Min, Max, lmax, l0, beta;

 if (NULL == (t_loop=(double *)R_alloc(IP, sizeof(double))))
 {
	 error("no memory available\n");
 }

 for(j=0; j<IP; j++)
  t_loop[j]=t_obs[j];

 point=currhead->next;
 Min=point->value[var];
 Max=point->value[var];
 while(point->next != point)
 {
  if(point->value[var] < Min)
   Min=point->value[var];
  if(point->value[var] > Max)
   Max=point->value[var];
  point=point->next;
 }



 point=currhead->next;
 while(point->next != point)
 {
  t_loop[var]=point->value[var];
  if(t_loop[var]==Min || t_loop[var]==Max)
  {
   lmax=1;
  }
  else
  {
   beta=LogXact_estimates(t_loop, var, currhead, 2);
   lmax=likelihood(t_loop, var, beta, currhead);
  }
  l0=likelihood(t_loop, var, beta0, currhead);
  point->LR=2*(log(lmax)-log(l0));
  // Rprintf("LR t=%4.0f: %4.8f\n",t_loop[var],point->LR);
  point=point->next;
  t_loop[var]=t_obs[var];
 }
 // getchar();
}

double LXcalc_pLR(double beta, double * t_obs, int var, struct t* currhead, int pmid)
{
 /* calculates p_value of greater/equal (greater/half equal) LR statistic at H0: beta=beta */
 
 double p_value=0, LRobs=0;
 long double sum=0;
 struct t* point;


 point=currhead->next;
 while(point->next != point)
 {
  sum+=point->counter;
  if(point->value[var]==t_obs[var])
   LRobs=point->LR;
  point=point->next;
 }


 LXcalc_prob(currhead, beta, var);

 //Rprintf("Entered LXcalc_pLR. beta=%4.4f, var=%d\n",beta,var);
 
 point=currhead->next;
 while(point->next != point)
 {
  //Rprintf("t=%4.0f probability=%4.4f LR=%4.4f\n",point->value[var],
  // point->probability, point->LR);
  if(point->LR > LRobs)
   p_value+=point->probability;
  if(point->LR == LRobs)
   p_value+=point->probability-(long double)pmid*point->probability/2;
  point=point->next;
 }

 return(p_value);
}

void LXconf_limits(double *t_obs, double *limits, double *p_value,
                   int var, double signif, struct t* currhead, int pmid)
{
 /*

 Calculates Confidence Limits of Estimators
 Stores them into limits[0] and limits[1]
 -999 denotes -INF
 +999 denotes +INF

  pmid: 0 = exact conf-limits
        1 = p-mid exact conf-limits

 */

 double Min, Max, currtheta, theta;
 long double  lowbound, upbound, lower, upper;
 struct t *point;
 double  psmalltail;
 int j;

 for(j=0; j<3; j++)
	 p_value[j]=0;
 LXcalc_prob(currhead, 0.0, var);
 point=currhead->next;
 while(point->next != point)
  {
  if(point->value[var]>t_obs[var])
	  p_value[2]+=point->probability;
  if(point->value[var]<t_obs[var])
	  p_value[1]+=point->probability;
  if(point->value[var]==t_obs[var])
  {
	  p_value[1]+=point->probability*(1.0-0.5*(double)pmid);
	  p_value[2]+=point->probability*(1.0-0.5*(double)pmid);
  }
  point=point->next;
  }
  psmalltail=p_value[1];
  if(p_value[2]<p_value[1])
	  psmalltail=p_value[2];
  p_value[0]=2.0*psmalltail;

 


 int ll=listlength(currhead);
 if(ll==1)
 {
	 for(int i=0; i<2; i++)
		 limits[i]=9999;
	 return;
 }
 else
 {

	Min=currhead->next->value[var];
	Max=Min;

	point=currhead->next;

	while(point->next != point)
	{
		if(point->value[var] > Max)
		Max=point->value[var];
		if(point->value[var] < Min)
			Min=point->value[var];
		point=point->next;
	}

	lowbound=-20.0;
	upbound=20.0;
	theta=upbound+1.0;
	currtheta=0.0;

	if(t_obs[var]==Max)
		limits[1]=999.0;
	else
	{
		do
		{
			lower=0.0;
			upper=0.0;
			LXcalc_prob(currhead, currtheta, var);
			point=currhead->next;
			while(point->next != point)
			{
				if(point->value[var]<t_obs[var])
					lower+=point->probability;
				if(point->value[var]==t_obs[var])
					lower+=point->probability-(long double)pmid*point->probability/2;
				point=point->next;
			}


			if((signif/2-lower)>PREC)
			{
				upbound=currtheta;
				currtheta=(lowbound+currtheta)/2;
			}
			else
				if((lower-signif/2)>PREC)
				{
					lowbound=currtheta;
					currtheta=(upbound+currtheta)/2;
				}
				else
					theta=currtheta;
		} while(theta !=currtheta);
		limits[1]=theta;
	}



	lowbound=-20.0;
	upbound=20.0;
	theta=upbound+1.0;
	currtheta=upbound/2;

	if(t_obs[var]==Min)
		limits[0]=-999.0;
	else
	{
		do	
		{
			upper=0.0;
			LXcalc_prob(currhead, currtheta, var);
			point=currhead->next;
			while(point->next != point)
			{
				if(point->value[var]>t_obs[var])
					upper+=point->probability;
				if(point->value[var]==t_obs[var])
					upper+=point->probability-(long double)pmid*point->probability/2;
    
			    point=point->next;
			}


			if((upper-signif/2)>PREC)
			{
				upbound=currtheta;
				currtheta=(lowbound+currtheta)/2;
			}
			else
				if((signif/2-upper)>PREC)
				{
					lowbound=currtheta;
					currtheta=(upbound+currtheta)/2;
				}
				else
					theta=currtheta;
		} while(theta !=currtheta);
		limits[0]=theta;
	}
	LXcalc_prob(currhead,0.0,var);
 }

}

void LXcl_scores(double *t_obs, double *limits, double *p_value, 
                 int var, double signif, 
                 struct t* currhead, double beta, int pmid, double *chi, double *z)
 {
 /*

 Calculates Confidence Limits of Estimators
 based on the exact distribution of the score statistic
 Stores them into limits[0] and limits[1]
 -999 denotes -INF
 +999 denotes +INF

 */

 double Min, Max, currtheta, theta, score_obs, upper_beta=0, zrun;
 long double  lowbound, upbound, upper;
 struct t *point;
 int it, hs, j;


 Min=currhead->next->value[var];
 Max=Min;
 for(j=0; j<3; j++)
	 p_value[j]=0.0;
 LXcalc_scores(currhead, 0.0, var, t_obs, &score_obs, z);
 chi[0]=score_obs;
 point=currhead->next;
 while(point->next != point)
  {
  if(point->value[var]>t_obs[var])
	  p_value[2]+=point->probability;
  if(point->value[var]<t_obs[var])
	  p_value[1]+=point->probability;
  if(point->value[var]==t_obs[var])
  {
	  p_value[1]+=point->probability*(1.0-0.5*(double)pmid);
	  p_value[2]+=point->probability*(1.0-0.5*(double)pmid);
  }
  if(point->score>score_obs)
	  p_value[0]+=point->probability;
  if(point->score==score_obs)
	  p_value[0]+=point->probability*(1.0-0.5*(double)pmid);
  point=point->next;
  }
 
 LXcalc_scores(currhead, beta, var, t_obs, &score_obs, &zrun);
 point=currhead->next;
 while(point->next != point)
  {
  if(point->score>=score_obs)
   upper_beta+=point->probability;
  point=point->next;
 }
 // Rprintf("\nscore p-value at beta[%d]: %4.4f\n",var,upper_beta);
   
 point=currhead->next;

 while(point->next != point)
  {
  if(point->value[var] > Max)
   Max=point->value[var];
  if(point->value[var] < Min)
   Min=point->value[var];
  point=point->next;
  }

 //Rprintf("Upper bound\n");


 lowbound=beta;
 upbound=beta+20.0;
 theta=upbound+1.0;
 currtheta=beta+5;

 if(t_obs[var]==Max)
  limits[1]=999.0;
 else
  {
  it=0;
  do
   {
   it++;
   hs=0;
   do
   {
    hs++;
    upper=0.0;
    LXcalc_scores(currhead, currtheta, var, t_obs, &score_obs, &zrun);
    point=currhead->next;
    while(point->next != point)
     {
     if(point->score>score_obs)
      upper+=point->probability;
     if(point->score==score_obs)
      upper+=point->probability-(long double)pmid*point->probability/2;
     point=point->next;
     }
    //Rprintf(" currtheta: %4.4f, upper: %4.4f",currtheta,upper);
    //getchar();
    if(upper>upper_beta)
    {
     upbound=currtheta;
     lowbound=beta;
     currtheta=(lowbound+upbound)/2;
    }
   } while(upper>=upper_beta && hs<10);

   
   /*
   if(var==2)
   {
    Rprintf("currtheta: %4.4f upper: %4.4f\n",currtheta,upper);
    getchar();
   }
   */

   if((signif-upper)>PREC)
    {
    upbound=currtheta;
    currtheta=(lowbound+currtheta)/2;
    }
   else
    if((upper-signif)>PREC)
     {
     lowbound=currtheta;
     currtheta=(upbound+currtheta)/2;
	  }
    else
     theta=currtheta;
   } while(theta !=currtheta && it<20);
  limits[1]=theta;
  if(it==20)
   limits[1]=currtheta;
  }

 
 //Rprintf("Lower bound\n");

 lowbound=beta-20.0;
 upbound=beta;
 theta=lowbound-1.0;
 currtheta=beta-5;

 if(t_obs[var]==Min)
  limits[0]=-999.0;
 else
  {
  it=0;
  do
   {
   it++;
   hs=0;
   do
   {
    hs++;
    upper=0.0;
    LXcalc_scores(currhead, currtheta, var, t_obs, &score_obs, &zrun);
    point=currhead->next;
    while(point->next != point)
     {
     if(point->score>score_obs)
      upper+=point->probability;
     if(point->score==score_obs)
      upper+=point->probability-(long double)pmid*point->probability/2;
     point=point->next;
	    }
    //Rprintf(" currtheta: %4.4f upper: %4.4f",currtheta,upper);
    //getchar();
    if(upper>upper_beta)
    {
     lowbound=currtheta;
     upbound=beta;
     currtheta=(upbound+lowbound)/2;
    }
   } while(upper>=upper_beta && hs<10);

   //Rprintf("currtheta: %4.4f upper: %4.4f\n",currtheta,upper);
   //getchar();

   if((upper-signif)>PREC)
    {
    upbound=currtheta;
    currtheta=(lowbound+currtheta)/2;
    }
   else
    if((signif-upper)>PREC)
     {
     lowbound=currtheta;
     currtheta=(upbound+currtheta)/2;
     }
    else
     theta=currtheta;
   } while(theta !=currtheta && it<20);
  limits[0]=theta;
  if(it==20)
   limits[0]=currtheta;
  }
 LXcalc_prob(currhead,0.0,var);

 }

void LXcl_LR(double *t_obs, double *limits,
                 int var, double signif, 
                 struct t* currhead, double beta, int IP, int pmid)
 {
 /*

 Calculates Confidence Limits of Estimators
 based on the exact distribution of the log(LR) statistic
 Stores them into limits[0] and limits[1]
 -999 denotes -INF
 +999 denotes +INF

 */

 double Min, Max, currtheta, theta, upper_beta=0;
 long double  lowbound, upbound,  upper;
 struct t *point;
 int it, hs;


 Min=currhead->next->value[var];
 Max=Min;

   
 point=currhead->next;

 while(point->next != point)
  {
  /*
  if(point->value[var]==t_obs[var])
   LR_obs=point->LR;
   */
  if(point->value[var] > Max)
   Max=point->value[var];
  if(point->value[var] < Min)
   Min=point->value[var];
  point=point->next;
  }

 //Rprintf("Upper bound\n");


 lowbound=beta;
 upbound=beta+20.0;
 theta=upbound+1.0;
 currtheta=beta+5;

 if(t_obs[var]==Max)
  limits[1]=999.0;
 else
  {
  it=0;
  do
   {
   it++;
   LXcalc_LR(t_obs, var, currhead, IP, beta);
   upper_beta=LXcalc_pLR(beta, t_obs, var, currhead, pmid);
   hs=0;
   do
   {
    hs++;
    LXcalc_LR(t_obs, var, currhead, IP, currtheta);
    upper=LXcalc_pLR(currtheta, t_obs, var, currhead, pmid);

    //Rprintf(" currtheta: %4.4f, upper: %4.4f",currtheta,upper);
    //getchar();
    if(upper>upper_beta)
    {
     upbound=currtheta;
     lowbound=beta;
     currtheta=(lowbound+upbound)/2;
    }
   } while(upper>=upper_beta && hs<10);

   
   /*
   if(var==2)
   {
    Rprintf("currtheta: %4.4f upper: %4.4f\n",currtheta,upper);
    getchar();
   }
   */

   if((signif-upper)>PREC)
    {
    upbound=currtheta;
    currtheta=(lowbound+currtheta)/2;
    }
   else
    if((upper-signif)>PREC)
     {
     lowbound=currtheta;
     currtheta=(upbound+currtheta)/2;
	  }
    else
     theta=currtheta;
   } while(theta !=currtheta && it<20);
  limits[1]=theta;
  if(it==20)
   limits[1]=currtheta;
  }

 
 //Rprintf("Lower bound\n");

 lowbound=beta-20.0;
 upbound=beta;
 theta=lowbound-1.0;
 currtheta=beta-5;

 if(t_obs[var]==Min)
  limits[0]=-999.0;
 else
  {
  it=0;
  do
   {
   it++;
   LXcalc_LR(t_obs, var, currhead, IP, beta);
   upper_beta=LXcalc_pLR(beta, t_obs, var, currhead, pmid);
   hs=0;
   do
   {
    hs++;
    LXcalc_LR(t_obs, var, currhead, IP, currtheta);
    upper=LXcalc_pLR(currtheta, t_obs, var, currhead, pmid);
    //Rprintf(" currtheta: %4.4f upper: %4.4f",currtheta,upper);
    //getchar();
    if(upper>upper_beta)
    {
     lowbound=currtheta;
     upbound=beta;
     currtheta=(upbound+lowbound)/2;
    }
   } while(upper>=upper_beta && hs<10);

   //Rprintf("currtheta: %4.4f upper: %4.4f\n",currtheta,upper);
   //getchar();

   if((upper-signif)>PREC)
    {
    upbound=currtheta;
    currtheta=(lowbound+currtheta)/2;
    }
   else
    if((signif-upper)>PREC)
     {
     lowbound=currtheta;
     currtheta=(upbound+currtheta)/2;
     }
    else
     theta=currtheta;
   } while(theta !=currtheta && it<20);
  limits[0]=theta;
  if(it==20)
   limits[0]=currtheta;
  }
 LXcalc_prob(currhead,0.0,var);

 }


double * group_data(double *array, int numvar, int obs, int *groups)
 {
 double * gdata;
 double * groupdata; 
 int i, k, j, equal;
 groups[0]=0;
 equal=0;
 if(NULL == (gdata=(double *)R_alloc(obs*numvar, sizeof(double))))
 {
	 error("no memory available\n");
 }
 for(i=0; i<obs*numvar; i++)
	 gdata[i]=0;
 /*
 for(j=1; j<numvar; j++)
  gdata[j]=Data(array, 0, j, numvar);
 gdata[numvar-1]=1;
   */
 for(i=0; i<obs; i++)     /* Data from CaseData  */
  {
   /*
   Rprintf("CaseData: %2.0f %2.0f  ",Data(array,i,1,numvar),
                              Data(array,i,2,numvar));
   */
  if(groups[0]!=0)
   for(k=0; k<groups[0]; k++) /* Data from GroupData */
   	{
   	equal=1;
   	for(j=1; j<numvar; j++)
   	 equal*=(Data(array, i, j, numvar)==Data(gdata, k, j, numvar));
   	if(equal)
   	 {
     /*
     Rprintf("matches Group: %2.0f %2.0f  w: %2.0f ",
                                    Data(gdata,k,1,numvar),
                                     Data(gdata,k,2,numvar),
                                     y(gdata,k,numvar));
     */
   	 gdata[k*numvar+numvar-1]++;
     /*
     Rprintf(" new w: %2.0f\n",y(gdata,k,numvar));
     */
   	 k=groups[0]+1;
   	 }
   	}
  if(!equal)
  	{
  	
  	for(j=1; j<numvar; j++)
  	 gdata[groups[0]*numvar+j-1]=Data(array, i, j, numvar);
  	gdata[groups[0]*numvar+numvar-1]=1.0;
	groups[0]++;
   /*
   Rprintf(" new group: %2.0f %2.0f  w: %2.0f\n",
          Data(gdata,groups[0],1,numvar),
          Data(gdata,groups[0],2,numvar),
          y(gdata,groups[0],numvar));
          */
  	}
  }
 /*groups[0]++;                groups is effective number of groups */
 /*
 for(k=0; k<groups[0]; k++)
  {
  for(j=1; j<numvar; j++)
   Rprintf("%4.0f ",Data(gdata, k, j, numvar));
  Rprintf("    Weight=%4.0f\n ", y(gdata, k, numvar));
  }
 */

 if (NULL == (groupdata=(double *)R_alloc(numvar*obs, sizeof(double))))
 {
	 error("no memory available");
 }
 
  for(i=0; i<groups[0]*numvar; i++)
   groupdata[i]=gdata[i];
 
 return(groupdata);
 }

double sort_score(double *group, int numvar, int var)
 {
 int j;
 double score=0.0;
 for(j=1; j<numvar; j++)
  if(j != var)
   score+=pow(2.0,(double)j)*group[j-1];
 return(score);
 }

double * sort_groups(double *gdata, int numvar, int groups, int var)
 {
 double * sorted;

 int i, j, k, l;

 if(NULL == (sorted=(double *)R_alloc(groups*numvar, sizeof(double))))
 {
	 error("no memory available");
 }
 for(i=0; i<groups; i++)
  for(j=1; j<=numvar; j++)
   sorted[i*numvar+j-1]=Data(gdata,i,j,numvar);

 for(i=1; i<groups; i++)
  for(k=0; k<i; k++)
   if(sort_score(&gdata[i*numvar],numvar,var) >
      sort_score(&sorted[k*numvar],numvar,var))
    {
    for(l=i-1; l>=k; l--)
     for(j=1; j<=numvar; j++)
      sorted[(l+1)*numvar+j-1]=sorted[l*numvar+j-1];
    for(j=1; j<=numvar; j++)
     sorted[k*numvar+j-1]=gdata[i*numvar+j-1];
    k=i;
    }
 /*
 for(k=0; k<groups; k++)
  {
  for(j=1; j<numvar; j++)
   Rprintf("%4.0f ",Data(sorted, k, j, numvar));
  Rprintf("    W=%4.0f  Sc=%5.0f\n ", y(sorted, k, numvar),
         sort_score(&sorted[k*numvar],numvar,var));
  }
 */

 return(sorted);
 }

int feasible(double * t_obs, double * t_obs0_strat, int * stratmap, double *val, int var, int numvar, int noint, int currobs, int resp, double * feasmap)
 {
 int feas;
 int j;
 feas=1;

 if(noint==0)
 {
	 if(val[0]+feasmap[2*currobs*numvar+0]<t_obs0_strat[stratmap[currobs]])
		 feas*=0;
	 if(val[0]+feasmap[2*currobs*numvar+numvar+0]>t_obs0_strat[stratmap[currobs]])
		 feas*=0;
 }
 if(feas==1)
	 for(j=1; j<numvar; j++)
  		 if(j != var)
		 {
	/*
	Rprintf("Max: %4.0f Min: %4.0f ",val[j]+u[j],val[j]+v[j]);
	*/
    		 if((val[j]+feasmap[2*(currobs)*numvar+j])<t_obs[j])
				feas*=0;
			 if((val[j]+feasmap[2*(currobs)*numvar+numvar+j])>t_obs[j])
  				feas*=0;
		 }

 return(feas);
}

long int intpowmod(long int value, int ex, int q)
 {
 if ((ex==0) && (q!=1)) 
  return(1);
 else
  if ((ex==0) && (q==1))
   return(0);
  else
   return(((value % q)*intpowmod(value,ex-1,q)) % q);
 }

int hash_f(double *val, int numvar, int numobs, int q, int noint)
 {
 int value, j, incr;
 long int D;
 D=(long int)numobs;
 /*value=(int)val[0] % q;*/
 value=(int)val[noint] % q;
 
 if(numvar>noint)
	for(j=noint+1; j<numvar; j++)
/*		for(j=1; j<numvar; j++)*/
	{
		incr=(((int)val[j] % q)*intpowmod(D,j-1,q)) % q;
		value=(value+incr) % q;
	}
 /*
 Rprintf("hash_f returns %d\n",value);
   */
 return(value);
 }

void compare_and_add(double *val, long double counter, int IP, int noint,
                     int numobs, long double * newlist, int q)
 {
 /*
   Uses hashing with linear probe
 */

 int match=1;
 int j;
 long double * addr;




 addr=newlist+(IP+1)*hash_f(val, IP, numobs, q/(IP+1), noint);

 do
  {
  if(addr[IP]==0.0)
   {
   for(j=noint; j<IP; j++)
    addr[j]=(long double)val[j];
	addr[IP]=counter;
   match=1;
   /*
   Rprintf("Field was empty.");
     */
   }
  else
   {
   /*
   Rprintf("Collision.");
     */
   match=1;
   for(j=noint; j<IP; j++)
    if(addr[j]!=(long double)val[j])
     match*=0;
   if(match)
    {
    /*
    Rprintf("Match.");
      */
    addr[IP]+=counter;
    }
   else
    addr+=IP+1;
   if(addr>(newlist+q-1))
    addr=newlist;
   }
  } while(match==0);

 }

void clear_t0(struct t * head)
{
	struct t* point;
	point=head->next;
	while(point!=point->next)
	{
		point->value[0]=0.0;
		point=point->next;
	}
}

struct t* MSA_recursiv(int stage, double * t_obs, double * t_obs0_strat, double * array,
							  int var, int groups, int numobs, int numvar, int noint, int * stratmap,
							  double * feasmap)
 {
 struct t* oldhead;
 struct t* point;
 struct t* newhead;
 struct t* newtail;
 double *val;
 long double *newlist;
 int k, i, j, q, new_q;
 if (NULL == (val=(double *)R_alloc(numvar, sizeof(double))))
 {
	 error("no memory available\n");
 }
 if(NULL == (newhead=(struct t *)R_alloc(1, sizeof(struct t)))) /* oder mit hashing auf ein Feld */
 {
	 error("no memory available\n");
 }
 if (NULL == (newtail=(struct t *)R_alloc(1, sizeof(struct t))))
 {
	 error("no memory available\n");
 }
 newhead->next=newtail;
 newtail->next=newtail;
 new_q=0;
 if(stage==0)
  /* initialization */


  {
	/*
	Rprintf("Stage 0.\n");
	  */
  for(k=noint; k<=(int)y(array,0,numvar); k++)
	{
	for(j=0; j<numvar; j++)
	 val[j]=Data(array,0, j, numvar)*k;

/*	t_obs[0]=t_obs0_strat[stage];*/
	if(feasible(t_obs, t_obs0_strat, stratmap, val, var, numvar, noint, 0, k, feasmap))
	{
	 new_q++;
	 if (NULL == (point=(struct t *)R_alloc(1, sizeof(struct t))))
	 {
		 error("no memory available");
	 }
	 if (NULL == (point->value=(double *)R_alloc(numvar, sizeof(double))))
	 {
		 error("no memory available");
	 }
	 for(j=noint; j<numvar; j++)
	  point->value[j]=val[j];
	 point->counter=coco((int)y(array,0,numvar),k);
	 point->next=newhead->next;
	 newhead->next=point;
	 }
	}
  newhead->counter=(long double)new_q;
  }
 else
  {
  /* recursive call */
	
   /*	Rprintf("Stage %d.\n",stage); */
	  
  /* initialize hash table */


  oldhead=MSA_recursiv(stage-1, t_obs, t_obs0_strat, array, var, groups, numobs,
							  numvar, noint, stratmap, feasmap);
  if (oldhead == NULL)
   error("no memory available\n");
/*  Rprintf("Stage %d of %d, %d entries,\n",stage,groups,(int)oldhead->counter);
*/

  if(stage>0)
	  if(stratmap[stage]!=stratmap[stage-1])
		  clear_t0(oldhead);
  
  

  new_q=(int)oldhead->counter*((int)y(array,stage,numvar)+1)
		  *(numvar+1);
/*  
  Rprintf("      reserving %d entries for next stage.\n",
   new_q/(numvar+1), numvar+1);
   */
    
  if (NULL == (newlist=(long double *)R_alloc(new_q, sizeof(long double))))
  {
	  error("no memory available\n");
  }

  for(i=numvar; i<new_q; i+=numvar+1)
   newlist[i]=0.0;

  while(oldhead->next->next != oldhead->next)
   {
   for(k=0; k<=(int)y(array,stage,numvar); k++)
    {

    for(j=noint; j<numvar; j++)
     {
     val[j]=oldhead->next->value[j]+k*Data(array,stage,j,numvar);
/*     
     Rprintf("%2.0f ",val[j]);
  */     
     }

/*    
    Rprintf("%6.0f ... ",oldhead->next->counter);
*/    
    if(feasible(t_obs, t_obs0_strat, stratmap, val, var, numvar, noint, stage, k, feasmap))
	  {
      compare_and_add(val,
       oldhead->next->counter*coco((int)y(array,stage,numvar),k),
       numvar, noint, numobs, newlist, new_q);
/*   
     Rprintf("feasible.\n");
*/     
     }
/*  
    else
     Rprintf("infeasible.\n");
*/    

    }

   point=oldhead->next;

   oldhead->next=oldhead->next->next;

   }

  /* read data from hash table into list */

  q=new_q;
  new_q=0;
  /*
  Rprintf("Feld hatte groesse %d * %d\n",q/(numvar+1), numvar+1);
    */
  for(i=0; i<q; i+=(numvar+1))
   {
   /*
   Rprintf("q-i: %d %d",q,i);

   for(j=0; j<=numvar; j++)
    Rprintf("%4.0f ",newlist[i+j]);
   Rprintf("\n");
   */
   if(newlist[i+numvar] != 0.0)
    {
    if (NULL == (point=(struct t *)R_alloc(1, sizeof(struct t))))
	{
		error("no memory available\n");
	}
    if (NULL == (point->value=(double *)R_alloc(numvar, sizeof(double))))
	{
		error("no memory available\n");
	}
    point->next=newhead->next;
    newhead->next=point;
    new_q++;
    for(j=noint; j<numvar; j++)
     point->value[j]=(double)newlist[i+j];

    point->counter=newlist[i+numvar];
    }
    
   }
  newhead->counter=(long double)new_q;
  /*
  Rprintf("Liste hat Laenge %d\n",new_q);
  */

  }
 /*
 point=newhead->next;
 while(point->next != point)
  {
  for(j=0; j<numvar; j++)
   Rprintf("%2.0f ",point->value[j]);
  Rprintf("  %10.0f\n",point->counter);
  point=point->next;
  }
 */

 return(newhead);

 }

struct t* derecode(struct t* head, double * t_obs,
                   int *recoded, int numvar)
 {
 struct t *point;
 int j;

 if(head->next==head->next->next)
  Rprintf("Derecode got empty list\n");

 point=head->next;
 while(point != point->next)
  {
  for(j=1; j<numvar; j++)
   if(recoded[j-1])
    point->value[j]=t_obs[0]-point->value[j];
  point=point->next;
  }
 for(j=1; j<numvar; j++)
  if(recoded[j-1])
   t_obs[j]=t_obs[0]-t_obs[j];

 return(head);
 }

void calc_prob(struct t *head, double theta, int var)
 {
 /*
   calc_prob is designed so that an overflow is avoided
 */

 long double sum, summand, d1, d2, d3;
 struct t* out;
 struct t* in;

 out=head->next;
 while(out->next !=out)
  {
  sum=0.0;
  in=head->next;
  while(in->next != in)
   {
   d1=log(in->counter);
   d2=log(out->counter);
   d3=theta*(in->value[var]-out->value[var]);
   if((d1-d2+d3) > LPREC)
    summand=1.0/PREC;
   else
    summand=exp(d1-d2+d3);
   sum+=summand;
   in=in->next;
   }
  out->probability=1.0/sum;
  out=out->next;
  }
 }

/* void read_distribution(struct t* currhead,double *t_obs,
                       long double *qsc_obs, int var, int IP, int noint)
 {
// FILE *distribution;
 struct t *point;
 struct t *newer;
 int j, numvar, currvar, lastvar;

 bool variable = false;

 distribution=fopen(DISTFILE,"r");
 point=currhead;
 point->next=currhead->next;
 fscanf(distribution, "%d", &numvar);
 if(numvar != IP)
  Rprintf(
  "Warning! Distribution and Parameter files have not the same IP\n");
 if(read_comment(distribution)==0)
  Rprintf("End of expected comment not found\n");
 for(j=0; j<numvar; j++)
  {
  fscanf(distribution,"%f",&t_obs[j]);
  if(read_comment(distribution)==0)
   Rprintf("End of expected comment not found\n");
  // read_comment(distribution); 
  }
 currvar=-1;
 lastvar=var;

 while(!feof(distribution))
  {
  fscanf(distribution,"%d ", &currvar);
  if(lastvar == var)
   {
   if (NULL == (newer=(struct t *)R_alloc(1, sizeof(struct t))))
   {
     error("no memory available\n");
   }
   if (NULL == (newer->value=(double *)R_alloc(numvar, sizeof(double))))
   {
     error ("no memory available\n");
   }
   variable = true;
   }
  for(j=noint; j<numvar; j++)
   fscanf(distribution,"%f ",&(newer->value[j]));
  fscanf(distribution,"%f %f %f",
         &(newer->counter), &(newer->probability), &(newer->score));

  if(currvar==var)
   {
   if(t_obs[var]==newer->value[var])
    qsc_obs[var]=newer->score;
   newer->next=point->next;
   point->next=newer;
   point=point->next;   
   }
  lastvar=currvar;

  if (variable)
  {
	  variable = false;
  }  
 fclose(distribution);
 }
} */


void empty_list(struct t* head)
{
	struct t* point;
	point=head->next;
	while(point != point->next)
	{
		head->next=point->next;
		point=head->next;
	}
}
	   
void xlr(double *by_data_array, int *Pmaxstrat, int *PIP, int *Pnoint, int *Poptions, int *Ptst, int *Psc, int *Plr, int *Ppmid,
	int *Pvstart, int *Pvstop, double *Palpha, int * n_j, int *Pn_by, double *ciout, double *estout, double *distout, int *Pausgabe)
{

	int maxstrat = *Pmaxstrat; 
		int IP = *PIP; 
		int noint = *Pnoint; 
		int options = *Poptions; 
		int tst = *Ptst; 
		int sc = *Psc; 
		int lr = *Plr; 
		int pmid = *Ppmid;
	int vstart = *Pvstart; 
	int vstop = *Pvstop;
	double alpha = *Palpha;
	int n_by = *Pn_by;

	int ausgabe = *Pausgabe;

	int distzaehl = 0;

	/* int maxby = 1;  */

 int iby, groups, pmid_r, i, j, N=0, var=0, 
	 index_i, index_ii, bound, index_strat, all_groups=0, istrat;

 double z, chi;
 double *data_array;

 double * * allgrdataarray;
 int * groupsperstratum;
 int * stratmap=(int *)R_alloc(1, sizeof(int));

 
 double *t_obs;
 double *t_obs0_strat;
 double *groupdata;
 double *bygrpdata=(double *)R_alloc(1, sizeof(double));
 double * feasmap=(double *)R_alloc(1, sizeof(double));
 long double *qsc_obs;
 struct t *dis;
 double beta;
 int method;
 /*
 int stratq, currq;
 */
 int *recoded;
 double *limits;
 double *p_value;
/* 
 long double * stratlist;
*/

struct t *currhead;
struct t *currtail;
struct t *strathead;
struct t *strattail;
 
 if (NULL == (limits=(double *)R_alloc(2, sizeof(double))))
 {
	 error("no memory available\n");
 }
 if (NULL == (p_value=(double *)R_alloc(3, sizeof(double))))
 {
	 error("no memory available\n");
 }

 chi = -1;
 z = -1;

 if (ausgabe == 1)
 {
 Rprintf("Exact Logistic Regression\n\n");
 Rprintf("Version 20121105 1327\n");
 Rprintf("(c) Georg Heinze 1997-2000/2012\n\n");
 //Rprintf("Number of BY-groups:     %d\n",maxby);
 }
 
 /*
 Options:
 0 ... read data and write estimates to ESTFILE and CIFILE
 1 ... read data and write conditional distributions to file
       DISTFILE
       do not calculate estimates
 2 ... read data and write conditional distributions to file
       DISTFILE
       and write estimates to ESTFILE and CIFILE
 3 ... read conditional distribution from file
       DISTFILE
		     and write estimates to ESTFILE and CIFILE

 tst, sc, lr:
 0 ... do not calculate
 1 ... calculate

 pmid:
 0 ... no pmid
 1 ... pmid
 2 ... both
 */


if (NULL == (t_obs=(double *)R_alloc(IP, sizeof(double))))
{
	error("no memory available\n");
}
if (NULL == (t_obs0_strat=(double *)R_alloc(maxstrat, sizeof(double))))
{
	error ("no memory available\n");
}
if (NULL == (qsc_obs=(long double *)R_alloc(IP, sizeof(long double))))
{
	error("no memory available\n");
}
if (NULL == (recoded=(int *)R_alloc(IP, sizeof(int))))
{
	error("no memory available\n");
}

if (NULL == (strathead=(struct t*)R_alloc(1, sizeof(struct t))))
{
	error("no memory available\n");
}
if (NULL == (strattail=(struct t*)R_alloc(1, sizeof(struct t))))
{
	error("no memory available\n");
}
strathead->next=strattail;
strattail->next=strattail;
if (NULL ==(currhead=(struct t*)R_alloc(1, sizeof(struct t))))
{
	error("no memory available\n");
}

if (NULL == (currtail=(struct t*)R_alloc(1, sizeof(struct t))))
{
	error("no memory available\n");
}
currhead->next=currtail;
currtail->next=currtail;

 
for(j=noint; j<IP; j++) 
	 t_obs[j]=0;

if (ausgabe == 1)
{
Rprintf("Number of variables:    %d\n\n",IP-noint);
}


iby = 1;
/*Rprintf("Processing BY-group %d of %d.\n", iby,maxby);*/


   if (NULL == (allgrdataarray=(double * *)R_alloc(maxstrat, sizeof(double *))))
   {
	   error("no memory available\n");
   }
   if (NULL == (groupsperstratum=(int *)R_alloc(maxstrat, sizeof(int))))
   {
	   error("no memory available\n");
   }

  for(j=noint; j<IP; j++) 
   t_obs[j]=0;
  for(i=0; i<n_by; i++)
	  for(j=noint; j<IP; j++)
	  {
		t_obs[j]+=y(by_data_array,i,IP)*Data(by_data_array,i,j,IP);	  
	  }
  if(options<3)
  {
		  index_strat=0;
		  all_groups=0;
		  for(istrat=0; istrat<maxstrat; istrat++)
		  {
			  N=n_j[(iby-1)*maxstrat+istrat];
			  data_array=&by_data_array[index_strat];
			  index_strat+=N*(IP);

			  groupdata=group_data(data_array,IP,N,&groups);
			  if(groupdata == NULL)
				  return;
			  t_obs0_strat[istrat]=0.0;
			  for(i=0; i<N; i++)
				  t_obs0_strat[istrat]+=y(data_array,i,IP);;
			  groupdata=sort_groups(groupdata,IP,groups,var);

			  all_groups+=groups;
			  allgrdataarray[istrat]=groupdata;

			  groupsperstratum[istrat]=groups;
				        
		  }
		  if (NULL == (bygrpdata=(double *)R_alloc(all_groups*IP, sizeof(double))))
		  {
			  error("no memory available\n");
		  }
		  if (NULL == (stratmap=(int *)R_alloc(all_groups, sizeof(int))))
		  {
			  error("no memory availble\n");
		  }

		  index_i=0;
		  index_ii=0;
		  for(istrat=0; istrat<maxstrat; istrat++)
		  {
			  for(i=0; i<groupsperstratum[istrat]; i++)
			  {
				  for(j=1; j<IP; j++)
				  {
					  bygrpdata[index_i]=Data(allgrdataarray[istrat],i,j,IP);
					  index_i++;
				  }
				  bygrpdata[index_i]=y(allgrdataarray[istrat],i,IP);

				  
				  stratmap[index_ii]=istrat;
				  index_i++;
				  index_ii++;
			  }

		  }

		  for(istrat=maxstrat-1; istrat>=0; istrat--)

		  /*   make feasmap */
          if (NULL == (feasmap=(double *)R_alloc(all_groups*2*IP, sizeof(double))))
		  {
			  error("no memory available\n");
		  }
          for(bound=0; bound<2; bound++)
              for(i=all_groups-1; i>=0; i--)
                  for(j=noint; j<IP; j++)
                  {
                      if(i==all_groups-1)
                      {
                          feasmap[2*IP*i+bound*IP+j]=0.0;
                      }
                      else
                          if(stratmap[i]!=stratmap[i+1] && j==0)
                              feasmap[2*i*IP+bound*IP]=0.0;
                          else
                              feasmap[2*i*IP+bound*IP+j]=feasmap[2*(i+1)*IP+bound*IP+j]+y(bygrpdata, i+1, IP)*
                                  (Data(bygrpdata, i+1, j, IP)-(bound*2-1)*absolute(Data(bygrpdata, i+1, j, IP)))/2.0;
                  }
  }
  int line = 1;
  int estzaehl = 0;

  for(var=vstart; var<vstop+1; var++)
  {
	  if(options < 3)
	  {
		  currhead=MSA_recursiv(all_groups-1,t_obs, t_obs0_strat,bygrpdata,var,all_groups,N,IP,noint,
				  stratmap, feasmap);
	  }
   

   	  if(options>=1 && options<=2)
	  {
		 // distribution=fopen(DISTFILE,"a");
		  //distout[distzaehl] = iby;
		  //distzaehl++;
		  distout[distzaehl] = (double)(var);
		  distzaehl++;
	   	 // fprintf(distribution,"%4d %2d  ",iby,IP-noint);
		  for(j=noint; j<IP; j++)
		  {
			  	distout[distzaehl] = (t_obs[j]);
				distzaehl++;
			   //fprintf(distribution,"%2.0f ",t_obs[j]);
		  }
		  distout[distzaehl] = (-1.0);
		  distzaehl++;
		  //fprintf(distribution,"%20.0f", -1.0);
		  //distout[distzaehl] = (-1.0);
		  //distzaehl++;
		  //fprintf(distribution," %12.9f", -1.0);
		  //distout[distzaehl] = (-1.0);
		  //distzaehl++;
		  //fprintf(distribution," %10.4f\n", -1.0);
		  dis=currhead->next;
	   	  while(dis->next != dis)
		  {
			  //distout[distzaehl] = (iby);
		      //distzaehl++;
			  distout[distzaehl] = (double)(var);
		      distzaehl++;
			  //fprintf(distribution,"%4d %2d  ",iby,var);
			  for(j=0; j<IP; j++)
			  {
			  distout[distzaehl] = (double)(dis->value[j]);
		      distzaehl++;
				//  fprintf(distribution,"%2.0f ",dis->value[j]);
			  }
			  distout[distzaehl] = (double)(dis->counter);
		      distzaehl++;
			  //fprintf(distribution,"%20.0f", dis->counter);
			  //distout[distzaehl] = (double)(dis->probability);
		      //distzaehl++;
			  //fprintf(distribution," %12.9f", dis->probability);
			  //distout[distzaehl] = (double)(dis->score);
		      //distzaehl++;
			  //fprintf(distribution," %10.4f\n", dis->score);
			  dis=dis->next;
		  }
		 // fclose(distribution);
   	   }
   	   else
   	   {
	//	   if(options==3)
		//	   read_distribution(currhead,t_obs, qsc_obs, var, IP, noint);
/*
   		   for(j=noint; j<IP; j++)
   			   printf("Observed t%d: %3.0f\n", j, t_obs[j]);
*/

	   }

  /* Estimation */
 
	   if(options != 1)
 	   {

		   for(method=1; method<=3; method++)
		   {
			   beta=LogXact_estimates(t_obs, var, currhead, method);
		//	   estfile=fopen(ESTFILE,"a");
			   switch(method)
			   {
			   case 1:// fRprintf(estfile, "%d ECMUE %d %8.4f\n",iby,var,beta);
						estout[estzaehl] = var;
						estout[estzaehl+1] = beta;
						estzaehl += 2;
				   break;
			   case 2: //fRprintf(estfile, "%d ECMLE %d %8.4f\n",iby,var,beta);
					    estout[estzaehl] = var;
						estout[estzaehl+1] = beta;
						estzaehl += 2;
				   break;
			   case 3: //fRprintf(estfile, "%d LX    %d %8.4f\n",iby,var,beta);
						estout[estzaehl] = var;
						estout[estzaehl+1] = beta;
						estzaehl += 2;
				   break;
			   }
		//	   fclose(estfile);
		   }
 
 /* confidence limits and p-values */
/*		Rprintf("Confidence limits for variable %d...",var);
*/

		   for(method=1; method <=3; method++)
			   for(pmid_r=0; pmid_r<=1; pmid_r++)
			   {
				   if((pmid_r==pmid) || (pmid==2))
				   {
					 //  cifile=fopen(CIFILE,"a");
					   switch(method)
					   {
					   case 1: if(tst==1)
							   {
								   LXconf_limits(t_obs, limits, p_value, var, alpha, currhead, pmid_r);
								  // fRprintf(cifile,"%d TST",iby);
//								   Rprintf("-TST");
							   }
						   break;
					   case 2: if(sc==1)
							   {
								   chi=-1;
								   z=-1;
								   LXcl_scores(t_obs, limits, p_value, var, alpha, currhead, beta, pmid_r, &chi, &z);
								  // fRprintf(cifile,"%d SC",iby);
//								   Rprintf("-SC");
							   }
						   break;
					   case 3: if(lr==1)
							   {
								   LXcl_LR(t_obs, limits, var, alpha, currhead, beta, IP, pmid_r);
								//   fRprintf(cifile,"%d LR",iby);
//								   Rprintf("-LR");
							   }
						   break;
					   }
					   if((tst==1 && method==1) || (sc==1 && method==2) || (lr==1 && method==3))
					   {
						   if(pmid_r==1)
						   {
//							   Rprintf("-PMID");
							  // fRprintf(cifile,"-PMID");
						   }
						   else
						   {
							//   fRprintf(cifile,"\t");
						   }
						 //  fRprintf(cifile," \t%d \t%f \t%f \t%f \t%f \t%f \t%f \t%f\n", var, limits[0], limits[1], 
						//	   p_value[0], p_value[1], p_value[2], chi, z);

						   ciout[(line-1)*8+0] = var;
						   ciout[(line-1)*8+1] = limits[0];
						   ciout[(line-1)*8+2] = limits[1];
						   if (p_value[0] > 1.0)
						   {
							p_value[0] = 1.0;
						   }
                           ciout[(line-1)*8+3] = p_value[0];
                           ciout[(line-1)*8+4] = p_value[1];
                           ciout[(line-1)*8+5] = p_value[2];
                           ciout[(line-1)*8+6] = chi;
						   ciout[(line-1)*8+7] = z;

						   line++;
					   }
					 //  fclose(cifile);
				   } 
			   }	
	   }
	   empty_list(currhead);
//	   Rprintf("\n");
  } /* var loop */
 
 /* if(options > 0 && options < 3)
   fclose(distribution); */
}

