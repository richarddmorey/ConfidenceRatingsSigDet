#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


#define DECORRCRITS 1
#define SET 1

int pcl_ran_mvgaussian(gsl_rng *rseed, gsl_vector *mean, gsl_matrix *Sigma, gsl_vector *output)
{
  int i;
  gsl_linalg_cholesky_decomp(Sigma);
  for (i = 0; i < output->size; i++){
    gsl_vector_set (output, i , gsl_ran_gaussian(rseed,1));      
  }
  gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, Sigma, output);
  gsl_blas_daxpy (1.0,mean, output);
  return 0;
} 

unsigned long int random_seed()
{

  unsigned int seed;
  FILE *devrandom;
  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    fprintf(stderr,"Cannot open /dev/random, setting seed to 0\n");
    seed = 0;
  } else {
    fread(&seed,sizeof(seed),1,devrandom);
    fclose(devrandom);
  }

  return(seed);

}

double ran_inv_gamma(const gsl_rng * r, double a, double b){
  return(1/gsl_ran_gamma(r, a, 1/b));
}

double ran_trunc_norm_upper(const gsl_rng * r, double mu, double sigma, double a){
  double u = gsl_rng_uniform_pos(r);
  double v = gsl_cdf_ugaussian_P((-a+mu)/sigma);
  return(-sigma * gsl_cdf_ugaussian_Pinv(u*v)+mu);
}

double ran_trunc_norm_lower(const gsl_rng * r, double mu, double sigma, double a){
  double u = gsl_rng_uniform_pos(r);
  double v = gsl_cdf_ugaussian_P((a-mu)/sigma);
  return(sigma * gsl_cdf_ugaussian_Pinv(u*v)+mu);
}


double ran_trunc_norm_both(const gsl_rng * r, double mu, double sigma, double a, double b){
  double low = gsl_cdf_ugaussian_P((a-mu)/sigma);
  double high = gsl_cdf_ugaussian_P((b-mu)/sigma);
  double u = gsl_ran_flat(r,low,high);
  return(sigma * gsl_cdf_ugaussian_Pinv(u)+mu);
}

int main(int argc, char *argv[]){
  
  unsigned long int seed;
  seed = random_seed();
  int ITERS,iter,i,j,k,K,I,J,nRespCat;
  char *fnEstimates,*fnChains,*fnInput;
  int yij,sij;
  int errs=0;
  char *name;

  double a,b,e,f,sigDecorr,mu_crits,sd_crits;
  double sig2mu,mumu;

  ////This is all data and prior entry code.

  if(argc<9){
    fprintf(stderr,"Arguments needed: #iterations nRespCat sig2_mu a b e f sigDecorr [analysis.name] [input.file]\n");
    fprintf(stderr,"\n****Not enough arguments.\n\n");
    return(1);
  }
  if(argc>9){
    name=argv[9];
  }else{
    name="analysis";
  }

  ITERS=atoi(argv[1]);
  nRespCat=atoi(argv[2]);
  sig2mu=atof(argv[3]);  
  a=atof(argv[4]);
  b=atof(argv[5]);
  e=atof(argv[6]);
  f=atof(argv[7]);
  sigDecorr=atof(argv[8]);

  //confirm and check command line arguments.
  printf("Name      : %s\n",name);
  printf("nRespCat = %i\n",nRespCat);
  printf("sig2mu   = %.2f\n",sig2mu);
  printf("a        = %.2f\n",a);
  printf("b        = %.2f\n",b);
  printf("e        = %.2f\n",e);
  printf("f        = %.2f\n",f);   
  printf("sigDecorr= %.2f\n",sigDecorr);
  

  if(nRespCat<=0){printf("\n*****nRespCat must be greater than 0.\n");errs++;}
  if(sig2mu<=0){printf("\n*****sig2mu must be greater than 0.\n");errs++;}
  if(a<=0){printf("\n*****a must be greater than 0.\n");errs++;}
  if(b<=0){printf("\n*****b must be greater than 0.\n");errs++;}
  if(e<=0){printf("\n*****e must be greater than 0.\n");errs++;}
  if(f<=0){printf("\n*****f must be greater than 0.\n");errs++;}
  if(sigDecorr<=0){printf("\n*****sigDecorr must be greater than 0.\n");errs++;}

  if(errs>0){printf("Exiting...\n\n");exit(1);}

  asprintf(&fnEstimates,"%s.est",name);
  asprintf(&fnChains,"%s.chn",name);

  if(argc>10){
    fnInput=argv[10];
  }else{
    asprintf(&fnInput,"%s.dat",name);
  }

  K=nRespCat-1;


  //open files.
  FILE *CHAINS;
  if ((CHAINS = fopen(fnChains,"w")) == NULL) {
    fprintf(stderr,"\n****Cannot open chain output file!\n\n");
    return(1);
  }

  FILE *ESTIMATES;
  if ((ESTIMATES = fopen(fnEstimates,"w")) == NULL) {
    fprintf(stderr,"\n****Cannot open estimate output file!\n\n");
    fclose(CHAINS);
    return(1);
  }

  FILE *INPUT;
  if ((INPUT = fopen(fnInput,"r")) == NULL) {
    fprintf(stderr,"\n****Cannot open input file %s!\n\n",fnInput);
    fclose(CHAINS);
    fclose(ESTIMATES);
    return(1);
  }
  
  //Now that we have opened the files, get the data.
  
  i=0;
  int newline=0,spaces=0;
  int got;
  while(!newline&&!feof(INPUT)){
    i++;
    got=fgetc(INPUT);
    //printf("Character %i: %i %c\n",i,got,got);
    if(got==32) spaces++;
    if(got==10) newline=1;
  }
  rewind(INPUT);

  if(!newline){
    printf("Could not determine number of conditions!\n");
    exit(1);
  }
  J=(spaces+1)/2;

  i=0;
  int totsig=0;
  int totnoi=0;

  printf("Reading data...\n");
  while(!feof(INPUT)){
    for(j=0;j<J;j++){
      if(j==(J-1)){
        fscanf(INPUT,"%i %i\n",&yij,&sij);
      }else{
        fscanf(INPUT,"%i %i ",&yij,&sij);
      }
      if((yij>=nRespCat)||(sij!=0&&sij!=1)){
        printf("\n****Invalid data for participant %i, item %i:  (yij=%i,sij=%i)\n\n",i+1,j+1,yij,sij);
        fclose(CHAINS);
        fclose(ESTIMATES);
        fclose(INPUT);
        return(1);
      }
      if(yij!=-1){
	totsig+=sij; 
	totnoi+=1-sij;
      }
    }
    i++;
  }
  rewind(INPUT);

  I=i;

  //printf("Read %d subjects and %d items, with %d total signal trials.\n",I,J,totsig);

  //initialize matrices
  
  printf("Reserving memory for matrices...\n");
  int Ys[totsig],Yn[totnoi],Subs[totsig],Subn[totnoi];
  gsl_vector *Ws,*Wn,*Ts,*Tn,*WsMean,*WnMean,*WsDev,*TsMean,*TnMean;
  gsl_matrix *Xs,*Xn,*Ss,*Sn,*TpSs,*TpSn,*XstXs,*XntXn,*Vs,*Vn;
  Xs = gsl_matrix_calloc(totsig, I+J+1);
  Xn = gsl_matrix_calloc(totnoi, I+J+1);
  XstXs = gsl_matrix_calloc(I+J+1, I+J+1);
  XntXn = gsl_matrix_calloc(I+J+1, I+J+1);
  Vs = gsl_matrix_calloc(I+J+1, I+J+1);
  Vn = gsl_matrix_calloc(I+J+1, I+J+1);
  gsl_matrix *Vn0 = gsl_matrix_calloc(I+J+1, I+J+1);
  gsl_matrix *Vs0 = gsl_matrix_calloc(I+J+1, I+J+1);
  Ws = gsl_vector_calloc(totsig);
  Wn = gsl_vector_calloc(totnoi);  
  Ts = gsl_vector_calloc(I+J+1);
  Tn = gsl_vector_calloc(I+J+1);
  WsMean = gsl_vector_calloc(totsig);
  WnMean = gsl_vector_calloc(totnoi);  
  WsDev = gsl_vector_calloc(totsig);
  TsMean= gsl_vector_calloc(I+J+1);
  TnMean= gsl_vector_calloc(I+J+1);
  gsl_vector *TsMean0= gsl_vector_calloc(I+J+1);
  gsl_vector *TnMean0= gsl_vector_calloc(I+J+1);
    


  gsl_vector *isAlph = gsl_vector_calloc(I+J+1);
  gsl_vector *isBeta = gsl_vector_calloc(I+J+1);
  
  TpSs = gsl_matrix_calloc(I+J+1,I+J+1);
  TpSn = gsl_matrix_calloc(I+J+1,I+J+1);
  gsl_matrix_set_identity(TpSs);
  gsl_matrix_set_identity(TpSn);

  Ss = gsl_matrix_calloc(totsig, totsig);
  gsl_matrix_set_identity(Ss);
  Sn = gsl_matrix_calloc(totnoi, totnoi);
  gsl_matrix_set_identity(Sn);
 
  gsl_permutation *gslPerm = gsl_permutation_alloc(I+J+1); 
  gsl_permutation_init(gslPerm);
 
  gsl_matrix_set(TpSs,0,0,1.0/sig2mu);
  gsl_matrix_set(TpSn,0,0,1.0/sig2mu);

 
  //end initialize matrices

  printf("Creating design matrix...\n");
  int y[I][J],s[I][J],sigindex=0,noiindex=0;
  i=0;
  while(!feof(INPUT)){
    for(j=0;j<J;j++){
      if(j==(J-1)){
	fscanf(INPUT,"%i %i\n",&yij,&sij);
      }else{
	fscanf(INPUT,"%i %i ",&yij,&sij);
      }
      y[i][j]=yij;
      s[i][j]=sij;
      if(yij!=-1){
	if(sij){
	  Ys[sigindex] = y[i][j];
	  Subs[sigindex] = i;
	  gsl_matrix_set(Xs,sigindex,0,1);
	  gsl_matrix_set(Xs,sigindex,i+1,1);
	  gsl_matrix_set(Xs,sigindex++,I+1+j,1);
	}else{
	  Yn[noiindex] = y[i][j];
	  Subn[noiindex] = i;
	  gsl_matrix_set(Xn,noiindex,0,1);
	  gsl_matrix_set(Xn,noiindex,i+1,1);
	  gsl_matrix_set(Xn,noiindex++,I+1+j,1);
	}
      }
    }
    i++;
  }
  
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xn,Xn,1,XntXn);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,Xs,Xs,1,XstXs);
 

  //gsl_vector_set(Tn,0,-.6);
  //gsl_vector_set(Ts,0,.6);

  
  //gsl_matrix_fprintf(stdout,Xs,"%f");

  //exit(0);

  printf("Read %i subjects values for %i items.\n\n",I,J);
  fclose(INPUT);
  ////End data entry.

  //Create starting values
  int keq0=0;
  double crits[I][K],critCand[I][K],w[I][J];
  double maxw[I][K+1];
  double minw[I][K+1]; 
  double sig2=1;  
  double sumWnSqr=0,sumWsSqr=0;
  int nSignal[I],totSignal=0;
  int nRespS[I][K+1];
  int nRespN[I][K+1];
  double sumwN[I];
  double sumwS[I];
  double zDecorr[I][K];
  double sig2N=1;
  double sig2beta0=.2,sig2beta1=.2,sig2alph0=.2,sig2alph1=.2;

  double bDecorr[I];
  int accDecorr[I],badCrits[I];
  double decorrRate=0,sumAlp20=0,sumAlp21=0,sumBet20=0,sumBet21=0;

  

  for(i=0;i<I;i++){
    accDecorr[i]=0;
    nSignal[i]=0;
    nRespN[i][K]=0;
    nRespS[i][K]=0;
    crits[i][0]=0;
    crits[i][K-1]=SET;
    for(k=0;k<K;k++){
      nRespS[i][k]=0;
      nRespN[i][k]=0;
      if(k!=0 && k<(K-1)){
	crits[i][k]=k*SET/((float)(K-1));
	//printf("i: %i k: %i crit: %f\n",i,k,crits[i][k]);
      }
    }
    for(j=0;j<J;j++){
      nSignal[i]+=s[i][j];
      totSignal+=s[i][j];
      if(y[i][j]!=-1){
	nRespN[i][y[i][j]]+=1-s[i][j];
	nRespS[i][y[i][j]]+=s[i][j];
      }
    }
  }
  
 


  //for(k=0;k<(K+1);k++)
  //{
  //  printf("nRespN[0][%d]=%d\n",k,nRespN[0][k]);
  //  printf("nRespS[0][%d]=%d\n",k,nRespS[0][k]);
  //}
  //end starting values
  

  //Initialize GSL
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  //gsl_rng_set(r,seed);
  gsl_rng_set(r,0);


  gsl_blas_dgemv(CblasNoTrans,1,Xs,Ts,0,WsMean);
  gsl_blas_dgemv(CblasNoTrans,1,Xn,Tn,0,WnMean);  

  
  printf("Total signal trials: %i\nTotal noise trials: %i\n\n",totsig,totnoi);
  
  //begin MCMC loop.
  i=0;
  printf("Starting %i MCMC iterations...\n",ITERS);
  for(iter=0;iter<ITERS;iter++){
    
    for(i=0;i<I;i++){
      sumwS[i]=0;
      sumwN[i]=0;
      for(k=0;k<(K+1);k++){
	minw[i][k]=GSL_POSINF;
	maxw[i][k]=GSL_NEGINF;
      }
    }
    //sig2=1;

    //sample latent variables
    
    for(i=0;i<totsig;i++){
      //printf("%d %d %f %f\n",Ys[i],Subs[i],crits[Subs[i]][0],crits[Subs[i]][1]);
      if( Ys[i]==0 ){
	gsl_vector_set(Ws,i,ran_trunc_norm_lower(r, gsl_vector_get(WsMean,i), sqrt(sig2), crits[Subs[i]][0]));
	if(gsl_vector_get(Ws,i)>maxw[Subs[i]][0]) { maxw[Subs[i]][0]=gsl_vector_get(Ws,i);}
      }else if( Ys[i]==K ){
      	gsl_vector_set(Ws,i,ran_trunc_norm_upper(r, gsl_vector_get(WsMean,i), sqrt(sig2), crits[Subs[i]][K-1]));
	if(gsl_vector_get(Ws,i)<minw[Subs[i]][K]) { minw[Subs[i]][K]=gsl_vector_get(Ws,i);}
      }else{
	gsl_vector_set(Ws,i,ran_trunc_norm_both(r, gsl_vector_get(WsMean,i), sqrt(sig2), crits[Subs[i]][Ys[i]-1],crits[Subs[i]][Ys[i]]));
	if(gsl_vector_get(Ws,i)>maxw[Subs[i]][Ys[i]]) { maxw[Subs[i]][Ys[i]]=gsl_vector_get(Ws,i);}
	if(gsl_vector_get(Ws,i)<minw[Subs[i]][Ys[i]]) { minw[Subs[i]][Ys[i]]=gsl_vector_get(Ws,i);}
      }
    }
    for(i=0;i<totnoi;i++){
      //if(i==511) printf("%d %d %f %f\n",Yn[i],Subn[i],crits[Subn[i]][Yn[i]-1],crits[Subs[i]][Yn[i]]);
      if( Yn[i]==0 ){
	gsl_vector_set(Wn,i,ran_trunc_norm_lower(r, gsl_vector_get(WnMean,i), sqrt(sig2N), crits[Subn[i]][0]));
	if(gsl_vector_get(Wn,i)>maxw[Subn[i]][0]) { maxw[Subn[i]][0]=gsl_vector_get(Wn,i);}
      }else if( Yn[i]==K ){
      	gsl_vector_set(Wn,i,ran_trunc_norm_upper(r, gsl_vector_get(WnMean,i), sqrt(sig2N), crits[Subn[i]][K-1]));
	if(gsl_vector_get(Wn,i)<minw[Subn[i]][K]) { minw[Subn[i]][K]=gsl_vector_get(Wn,i);}
      }else{
	gsl_vector_set(Wn,i,ran_trunc_norm_both(r, gsl_vector_get(WnMean,i), sqrt(sig2N), crits[Subn[i]][Yn[i]-1],crits[Subn[i]][Yn[i]]));
	if(gsl_vector_get(Wn,i)>maxw[Subn[i]][Yn[i]]) { maxw[Subn[i]][Yn[i]]=gsl_vector_get(Wn,i);}
	if(gsl_vector_get(Wn,i)<minw[Subn[i]][Yn[i]]) { minw[Subn[i]][Yn[i]]=gsl_vector_get(Wn,i);}
      }
    }
    //printf("***\n");
    //gsl_vector_fprintf(stderr,Wn,"%f");
    //gsl_vector_fprintf(stderr,WnMean,"%f");
    //for(i=0;i<totnoi;i++) printf("%i : %f\n",i,gsl_vector_get(WnMean,i));


    
    gsl_blas_daxpy(-1,Ws,WsMean);
    gsl_blas_daxpy(-1,Wn,WnMean);

    gsl_blas_ddot(WsMean,WsMean,&sumWsSqr);
    gsl_blas_ddot(WnMean,WnMean,&sumWnSqr);

    //gsl_vector_fprintf(stderr,WnMean,"%f");
    //for(i=0;i<totnoi;i++) printf("%i : %f\n",i,gsl_vector_get(Wn,i));

    //printf("\nSS: %f,%f,%f,%d\n--------\n",sumWnSqr,a,b,totnoi);
    

    //sample sigma^2
    sig2=ran_inv_gamma(r, a+.5*totsig,b+(.5*sumWsSqr));
    //sig2=1;
    //printf("IG: %f   %f\n",a0+.5*totSignal,(1.0)*b0+(.5*sumWSqr));
    //sig2=1.3;
    fwrite(&sig2,sizeof(double),1,CHAINS);

    sig2N=ran_inv_gamma(r, a+.5*totnoi,b+(.5*sumWnSqr));

    fwrite(&sig2N,sizeof(double),1,CHAINS);

    // set up variances for full conditional on parameter vector 
    for(i=0;i<(I+J+1);i++){
      if(i>I){
	gsl_matrix_set(TpSs,i,i,1.0/sig2beta1);
	gsl_matrix_set(TpSn,i,i,1.0/sig2beta0);
      }else if(i>0){	
  	gsl_matrix_set(TpSs,i,i,1.0/sig2alph1);
	gsl_matrix_set(TpSn,i,i,1.0/sig2alph0);
      }
      //printf("i: %d  |  %f\n",i,gsl_matrix_get(TpSs,i,i));
    }
    
   
    
    gsl_matrix_memcpy(Vs,TpSs);
    gsl_matrix_memcpy(Vn,TpSn);
    
    gsl_matrix_scale(Vs,sig2);
    gsl_matrix_scale(Vn,sig2N);
    gsl_matrix_add(Vn,XntXn);
    gsl_matrix_add(Vs,XstXs);
    gsl_matrix_scale(Vs,1.0/sig2);
    gsl_matrix_scale(Vn,1.0/sig2N);

    //for(i=0;i<(I+J+1);i++){
    // for(j=0;j<(I+J+1);j++){
    //	printf("%f ",gsl_matrix_get(Vs,i,j));
    //  }
    //  printf("\n");
    //}


    gsl_linalg_LU_decomp(Vn,gslPerm,&i);
    gsl_linalg_LU_invert(Vn,gslPerm,Vn0);
    
    gsl_linalg_LU_decomp(Vs,gslPerm,&i);
    gsl_linalg_LU_invert(Vs,gslPerm,Vs0);

    //for(i=0;i<(I+J+1);i++){
    //  for(j=0;j<(I+J+1);j++){
    //	printf("%f ",gsl_matrix_get(Vs0,i,j));
    //  }
    //  printf("\n");
    //}

    

    //sample parameter vectors
      
    gsl_blas_dgemv(CblasTrans,1.0/sig2,Xs,Ws,0,TsMean);
    gsl_blas_dgemv(CblasNoTrans,1,Vs0,TsMean,0,TsMean0);

    gsl_blas_dgemv(CblasTrans,1.0/sig2N,Xn,Wn,0,TnMean);
    gsl_blas_dgemv(CblasNoTrans,1,Vn0,TnMean,0,TnMean0);

    pcl_ran_mvgaussian(r, TsMean0, Vs0, Ts);
    pcl_ran_mvgaussian(r, TnMean0, Vn0, Tn);

 
    gsl_vector_fwrite(CHAINS,Ts);
    gsl_vector_fwrite(CHAINS,Tn);
   

    //Sample variances of parameters
    sumAlp20=0;
    sumAlp21=0;
    sumBet20=0;
    sumBet21=0;
    for(i=0;i<(I+J+1);i++){
      if(i>I){
	sumBet20+=pow(gsl_vector_get(Tn,i),2);
	sumBet21+=pow(gsl_vector_get(Ts,i),2);
      }else if(i>0){	
	sumAlp20+=pow(gsl_vector_get(Tn,i),2);
	sumAlp21+=pow(gsl_vector_get(Ts,i),2);
      }
    }
    //printf("%f\n",sumBet20);
    
    //sample sigma^2
    sig2alph0=ran_inv_gamma(r, e+.5*I,f+(.5*sumAlp20));
    sig2alph1=ran_inv_gamma(r, e+.5*I,f+(.5*sumAlp21));
    sig2beta0=ran_inv_gamma(r, e+.5*J,f+(.5*sumBet20));
    sig2beta1=ran_inv_gamma(r, e+.5*J,f+(.5*sumBet21));

    //printf("IG: %f   %f\n",a0+.5*totSignal,(1.0)*b0+(.5*sumWSqr));
    //sig2=1.3;
    fwrite(&sig2alph0,sizeof(double),1,CHAINS);
    fwrite(&sig2alph1,sizeof(double),1,CHAINS);
    fwrite(&sig2beta0,sizeof(double),1,CHAINS);
    fwrite(&sig2beta1,sizeof(double),1,CHAINS);

    gsl_blas_dgemv(CblasNoTrans,1,Xs,Ts,0,WsMean);
    gsl_blas_dgemv(CblasNoTrans,1,Xn,Tn,0,WnMean);  

    //gsl_vector_fprintf(stderr,WsMean,"%f");
    //printf("\n\n");
    //gsl_vector_fprintf(stderr,WnMean,"%f");
    //printf("\n\n");
    //gsl_vector_fprintf(stderr,WsMean,"%g");

    
    //printf("\ns2=%f\n%f %f,%f %f\n",sig2,sumwN[i]/(1.0*(J-nSignal[i])),1.0/sqrt(1.0*(J-nSignal[i])),sumwS[i]/(1.0*nSignal[i]),sqrt(sig2/(1.0*(nSignal[i]))));

      //sample criteria
    for(i=0;i<I;i++){
      //zDecorr[i]=gsl_ran_gaussian(r,sigDecorr);
      bDecorr[i]=1;
      badCrits[i]=0;
      crits[i][0]=0; crits[i][K-1]=SET;
      for(k=0;k<K;k++){
	zDecorr[i][k]=gsl_ran_gaussian(r,sigDecorr);
	//printf("i=%d k=%i max=%f min=%f\n",i,k,maxw[i][k],minw[i][k+1]);
	//if((k!=0) && (k!=(K-1))) crits[i][k]=(gsl_rng_uniform_pos(r)*(minw[i][k+1]-maxw[i][k])+maxw[i][k]);
	//crits[i][k]=ran_trunc_norm_both(r,mu_crits,sd_crits,maxw[i][k],minw[i][k+1])*(k!=keq0);
	//printf("newcrit: %f\n",crits[i][k]);
	//for decorrelating
	critCand[i][k]=crits[i][k]+zDecorr[i][k]*(k!=0 && k!=(K-1));
	//bDecorr[i]=bDecorr[i]*exp(-.5*pow((critCand[i][k]-mu_crits)/sd_crits,2)+.5*pow((crits[i][k]-mu_crits)/sd_crits,2));
	//printf("candcrit: %f\n",critCand[i][k]);
	if(k>0 && critCand[i][k]<critCand[i][k-1] ) badCrits[i]=1;
	//printf("CRITS: %d %d d:%f s:%f  :   %f %f\n",iter,k,ds[i],sig2,crits[i][k],critCand[i][k]);
      }       
    }
    
    if(DECORRCRITS){
      for(i=0;i<totsig;i++){
	if(!badCrits[Subs[i]])
	  if(Ys[i]==0){
	    bDecorr[Subs[i]]=bDecorr[Subs[i]]*(
					       gsl_cdf_ugaussian_P((critCand[Subs[i]][0]-gsl_vector_get(WsMean,i))/sqrt(sig2)) /
					       gsl_cdf_ugaussian_P((crits[Subs[i]][0]   -gsl_vector_get(WsMean,i))/sqrt(sig2))
					       );
	  }else if(Ys[i]==K){
	    bDecorr[Subs[i]]=bDecorr[Subs[i]]*(
					       (1-gsl_cdf_ugaussian_P((critCand[Subs[i]][K-1]-gsl_vector_get(WsMean,i))/sqrt(sig2))) /
					       (1-gsl_cdf_ugaussian_P((crits[Subs[i]][K-1]   -gsl_vector_get(WsMean,i))/sqrt(sig2)))
					       );
	  }else{
	    bDecorr[Subs[i]]=bDecorr[Subs[i]]*(
					       (gsl_cdf_ugaussian_P((critCand[Subs[i]][Ys[i]]-gsl_vector_get(WsMean,i))/sqrt(sig2)) - gsl_cdf_ugaussian_P((critCand[Subs[i]][Ys[i]-1]-gsl_vector_get(WsMean,i))/sqrt(sig2))) /
					       (gsl_cdf_ugaussian_P((crits[Subs[i]][Ys[i]]   -gsl_vector_get(WsMean,i))/sqrt(sig2)) - gsl_cdf_ugaussian_P((crits[Subs[i]][Ys[i]-1]   -gsl_vector_get(WsMean,i))/sqrt(sig2)))
					       );
	  }
      }
      
      for(i=0;i<totnoi;i++){
	if(!badCrits[Subn[i]])
	  if(Yn[i]==0){
	    bDecorr[Subn[i]]=bDecorr[Subn[i]]*(
					       gsl_cdf_ugaussian_P((critCand[Subn[i]][0]-gsl_vector_get(WnMean,i))/sqrt(sig2N)) /
					       gsl_cdf_ugaussian_P((crits[Subn[i]][0]   -gsl_vector_get(WnMean,i))/sqrt(sig2N))
					       );
	  }else if(Yn[i]==K){
	    bDecorr[Subn[i]]=bDecorr[Subn[i]]*(
					       (1-gsl_cdf_ugaussian_P((critCand[Subn[i]][K-1]-gsl_vector_get(WnMean,i))/sqrt(sig2N))) /
					       (1-gsl_cdf_ugaussian_P((crits[Subn[i]][K-1]   -gsl_vector_get(WnMean,i))/sqrt(sig2N)))
					       );
	  }else{
	    bDecorr[Subn[i]]=bDecorr[Subn[i]]*(
					       (gsl_cdf_ugaussian_P((critCand[Subn[i]][Yn[i]]-gsl_vector_get(WnMean,i))/sqrt(sig2N)) - gsl_cdf_ugaussian_P((critCand[Subn[i]][Yn[i]-1]-gsl_vector_get(WnMean,i))/sqrt(sig2N))) /
						(gsl_cdf_ugaussian_P((crits[Subn[i]][Yn[i]]   -gsl_vector_get(WnMean,i))/sqrt(sig2N)) - gsl_cdf_ugaussian_P((crits[Subn[i]][Yn[i]-1]   -gsl_vector_get(WnMean,i))/sqrt(sig2N)))
					       );
	  }
      }
      
      //for(i=0;i<I;i++)
      //for(k=0;k<K;k++)
      //printf("m: %d, i: %d, k: %d,  z: %f crit: %f  cand: %f, BAD: %d\n",iter,i,k,zDecorr[i][k],crits[i][k],critCand[i][k],badCrits[i]);
      
 
      for(i=0;i<I;i++){
	if(!badCrits[i]){
	  
	  accDecorr[i]=gsl_ran_bernoulli(r,(bDecorr[i]>1)?1:bDecorr[i]);
	  //printf("BDECORR: m%d i%d    bdecorr %f\n",iter,i,bDecorr[i]);
	  if(accDecorr[i]){
	    decorrRate+=1/(1.0*I*ITERS);
	    
	    for(k=0;k<K;k++){
	      //printf("m: %d, i: %d, k: %d,  z: %f crit: %f  cand: %f\n",iter,i,k,zDecorr[i][k],crits[i][k],critCand[i][k]);
	      crits[i][k]=critCand[i][k];
	    }
	  }
	}
      }
    }
    
    fwrite(crits,sizeof(double),I*K,CHAINS);
    

    if(!((iter+1)%25)) { printf("Iteration %i\n",iter+1); }
      
    
  }//End MCMC loop
  

  printf("\nDone. Decorrelating step acceptance rate: %f\n\n",decorrRate);
  //close files
  fclose(ESTIMATES);
  fclose(CHAINS);
  return 0;
}
