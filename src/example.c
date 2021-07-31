#include "d3_qpgf_d1.h"

void normal_sum(double complex *qG,double complex *dqG,double k,double kx,double d,double *r);


int main()
{
  QPDT qd;
  double lambda,k,theta,phi,kx,r[3],d,eps;
  double complex qG,dqG[3];
  int erc; // returned code

  eps=1.0e-15; // relative error requested
  lambda=1.0;  // wave length
  theta=0.3;   // for wave vector 
  phi=0.5;     // kx=k*sin(theta)*cos(phi), ky=k*sin(theta)*sin(phi), kz=k*cos(theta)
  d=0.5;       // periodic number, lattice vector is (d,0,0)

  k=2.0*M_PI/lambda;
  kx=k*sin(theta)*cos(phi);
  qd.k=k;
  qd.d=d;
  qd.kx=kx;
  printf("eps= %g\n",eps);
  printf("k  = %g\n",k);
  printf("k_x= %g\n",kx);
  printf("d  = %g\n\n",d);

  r[0]=-2.1;
  r[1]=0.3;
  r[2]=-0.35;
  printf("r  =(%5g, %5g,%5g)\n",r[0],r[1],r[2]);

  printf("--- normal summation ---\n");
  normal_sum(&qG,dqG,k,kx,d,r);

  erc=d3hm_qpgf_d1_fd(&qG,dqG,r,eps,&qd);
  printf("--- Fourier domain summation ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  erc=d3hm_qpgf_d1_ew(&qG,dqG,r,eps,&qd);
  printf("--- Ewald method ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  erc=d3hm_qpgf_d1(&qG,dqG,r,eps,&qd);
  printf("--- Auto ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  r[0]=0.1;
  r[1]=0.04;
  r[2]=0.01;
  printf("\nr  =(%5g, %5g,%5g)\n",r[0],r[1],r[2]);

  printf("--- normal summation ---\n");
  normal_sum(&qG,dqG,k,kx,d,r);

  erc=d3hm_qpgf_d1_fd(&qG,dqG,r,eps,&qd);
  printf("--- Fourier domain summation ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  erc=d3hm_qpgf_d1_ew(&qG,dqG,r,eps,&qd);
  printf("--- Ewald method ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  erc=d3hm_qpgf_d1(&qG,dqG,r,eps,&qd);
  printf("--- Auto ---\n");
  printf("qG    =% 15.14e %+15.14e I,erc=%d\n",creal(qG),cimag(qG),erc);
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));

  return 0;
}


void normal_sum(double complex *qG,double complex *dG,double k,double kx,double d,double *r)
{
  double ar;
  int lmax,i;

  lmax=1000000;

  ar=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  *qG=cexp(I*k*ar)/ar;
  dG[0]=(I*k*ar-1.0)*r[0]/(ar*ar*ar)*cexp(I*k*ar);
  dG[1]=(I*k*ar-1.0)*r[1]/(ar*ar*ar)*cexp(I*k*ar);
  dG[2]=(I*k*ar-1.0)*r[2]/(ar*ar*ar)*cexp(I*k*ar);
  for(i=1;i<lmax;i++){
  	ar=sqrt(pow(r[0]+(double)i*d,2)+pow(r[1],2)+pow(r[2],2));
  	*qG+=cexp(I*k*ar)/ar*cexp(+I*(double)i*kx*d);
  	dG[0]+=(I*k*ar-1.0)*(r[0]+(double)i*d)/(ar*ar*ar)*cexp(I*k*ar)*cexp(+I*(double)i*kx*d);
  	dG[1]+=(I*k*ar-1.0)*r[1]/(ar*ar*ar)*cexp(I*k*ar)*cexp(+I*(double)i*kx*d);
  	dG[2]+=(I*k*ar-1.0)*r[2]/(ar*ar*ar)*cexp(I*k*ar)*cexp(+I*(double)i*kx*d);

  	ar=sqrt(pow(r[0]-(double)i*d,2)+pow(r[1],2)+pow(r[2],2));
  	*qG+=cexp(I*k*ar)/ar*cexp(-I*(double)i*kx*d);
  	dG[0]+=(I*k*ar-1.0)*(r[0]-(double)i*d)/(ar*ar*ar)*cexp(I*k*ar)*cexp(-I*(double)i*kx*d);
  	dG[1]+=(I*k*ar-1.0)*r[1]/(ar*ar*ar)*cexp(I*k*ar)*cexp(-I*(double)i*kx*d);
  	dG[2]+=(I*k*ar-1.0)*r[2]/(ar*ar*ar)*cexp(I*k*ar)*cexp(-I*(double)i*kx*d);

  	if(i%(lmax/2)==0){
  		printf("|l|=%d\n",i);
  		printf("qG    =% 15.14e %+15.14e I\n",creal(*qG)/(4.0*M_PI),cimag(*qG)/(4.0*M_PI));
  		printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dG[0])/(4.0*M_PI),cimag(dG[0])/(4.0*M_PI));
  		printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dG[1])/(4.0*M_PI),cimag(dG[1])/(4.0*M_PI));
  		printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dG[2])/(4.0*M_PI),cimag(dG[2])/(4.0*M_PI));
  	}
  }

  *qG/=4.0*M_PI;
  dG[0]/=4.0*M_PI;
  dG[1]/=4.0*M_PI;
  dG[2]/=4.0*M_PI;

  printf("|l|=%d\n",i);
  printf("qG    =% 15.14e %+15.14e I\n",creal(*qG),cimag(*qG));
  printf("dqG/dx=% 15.14e %+15.14e I\n",creal(dG[0]),cimag(dG[0]));
  printf("dqG/dy=% 15.14e %+15.14e I\n",creal(dG[1]),cimag(dG[1]));
  printf("dqG/dz=% 15.14e %+15.14e I\n",creal(dG[2]),cimag(dG[2]));
}

