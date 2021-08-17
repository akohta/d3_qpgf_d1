// quasi-periodic Green's function in cylindrical coordinate system
#include "d3_qpgf_d1_cs.h"

// verification of second derivative using central difference
void central_diff(double complex *d2qG,double *cr,QPDT *qd);

int main()
{
  QPDT qd;
  double lambda,k,theta,phi,kx,r[3],cr[3],d,eps;
  double complex qG,dqG[3],d2qG[3],fv;
  int err;

  eps=1.0e-15; // relative error requested
  lambda=1.0;  // wave length
  theta=0.3;   // wave vector parameter
  phi=0.5;     // kx=k*sin(theta)*cos(phi), ky=k*sin(theta)*sin(phi), kz=k*cos(theta)
  d=0.5;       // periodic number, lattice vector is (d,0,0)

  k=2.0*M_PI/lambda;
  kx=k*sin(theta)*cos(phi);
  qd.k=k;
  qd.d=d;
  qd.kx=kx;
  printf("eps=%g\n",eps);
  printf("k  = %g\n",k);
  printf("k_x= %g\n",kx);
  printf("d  = %g\n\n",d);
  
  printf("--- first derivative of quasi-periodic Green's function ---\n"); 
  r[0]=-2.1;
  r[1]=0.3;
  r[2]=-0.35;
  cr[0]=r[0];
  cr[1]=sqrt(r[1]*r[1]+r[2]*r[2]);
  cr[2]=atan2(r[2],r[1]);
  printf("(x,y,z) = (% 5g, % 5g,% 5g), (x,rho,phi) = (% 5g,% 5g,% 5g)\n",r[0],r[1],r[2],cr[0],cr[1],cr[2]);
  err=d3hm_qpgf_d1_cv1(&qG,dqG,cr,eps,&qd);
  printf("qG            =% 15.14e %+15.14e I,err=%d\n",creal(qG),cimag(qG),err);
  printf("dqG/dx        =% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("1/rho*dqG/drho=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dphi      =% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));
  fv=dqG[1]*r[1]; // dqG/dy=(dqG/drho)*(drho/dy)=(1/rho*dqG/drho)*y
  printf("dqG/dy        =% 15.14e %+15.14e I (verification)\n",creal(fv),cimag(fv));
  fv=dqG[1]*r[2]; // dqG/dz=(dqG/drho)*(drho/dz)=(1/rho*dqG/drho)*z
  printf("dqG/dz        =% 15.14e %+15.14e I (verification)\n",creal(fv),cimag(fv)); 

  r[0]=0.1;
  r[1]=0.04;
  r[2]=0.01;
  cr[0]=r[0];
  cr[1]=sqrt(r[1]*r[1]+r[2]*r[2]);
  cr[2]=atan2(r[2],r[1]);
  printf("(x,y,z) = (% 5g, % 5g,% 5g), (x,rho,phi) = (% 5g,% 5g,% 5g)\n",r[0],r[1],r[2],cr[0],cr[1],cr[2]);
  err=d3hm_qpgf_d1_cv1(&qG,dqG,cr,eps,&qd);
  printf("qG            =% 15.14e %+15.14e I,err=%d\n",creal(qG),cimag(qG),err);
  printf("dqG/dx        =% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("1/rho*dqG/drho=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dphi      =% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2]));
  fv=dqG[1]*r[1]; // dqG/dy=(dqG/drho)*(drho/dy)=(1/rho*dqG/drho)*y
  printf("dqG/dy        =% 15.14e %+15.14e I (verification)\n",creal(fv),cimag(fv));
  fv=dqG[1]*r[2]; // dqG/dz=(dqG/drho)*(drho/dz)=(1/rho*dqG/drho)*z
  printf("dqG/dz        =% 15.14e %+15.14e I (verification)\n",creal(fv),cimag(fv)); 


  printf("\n--- second derivative of quasi-periodic Green's function ---\n"); 
  r[0]=-2.1;
  r[1]=0.3;
  r[2]=-0.35;
  cr[0]=r[0];
  cr[1]=sqrt(r[1]*r[1]+r[2]*r[2]);
  cr[2]=atan2(r[2],r[1]);
  printf("(x,y,z) = (% 5g, % 5g,% 5g), (x,rho,phi) = (% 5g,% 5g,% 5g)\n",r[0],r[1],r[2],cr[0],cr[1],cr[2]);
  err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,cr,eps,&qd);
  printf("qG            =% 15.14e %+15.14e I,err=%d\n",creal(qG),cimag(qG),err);
  printf("dqG/dx        =% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("1/rho*dqG/drho=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dphi      =% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2])); 
  printf("d^2qG/dx^2    =% 15.14e %+15.14e I\n",creal(d2qG[0]),cimag(d2qG[0]));
  printf("d^2qG/dphidx  =% 15.14e %+15.14e I\n",creal(d2qG[1]),cimag(d2qG[1]));
  printf("d^2qG/dphi^2  =% 15.14e %+15.14e I\n",creal(d2qG[2]),cimag(d2qG[2]));
  //verification
  central_diff(d2qG,cr,&qd);
  
  r[0]=0.1;
  r[1]=0.04;
  r[2]=0.01;
  cr[0]=r[0];
  cr[1]=sqrt(r[1]*r[1]+r[2]*r[2]);
  cr[2]=atan2(r[2],r[1]);
  printf("(x,y,z) = (% 5g, % 5g,% 5g), (x,rho,phi) = (% 5g,% 5g,% 5g)\n",r[0],r[1],r[2],cr[0],cr[1],cr[2]);
  err=d3hm_qpgf_d1_cv2(&qG,dqG,d2qG,cr,eps,&qd);
  printf("qG            =% 15.14e %+15.14e I,err=%d\n",creal(qG),cimag(qG),err);
  printf("dqG/dx        =% 15.14e %+15.14e I\n",creal(dqG[0]),cimag(dqG[0]));
  printf("1/rho*dqG/drho=% 15.14e %+15.14e I\n",creal(dqG[1]),cimag(dqG[1]));
  printf("dqG/dphi      =% 15.14e %+15.14e I\n",creal(dqG[2]),cimag(dqG[2])); 
  printf("d^2qG/dx^2    =% 15.14e %+15.14e I\n",creal(d2qG[0]),cimag(d2qG[0]));
  printf("d^2qG/dphidx  =% 15.14e %+15.14e I\n",creal(d2qG[1]),cimag(d2qG[1]));
  printf("d^2qG/dphi^2  =% 15.14e %+15.14e I\n",creal(d2qG[2]),cimag(d2qG[2]));
  //verification
  central_diff(d2qG,cr,&qd);
  return 0;
}

void central_diff(double complex *d2qG,double *cr,QPDT *qd)
{
  double complex fv,dqGp[3],dqGm[3];
  double h,ct[3];
  
  h=qd->d*1.0e-7;
  
  printf("verification of second derivative using central difference\n");
  
  ct[0]=cr[0]+h;
  ct[1]=cr[1];
  ct[2]=ct[2];
  d3hm_qpgf_d1_cv1_dqG(dqGp,ct,1.0e-15,qd);
  dqGp[1]*=ct[1];
  ct[0]=cr[0]-h;
  d3hm_qpgf_d1_cv1_dqG(dqGm,ct,1.0e-15,qd);
  dqGm[1]*=ct[1];
  fv=(dqGp[0]-dqGm[0])/(2.0*h);
  printf("d^2qG/dx^2    =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
  fv=(dqGp[1]-dqGm[1])/(2.0*h);
  printf("d^2qG/dphidx  =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));

  ct[0]=cr[0];
  ct[1]=cr[1]+h;
  ct[2]=ct[2];
  d3hm_qpgf_d1_cv1_dqG(dqGp,ct,1.0e-15,qd);
  dqGp[1]*=ct[1];
  ct[1]=cr[1]-h;
  d3hm_qpgf_d1_cv1_dqG(dqGm,ct,1.0e-15,qd);
  dqGm[1]*=ct[1];
  fv=(dqGp[1]-dqGm[1])/(2.0*h);
  printf("d^2qG/dphi^2  =% 15.14e %+15.14e I\n",creal(fv),cimag(fv));
}

