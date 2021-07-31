#include "d3_qpgf_d1.h"


int d3hm_qpgf_d1_fd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd)
{
	int d3hm_qpgf_d1_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd);

	double complex c;
	double rt[3],tmp;
	int err,l;

	rt[0]=r[0];
	rt[1]=r[1];
	rt[2]=r[2];

	l=(int)floor(r[0]/qd->d+0.5);
	rt[0]-=(double)l*qd->d;

	err=d3hm_qpgf_d1_fd_sd(qG,dqG,rt,eps,qd);
	if(l!=0){
		tmp=(double)l*qd->kx*qd->d;
		c=cos(tmp)-I*sin(tmp);
		*qG*=c;
		dqG[0]*=c;
		dqG[1]*=c;
		dqG[2]*=c;
	}

	return err;
}

int d3hm_qpgf_d1_ew(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd)
{
	int d3hm_qpgf_d1_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd);

	double complex c;
	double rt[3],tmp;
	int err,l;

	rt[0]=r[0];
	rt[1]=r[1];
	rt[2]=r[2];

	l=(int)floor(r[0]/qd->d+0.5);
	rt[0]-=(double)l*qd->d;

	err=d3hm_qpgf_d1_ew_sd(qG,dqG,rt,eps,qd);
	if(l!=0){
		tmp=(double)l*qd->kx*qd->d;
		c=cos(tmp)-I*sin(tmp);
		*qG*=c;
		dqG[0]*=c;
		dqG[1]*=c;
		dqG[2]*=c;
	}

	return err;
}

int d3hm_qpgf_d1(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd)
{

	double kd,k2,bN2,aN,ap,bp2,am,bm2,FN,Fp,Fm,ro,abn;
	int N,err;

	N=(int)round(qd->kx*qd->d/(2.0*M_PI));

	kd=2.0*M_PI/qd->d;
	k2=qd->k*qd->k;

	ap=kd*(double)(N+S_LIMIT)-qd->kx;
	bp2=k2-ap*ap;
	am=kd*(double)(N-S_LIMIT)-qd->kx;
	bm2=k2-am*am;

	if(bp2>0.0 || bm2>0.0){ // Ewald method
		err=d3hm_qpgf_d1_ew(qG,dqG,r,eps,qd);
		if(err<0) err-=20;
	}

	ro=sqrt(r[1]*r[1]+r[2]*r[2]);
	aN=kd*(double)N-qd->kx;
	bN2=k2-aN*aN;

	abn=sqrt(fabs(bN2));
	FN=0.5*M_PI*cabs(besj0(abn*ro)+I*besy0(abn*ro));
	abn=sqrt(fabs(bp2));
	Fp=fabs(bessk0(abn*ro));
	abn=sqrt(fabs(bm2));
	Fm=fabs(bessk0(abn*ro));

	if(Fp < eps*FN && Fm < eps*FN){ // Fourier
		err=d3hm_qpgf_d1_fd(qG,dqG,r,eps,qd);
		if(err<0) err-=10;
	}
	else { // Ewald
		err=d3hm_qpgf_d1_ew(qG,dqG,r,eps,qd);
		if(err<0) err-=20;
	}

	return err;
}


/////////////////////////////////////////////////////////////////////////////////////
int d3hm_qpgf_d1_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd)
{

	double complex ekx,ten,i_ten,eiaxN,eiaxp,eiaxm,ch,tc,tf,oG,odG0,odG1;
	double arg,ro,k2,kd,kdx,bn2,aN,anp,anm,bn,bnr,ct,st;
	int i,N,flg[3],n,sig;

	ro=sqrt(r[1]*r[1]+r[2]*r[2]);
	k2=qd->k*qd->k;

	N=(int)round(qd->kx*qd->d/(2.0*M_PI));
	ch=M_PI*I/2.0;

	arg=-qd->kx*r[0];
	ekx=cos(arg)+I*sin(arg);

	kd=2.0*M_PI/qd->d;
	kdx=kd*r[0];
	ct=cos(kdx);
	st=sin(kdx);
	ten=ct+I*st;
	i_ten=ct-I*st;

	for(i=0;i<3;i++) flg[i]=0.0;

	n=N;
	aN=kd*(double)n-qd->kx;
	arg=kdx*(double)n;
	eiaxN=(cos(arg)+I*sin(arg))*ekx;
	bn2=k2-aN*aN;
	if(bn2>0.0){
		bn=sqrt(bn2);
		bnr=bn*ro;
		tc=eiaxN*ch;
		tf=tc*(besj0(bnr)+I*besy0(bnr));
		*qG=tf;
		dqG[0]=aN*tf;
		dqG[1]=-bn*tc*(besj1(bnr)+I*besy1(bnr));
	}
	else if(bn2<0.0){
		bn=sqrt(-bn2);
		bnr=bn*ro;
		tf=eiaxN*bessk0(bnr);
		*qG=tf;
		dqG[0]=aN*tf;
		dqG[1]=-bn*eiaxN*bessk1(bnr);
	}
	else return -2; // beta_n=0.0;

	eiaxp=eiaxN;
	eiaxm=eiaxN;
	anp=aN;
	anm=aN;
	for(n=0;n<FD_LIMIT;n++){
		oG=*qG;
		odG0=dqG[0];
		odG1=dqG[1];
		sig=0;

		eiaxp*=ten;
		eiaxm*=i_ten;
		anp+=kd;
		anm-=kd;

		// +
		bn2=k2-anp*anp;
		if(bn2>0.0){
			bn=sqrt(bn2);
			bnr=bn*ro;
			tc=eiaxp*ch;
			tf=tc*(besj0(bnr)+I*besy0(bnr));
			*qG+=tf;
			dqG[0]+=anp*tf;
			dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
			sig+=1;
		}
		else if(bn2<0.0){
			bn=sqrt(-bn2);
			bnr=bn*ro;
			tf=eiaxp*bessk0(bnr);
			*qG+=tf;
			dqG[0]+=anp*tf;
			dqG[1]+=-bn*eiaxp*bessk1(bnr);
		}
		else return -2; // beta_n=0.0;

		// -
		bn2=k2-anm*anm;
		if(bn2>0.0){
			bn=sqrt(bn2);
			bnr=bn*ro;
			tc=eiaxm*ch;
			tf=tc*(besj0(bnr)+I*besy0(bnr));
			*qG+=tf;
			dqG[0]+=anm*tf;
			dqG[1]+=-bn*tc*(besj1(bnr)+I*besy1(bnr));
			sig+=1;
		}
		else if(bn2<0.0){
			bn=sqrt(-bn2);
			bnr=bn*ro;
			tf=eiaxm*bessk0(bnr);
			*qG+=tf;
			dqG[0]+=anm*tf;
			dqG[1]+=-bn*eiaxm*bessk1(bnr);
		}
		else return -2; // beta_n=0.0;

		if(sig==0){
			if(flg[0]==0){
				if(  2.0*fabs(creal(oG)-creal(*qG)) < eps*(fabs(creal(oG))+fabs(creal(*qG)))
					&& 2.0*fabs(cimag(oG)-cimag(*qG)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG))) ) flg[0]=1;
			}
			if(flg[1]==0){
				if(  2.0*fabs(creal(odG0)-creal(dqG[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG[0])))
					&& 2.0*fabs(cimag(odG0)-cimag(dqG[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG[0]))) ) flg[1]=1;
			}
			if(flg[2]==0){
				if(  2.0*fabs(creal(odG1)-creal(dqG[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG[1])))
					&& 2.0*fabs(cimag(odG1)-cimag(dqG[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG[1]))) ) flg[2]=1;
			}
		}

		if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
	}

	if(n==FD_LIMIT) return -1;

	tc=1.0/(2.0*M_PI*qd->d);
	*qG*=tc;
	dqG[0]*=I*tc;

	dqG[2]=dqG[1];
	tc/=ro;
	dqG[1]*=r[1]*tc;
	dqG[2]*=r[2]*tc;

	return 2*n+1;
}


int d3hm_qpgf_d1_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPDT *qd)
{
	int d3hm_ewald_qG1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPDT *qd);
	int d3hm_ewald_qG2(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPDT *qd);

	double complex tqG,tdqG[3];
	double veps,r02;
	int err1,err2;

	veps=qd->d/(2.0*sqrt(M_PI));

	// veps
	r02=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
	if( (-r02/(4.0*veps*veps)+qd->k*qd->k*veps*veps) > EW_EXPM)
		veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r02+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

	// qG1
	err1=d3hm_ewald_qG1(&tqG,tdqG,r,eps,veps,qd);
	*qG=tqG;
	dqG[0]=tdqG[0];
	dqG[1]=tdqG[1];
	dqG[2]=tdqG[2];
	if(err1==-1) return -1;

	// qG2
	err2=d3hm_ewald_qG2(&tqG,tdqG,r,eps,veps,qd);
	*qG+=tqG;
	dqG[0]+=tdqG[0];
	dqG[1]+=tdqG[1];
	dqG[2]+=tdqG[2];
	if(err2==-1) return -2;
	if(err2==-2) return -3; // beta_n==0.0

	return err1+err2;
}

int d3hm_ewald_qG1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPDT *qd)
{
	double complex ekd,i_ekd,eLkd,tw,elkdp,elkdm,tce,tcf,oG,odG0,odG1;
	double kd,ro2,rl2,rl,i_rl,xLd,arg,ke,ke2,i2e,mi2e2,xldp,xldm,st,ct,td,cpe;
	int L,l,flg[3];

	kd=qd->kx*qd->d;
	st=sin(kd);
	ct=cos(kd);
	ekd=ct+I*st;
	i_ekd=ct-I*st;

	ro2=r[1]*r[1]+r[2]*r[2];
	ke=qd->k*veps;
	ke2=ke*ke;
	i2e=1.0/(2.0*veps);
	mi2e2=-i2e*i2e;
	cpe=1.0/(sqrt(M_PI)*veps);

	if(fabs(r[0])<qd->d*sqrt(eps)) L=1; // avoid symmetric cancellation
	else L=0;
	flg[0]=0;	flg[1]=0;	flg[2]=0;

	arg=kd*(double)L;
	eLkd=cos(arg)+I*sin(arg);
	xLd=r[0]+(double)L*qd->d;
	rl2=xLd*xLd+ro2;
	rl=sqrt(rl2);
	i_rl=1.0/rl;
	tw=Faddeeva_w(ke+I*rl*i2e,eps);
	tce=eLkd*i_rl*exp(rl2*mi2e2+ke2);
	tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
	*qG1=tce*creal(tw);
	dqG1[0]=xLd*tcf;
	dqG1[1]=tcf;

	elkdp=eLkd;
	elkdm=eLkd;
	xldp=xLd;
	xldm=xLd;
	for(l=0;l<EW_LIMIT;l++){
		oG=*qG1;
		odG0=dqG1[0];
		odG1=dqG1[1];
		elkdp*=ekd;
		elkdm*=i_ekd;
		xldp+=qd->d;
		xldm-=qd->d;

		// +
		rl2=xldp*xldp+ro2;
		rl=sqrt(rl2);
		i_rl=1.0/rl;
		tw=Faddeeva_w(ke+I*rl*i2e,eps);
		tce=elkdp*i_rl*exp(rl2*mi2e2+ke2);
		tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
		*qG1+=tce*creal(tw);
		dqG1[0]+=xldp*tcf;
		dqG1[1]+=tcf;

		//-
		rl2=xldm*xldm+ro2;
		rl=sqrt(rl2);
		i_rl=1.0/rl;
		tw=Faddeeva_w(ke+I*rl*i2e,eps);
		tce=elkdm*i_rl*exp(rl2*mi2e2+ke2);
		tcf=i_rl*tce*(i_rl*creal(tw)-qd->k*cimag(tw)+cpe);
		*qG1+=elkdm/rl*exp(rl2*mi2e2+ke2)*creal(tw);
		dqG1[0]+=xldm*tcf;
		dqG1[1]+=tcf;

		if(flg[0]==0){
			if(  2.0*fabs(creal(oG)-creal(*qG1)) < eps*(fabs(creal(oG))+fabs(creal(*qG1)))
				&& 2.0*fabs(cimag(oG)-cimag(*qG1)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG1))) ) flg[0]=1;
		}
		if(flg[1]==0){
			if(  2.0*fabs(creal(odG0)-creal(dqG1[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG1[0])))
				&& 2.0*fabs(cimag(odG0)-cimag(dqG1[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG1[0]))) ) flg[1]=1;
		}
		if(flg[2]==0){
			if(  2.0*fabs(creal(odG1)-creal(dqG1[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG1[1])))
				&& 2.0*fabs(cimag(odG1)-cimag(dqG1[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG1[1]))) ) flg[2]=1;
		}

		if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
	}

	td=1.0/(4.0*M_PI);
	*qG1*=td;
	dqG1[0]*=-td;
	dqG1[2]=dqG1[1];
	dqG1[1]*=-td*r[1];
	dqG1[2]*=-td*r[2];

	if(l!=EW_LIMIT) return 2*l+1;
	else return -1;
}

int d3hm_ewald_qG2(double complex *qG2,double complex *dqG2,double *r,double eps,double veps,QPDT *qd)
{
	int sum_expint(double complex *ret0,double complex *ret1,double ro2,double veps,double bn2,double eps);

	double complex ten,i_ten,eiaxN,tc0,tc1,eiaxp,eiaxm,oG,odG0,odG1;
	double ro2,kd,k2,ve2,arg,ct,st,aN,bn2,anp,anm;
	int N,n,err,flg[3],sig;

	flg[0]=0;	flg[1]=0;	flg[2]=0;
	N=(int)round(qd->kx*qd->d/(2.0*M_PI));

	ro2=r[1]*r[1]+r[2]*r[2];
	kd=2.0*M_PI/qd->d;
	k2=qd->k*qd->k;
	ve2=veps*veps;

	arg=kd*r[0];
	st=sin(arg);	ct=cos(arg);
	ten=ct+I*st;
	i_ten=ct-I*st;

	n=N;
	aN=(double)n*kd-qd->kx;
	arg=aN*r[0];
	eiaxN=cos(arg)+I*sin(arg);
	bn2=k2-aN*aN; if(bn2==0.0) return -2;
	err=sum_expint(&tc0,&tc1,ro2,veps,bn2,eps);	if(err<0) return -1;
	*qG2=eiaxN*tc0;
	dqG2[0]=aN*eiaxN*tc0;
	dqG2[1]=eiaxN*tc1;

	eiaxp=eiaxN;
	eiaxm=eiaxN;
	anp=aN;
	anm=aN;
	for(n=0;n<EW_LIMIT;n++){
		oG=*qG2;
		odG0=dqG2[0];
		odG1=dqG2[1];
		eiaxp*=ten;
		eiaxm*=i_ten;
		anp+=kd;
		anm-=kd;
		sig=0;

		// +
		bn2=k2-anp*anp;	if(bn2==0.0) return -2;
		if(bn2<0.0) sig++;
		err=sum_expint(&tc0,&tc1,ro2,veps,bn2,eps);	if(err<0) return -1;
		*qG2+=eiaxp*tc0;
		dqG2[0]+=anp*eiaxp*tc0;
		dqG2[1]+=eiaxp*tc1;

		// -
		bn2=k2-anm*anm;	if(bn2==0.0) return -2;
		if(bn2<0.0) sig++;
		err=sum_expint(&tc0,&tc1,ro2,veps,bn2,eps);	if(err<0) return -1;
		*qG2+=eiaxm*tc0;
		dqG2[0]+=anm*eiaxm*tc0;
		dqG2[1]+=eiaxm*tc1;

		if(sig>0){
			if(flg[0]==0){
				if(  2.0*fabs(creal(oG)-creal(*qG2)) < eps*(fabs(creal(oG))+fabs(creal(*qG2)))
					&& 2.0*fabs(cimag(oG)-cimag(*qG2)) < eps*(fabs(cimag(oG))+fabs(cimag(*qG2))) ) flg[0]=1;
			}
			if(flg[1]==0){
				if(  2.0*fabs(creal(odG0)-creal(dqG2[0])) < eps*(fabs(creal(odG0))+fabs(creal(dqG2[0])))
					&& 2.0*fabs(cimag(odG0)-cimag(dqG2[0])) < eps*(fabs(cimag(odG0))+fabs(cimag(dqG2[0]))) ) flg[1]=1;
			}
			if(flg[2]==0){
				if(  2.0*fabs(creal(odG1)-creal(dqG2[1])) < eps*(fabs(creal(odG1))+fabs(creal(dqG2[1])))
					&& 2.0*fabs(cimag(odG1)-cimag(dqG2[1])) < eps*(fabs(cimag(odG1))+fabs(cimag(dqG2[1]))) ) flg[2]=1;
			}
		}

		if(flg[0]!=0 && flg[1]!=0 && flg[2]!=0) break;
	}

	tc0=1.0/(4.0*M_PI*qd->d);
	tc1=tc0/(2.0*ve2);
	*qG2*=tc0;
	dqG2[0]*=I*tc0;
	dqG2[2]=dqG2[1];
	dqG2[1]*=r[1]*tc1;
	dqG2[2]*=r[2]*tc1;

	if(n!=EW_LIMIT) return 2*n+1;
	else return -1;
}

int sum_expint(double complex *ret0,double complex *ret1,double ro2,double veps,double bn2,double eps)
{
	double complex E1,E0,tG0,tG0o,tG1,tG1o;
	double sig,ipm,rpm,z,im,dp,rpm1,ipm1;
	int m;

	dp=ro2/(4.0*veps*veps);
	z=-bn2*veps*veps;

	sig=1.0;
	ipm=1.0;
	rpm=1.0;
	if(bn2>0.0) E1=I*M_PI-expint_i(-z);
	else if(bn2<0.0) E1=expint(1,z);
	else return -3;
	tG0=sig*ipm*rpm*E1;
	tG1=0.0;
	for(m=1;m<EW_LIMIT;m++){
		tG0o=tG0;
		tG1o=tG1;
		sig*=-1.0;
		im=1.0/(double)m;
		ipm1=ipm;
		ipm*=im;
		rpm1=rpm;
		rpm*=dp;
		E0=E1;
		E1=im*(exp(-z)-z*E0);
		tG0+=sig*ipm*rpm*E1;
		tG1+=sig*ipm1*rpm1*E1;
		if(2.0*cabs(tG0o-tG0)/(cabs(tG0o)+cabs(tG0))<eps) break;
	}

	*ret0=tG0;
	*ret1=tG1;

	if(m==EW_LIMIT) return -1;
	else return 2*m-1;
}
