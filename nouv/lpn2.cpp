#include <cstdio>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

const int N=256;
const double p=0.125;
const int M=1<<10;

double logBinomSecret[N+1];
double logBinomExemple[M+1];
double nbSecret[N+1];
double nbGaussienne[M+1];

void calcLog(double *t,int N,double p){
	t[0]=N*log(1-p);
	for (int i=1;i<=N;i++)
		t[i]=t[i-1]+log(p/(1-p));
}

double log1M(double a){
	if (a>1e-3) return 1-a<1e-8 ? -50 : log(1-a);
	return -a-a*a/2-a*a*a/3;
}

double S(int N,int p){
	double nb=1,total=0;
	for (int poids=0;poids<=p;poids++){
		total+=nb;
		nb*=(N-poids)/(poids+1.);
	}
	return total;
}

double binom(int n,int p){
	if (p<0 || p>n)
		return 0;
	double r=1;
	for (int i=1;i<=p;i++)
		r=(r*(n-i+1.))/i;
	return r;
}

double logAdd(double a,double b){
	// log(exp(a)+exp(b))
	double maxi=max(a,b),mini=min(a,b);
	return maxi+log1M(-exp(max(-50.,mini-maxi)));
}

double logb2(double x){
	return log(x)/log(2);
}

double probaWalsh(int N,int P,int M,double bruit){
	calcLog(logBinomSecret,N,p);
	calcLog(logBinomExemple,M,bruit);
	double nbBinom=-M*log(2);
	for (int i=0;i<=M;i++){
		nbGaussienne[i]=i>0 ? logAdd(nbGaussienne[i-1],nbBinom) : nbBinom;
		if (i<M) nbBinom+=log((M-i)/(i+1.));
	}
	double proba=0;
	double pSecret=0;
	double probaMaxi=0;
	for (int p=0;p<=P;p++)
		nbSecret[p]=binom(N,p);
	for (int poids=0;poids<=P;poids++){
		const int q=-(M*log(1-bruit)-logBinomSecret[0]-logBinomExemple[int(bruit*M)]+logBinomSecret[poids])/log(bruit/(1-bruit));
		if (nbGaussienne[q]>-2+log(nbSecret[poids])) // on s'assurre d'être capable de distinguer le secret
			return 1e-100;
	}
	for (int nE=0;nE<=M;nE++){
		for (int poids=0;poids<=P;poids++){
			const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
			double cur=vrais+pSecret+log(nbSecret[poids]);
			if (cur>0) printf("%d(%d) %d %f : %f+%f  %f\n",nE,M,poids,cur,vrais,pSecret,logBinomExemple[nE]),exit(0);
			probaMaxi+=exp(max(-50.,cur));
			for (int p2=0;p2<=P;p2++){
				const int q=-(M*log(1-bruit)-vrais+logBinomSecret[p2])/log(bruit/(1-bruit));
				const double pAuDessus=exp(nbGaussienne[max(0,min(M,q))]);
				cur+=log1M(pAuDessus)*(nbSecret[p2]-(poids==p2));
//					printf("%d %d %f %f %f\n",nE,poids,vrais,cur,pAuDessus);
			}
			proba+=exp(max(-50.,cur));
		}
		if (nE<M) pSecret+=log((M-nE)/(nE+1.));
	}
//	printf("M=%d b=%f => p=%f pMaxi=%f\n",M,bruit,proba,probaMaxi);
	return proba/probaMaxi;
}

int main(){
	srand(43);
	double ancTravail=pow(2,N);
	char param[100]="aucun\n";
	int P3=0;
	double pCorrect=1;
	for (P3=1;2*pCorrect-1>0.5;P3++){
		pCorrect=0;
		for (int i=0;i<=P3;i+=2)
			pCorrect+=binom(P3,i)*pow(p,i)*pow(1-p,P3-i);
		printf("%d %f\n",P3,2*pCorrect-1);
	}
	P3--;
	for (int m=10;;m++){
		double ancTravail2=pow(2,N);
		double minTravail=ancTravail2;
		for (int w=0;;w++){
			const int n=N-m*w;
			const double tau=pow(1-2*p,1<<w);
			if (tau*pow(2,m)*tau<0.01)
				break;
			double ancTravail3=pow(2,N);
			double nbOp=3*N*N+pow(2,m)*(w+1)*N;
			for (int n2=min(n,max(0,n-2*m));n2<=n;n2++){ /* facteur 2 à régler… On s'attend à 1/p+o(1) !*/
				const int L=n==n2 ? m : m-n+n2+logb2(S(n-n2,min(n-n2,P3)));
//				printf("%d %d\n",n2,L);
				if (pow(2,L)*tau*tau<1)
					continue;
				double ancTravail4=pow(2,N);
				double probaSecretCur=L*log(1-p),probaSecret=0;
				int P1=0;
				for (P1=0;probaSecret<0.9;P1++){
					probaSecret+=exp(probaSecretCur);
					if (P1<L) probaSecretCur+=log(p/(1-p)*(L-P1)/(P1+1.));
				}
				double nbOp2=0;
				for (int l=0;l<=L;l++)
					nbOp2+=pow(2,l)*S(L-l,P1)*(L-l+1);
				double nbSecret=1;
				double probaSecret2=0;
				for (int P=0;P<=n2-L;P++){
//					printf("P2=%d\n",P);
					probaSecret2+=nbSecret*pow(p,P)*pow(1-p,n2-L-P);
					const double probaDetect=M>pow(2,L) ? probaWalsh(n2,P,1<<L,(1-tau)/2) : probaWalsh(n2,P,M,max(1e-5,(1-tau*sqrt(pow(2,L)/M))/2));
					const double travail=(nbOp+nbOp2*nbSecret)/probaSecret/probaSecret2/probaDetect;
//					printf("P2=%d %g %f %g\n",P,probaSecret2,logb2(travail),probaDetect);
					if (travail>ancTravail4*4)
						break;
					ancTravail4=min(ancTravail4,travail);
					nbSecret*=(n2-L-P)/(P+1.);
					if (travail<minTravail){
						minTravail=travail;
						sprintf(param,"m=%d w=%d n2=%d (L=%d : Ltau²=%f) P1=%d P2=%d (%d perdus): %f (%f)\n",m,w,n2,L,pow(2,L)*tau*tau,P1,P,n-n2,logb2(travail),probaDetect);
					}
				}
//				printf("m=%d w=%d n2=%d (L=%d) P1=%d : %f\n",m,w,n2,L,P1,logb2(ancTravail4));
				if (ancTravail4>ancTravail3*64)
					break;
				ancTravail3=min(ancTravail3,ancTravail4);
			}
			if (ancTravail3>ancTravail2*8)
				break;
			ancTravail2=min(ancTravail2,ancTravail3);
		}
		printf("%s",param);
		if (ancTravail2>ancTravail*4)
			break;
		ancTravail=min(ancTravail,ancTravail2);
	}
	return 0;
}
