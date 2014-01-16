#include <cstdio>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

const int N=50,L=20,M=1<<16; // M|2^L
const double p=0.2,bruit=0.497;

double logBinomSecret[N+1];
double logBinomExemple[M+1];
pair<double,double> vraisTotale[2*N+1][M+2];
double nbSecret[N+1];
double nbGaussienne[M+1];

void calcLog(double *t,int N,double p,int pas=1){
	t[0]=pas*N*log(1-p);
	for (int i=1;i<=N;i++)
		t[i]=t[i-1]+pas*log(p/(1-p));
}

double log1M(double a){ // log(1-a)
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

double H(double x){
	if (x<1e-6 || 1-x<1e-6)
		return 0;
	return x*log(x)+(1-x)*log(1-x);
}

double logAdd(double a,double b){
	// log(exp(a)+exp(b))
	double maxi=max(a,b),mini=min(a,b);
	return maxi+log1M(-exp(max(-50.,mini-maxi)));
}

int main(){
	srand(43);
	calcLog(logBinomSecret,N,p);
	calcLog(logBinomExemple,M,bruit,(1<<L)/M);
	for (int P=0;P<=N/2;P++){
		double travail=0;
		for (int l=0;l<=L;l++)
			travail+=S(N-l,P)*pow(2,l);
		nbSecret[P]=binom(N,P);
		double proba=0;
		double pSecret=0;
		double probaMaxi=0;
		for (int nE=0;nE<=M;nE++){
			for (int poids=0;poids<=P;poids++){
				const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
				double cur=logBinomSecret[poids]+logBinomExemple[nE]*(M/pow(2,L))+pSecret+log(binom(N,poids));
				probaMaxi+=exp(max(-50.,cur));
				for (int p2=0;p2<=P;p2++){
					const double q=-(pow(2,L)*log(1-bruit)-vrais+logBinomSecret[p2])/log(bruit/(1-bruit));
					const double pAuDessus=q>=pow(2,L-1) ? 1 : pow(2,pow(2,L)*(H(ceil(q)/pow(2,L))-1));
//					const double pAuDessus=exp(nbGaussienne[min(M,int((q*M)/pow(2,L)))]*pow(2,L)/M);
					cur+=log1M(pAuDessus)*(nbSecret[p2]-(poids==p2));
//					printf("%d %d %f %f %f (%f)\n",nE,poids,vrais,cur,pAuDessus,q);
				}
				if (cur>50) // anti-Heisenbug
					return 0;
				proba+=exp(max(-50.,cur));
			}
			if (nE<M) pSecret+=log((M-nE)/(nE+1.));
		}
		printf("%g %g %g\n",proba,travail/proba,probaMaxi);
//		printf("%d\n",P);
//		return 0;
	}
	return 0;
	puts("");
	puts("");
	for (int P=0;P<=N-L;P++){
		double travail=pow(2,L)*L*S(N-L,P);
		for (int poids=0;poids<=P+L;poids++){
			double curP=0;
			for (int p1=max(poids-P,0);p1<=L;p1++)
				curP+=binom(L,p1)*binom(N-L,poids-p1);
			double curLogP=log(curP)-M*log(2);
			for (int nE=0;nE<=M;nE++){
				const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
				vraisTotale[poids][nE]=make_pair(vrais,exp(curLogP));
//				printf("%f %f\n",vrais,exp(curLogP));
				if (nE<M) curLogP+=log((M-nE)/(nE+1.));
			}
//			puts("");
			sort(vraisTotale[poids],vraisTotale[poids]+M+1);
			for (int i=M-1;i>=0;i--)
				vraisTotale[poids][i].second+=vraisTotale[poids][i+1].second;
/*			for (int i=0;i<M+1;i++)
				printf("%f %f\n",vraisTotale[poids][i].first,vraisTotale[poids][i].second);
			puts("");*/
		}
		double proba=0;
		double pSecret=0;
		double probaMaxi=0;
		for (int nE=0;nE<=M;nE++){
			for (int poids=0;poids<=P+L;poids++){
				const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
				double cur=log(vraisTotale[poids][0].second)+poids*log(p)+(N-poids)*log(1-p)+logBinomExemple[nE]*M/pow(2,L)+pSecret;
				probaMaxi+=exp(max(-50.,cur));
				for (int p2=0;p2<=P+L;p2++){
					const double pAuDessus=upper_bound(vraisTotale[p2],vraisTotale[p2]+M+1,make_pair(vrais,0.))->second/vraisTotale[p2][0].second;
					cur+=log1M(pAuDessus)*(vraisTotale[p2][0].second-(poids==p2));
				}
//				printf("%d %d %f %f\n",nE,poids,vrais,cur);
				proba+=exp(max(-50.,cur));
			}
			if (nE<M) pSecret+=log((M-nE)/(nE+1.));
		}
		printf("%g %g %g\n",proba,travail/proba,proba/probaMaxi);
//		return 0;
	}
	return 0;
}
