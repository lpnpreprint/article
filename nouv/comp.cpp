#include <cstdio>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

const int N=391,L=14,M=1<<14; // M|2^L
const double p=0.125,bruit=0.494988702;

double logBinomSecret[N+1];
double logBinomExemple[M+1];
pair<double,double> vraisTotale[2*N+1][M+2];
double tailleListe[N+1];

void calcLog(double *t,int N,double p,int pas=1){
	t[0]=N*log(1-p);
	for (int i=1;i<=N;i++)
		t[i]=t[i-1]+pas*log(p/(1-p));
}

double logPuiss(double a,double b){
	if (a>1e-3) return 1-a<1e-8 ? -50 : b*log(1-a);
	return b*(-a-a*a/2-a*a*a/3);
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

int main(){
	srand(43);
	calcLog(logBinomSecret,N,p);
	calcLog(logBinomExemple,M,bruit,(1<<L)/M);
	double curProba=1;
	for (int poids=0;poids<=N;poids++){
		double curLogP=log(curProba)-M*log(2);
		for (int nE=0;nE<=M;nE++){
			const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
			vraisTotale[poids][nE]=make_pair(vrais,exp(curLogP));
			if (nE<M) curLogP+=log((M-nE)/(nE+1.));
		}
		curProba*=(N-poids)/(poids+1.);
		sort(vraisTotale[poids],vraisTotale[poids]+M+1);
		for (int i=M-1;i>=0;i--)
			vraisTotale[poids][i].second+=vraisTotale[poids][i+1].second;
	}
	for (int P=0;P<=N/20;P++){
		double travail=0;
		for (int l=0;l<=L;l++){
			tailleListe[l]=S(N-l,P);
			travail+=S(N-l,P)*pow(2,l);
		}
		double proba=0;
		double pSecret=0;
		double probaMaxi=0;
		for (int nE=0;nE<=M;nE++){
			for (int poids=0;poids<=P;poids++){
				const double vrais=logBinomExemple[nE]+logBinomSecret[poids];
				double cur=vrais+pSecret+log(binom(N,poids));
				probaMaxi+=exp(max(-50.,cur));
				for (int p2=0;p2<=P;p2++){
					const double pAuDessus=upper_bound(vraisTotale[p2],vraisTotale[p2]+M+1,make_pair(vrais,0.))->second/vraisTotale[p2][0].second;
					cur+=logPuiss(pAuDessus,vraisTotale[p2][0].second-(poids==p2));
//					printf("%d %d %f %f %f\n",nE,poids,vrais,cur,pAuDessus);
				}
				proba+=exp(max(-50.,cur));
			}
			if (nE<M) pSecret+=log((M-nE)/(nE+1.));
		}
		printf("%g %g %g\n",proba,travail/proba,proba/probaMaxi);
//		printf("%d\n",P);
//		return 0;
	}
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
				double cur=log(vraisTotale[poids][0].second)+poids*log(p)+(N-poids)*log(1-p)+logBinomExemple[nE]+pSecret;
				probaMaxi+=exp(max(-50.,cur));
				for (int p2=0;p2<=P+L;p2++){
					const double pAuDessus=upper_bound(vraisTotale[p2],vraisTotale[p2]+M+1,make_pair(vrais,0.))->second/vraisTotale[p2][0].second;
					cur+=logPuiss(pAuDessus,vraisTotale[p2][0].second-(poids==p2));
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
