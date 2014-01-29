
#include <cstdio>
#include <cmath>
#include <set>
#include <algorithm>

using namespace std;

const int n=1<<15;
double q=n*(n+0.),Ve,Vs;

double log2(double x){
	return log(x)/log(2);
}

double sommeCarre(double x){
	return x*(x+1)*(2*x+1)/6;
}

double var(int L,int l){
	const int p=ceil(q/pow(2,(L+0.)/l));
	return (sommeCarre(p/2)+sommeCarre((p-1)/2))/p*l*Vs;
}

double trouveR(double y){ // r+ln(1-r)=-y via Newton
	double maxx=1-exp(-y-1);
	double x=min(sqrt(2*y),maxx);
	if (x>1e-3 && y<30){
		x=min(maxx,x-(x+log(1-x)+y)/(1-1/(1-x)));
		x=min(maxx,x-(x+log(1-x)+y)/(1-1/(1-x)));
		x=min(maxx,x-(x+log(1-x)+y)/(1-1/(1-x)));
		x=min(maxx,x-(x+log(1-x)+y)/(1-1/(1-x)));
		if (abs(x+log(1-x)+y)/y>1e-2)
			printf("Err : y=%f x=%f err=%f\n",y,x,x+log(1-x)+y);
	}
	return x;
}

double eval(int *longueur,int a,int L,int N){
	double somme=0;
	int curL=0;
	for (int i=0;i<a;i++){
		curL+=longueur[i];
		somme=(somme+var(L,longueur[i]))*2*(1-trouveR(2*N*log(2)/curL));
	}
	return somme+pow(2,a)*Ve;
}

bool possible(int a,int L,int N){
	int *longueur=new int[a];
	for (int i=0;i<a;i++){
		if (i<a-1) longueur[i]=n/(a+1);
		else longueur[i]=n-(a-1)*(n/(a+1))-(L+N)/log2(q)+2;
	}
	double cur=eval(longueur,a,L,N);
	int nbEchec=0;
	for (int iter=0;iter<n*n;iter++){
		int i=rand()%a,j=rand()%a;
		if (i==j || longueur[i]==0)
			continue;
		longueur[i]--;
		longueur[j]++;
		const double nCur=eval(longueur,a,L,N);
		if (2*M_PI*M_PI*nCur/log(2)/q/q+4<L+N){ // on peut distinguer
/*			int curN=0;
			for (int i=0;i<a;i++){
				curN+=longueur[i];
				printf("%d %f %f\n",longueur[i],trouveR(2*N*log(2)/curN),2*sqrt(N*log(2)/curN));
			}
			printf("\n%f : %f %f\n",nCur,2*M_PI*M_PI*nCur/log(2)/q/q,var(L,10));*/
			delete []longueur;
			return true;
		}
		if (nCur<cur){
			cur=nCur;
			nbEchec=0;
		}
		else{
			longueur[i]++;
			longueur[j]--;
			nbEchec++;
			if (nbEchec==n*5)
				break;
		}
	}
//	printf("%d %d : %f\n",a,L,cur/q/q);
	delete []longueur;
	return false;
}

pair<int,int> trouveL(int debuta,int fina,int debutL,int finL,int N,int pasVoulu=-1){
	if (debutL+1==finL)
		return make_pair(debuta,finL);
	int pas=pasVoulu>0 ? pasVoulu : max((finL-debutL)/(fina-debuta+1),1);
//	printf("[%d;%d[ ]%d;%d] : %d\n",debuta,fina,debutL,finL,pas);
	for (int a=debuta;a<fina;a++)
		while (finL>pas && possible(a,finL-pas,N)){
			finL-=pas;
//			printf("%d,%d\n",a,finL);
			debuta=a;
		}
	while (debuta+1<fina && !possible(fina-1,finL,N)) fina--;
	return trouveL(debuta,fina,max(0,finL-pas),finL,N);
}


int main(){
	Ve=q*n/pow(log2(n),4);
	Vs=1/4.;
	int a=max(1.,2*log2(n)-2);
	int L=2*n;	
	int N=0;
	srand(0);
//	printf("%d\n",possible(29,1015,150));
//	return 0;
	for (N=400;N<500;N++){
	pair<int,int> param=trouveL(a,log2(L*q*q/Ve/2/M_PI/M_PI*log(2)),0,L,N);
	a=param.first;
	L=param.second;
	const double travail=L+log2(n*(a+1.)*q)+2*N;
	printf("%.0f %d %d %d %f %f\n",log2(n),N,a,L,travail,travail/n);
	}
	return 0;
}
