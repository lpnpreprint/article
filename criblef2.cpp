#include <cstdio>
#include <bitset>
#include <algorithm>
#include <cmath>

using namespace std;

const int N=100;
const int K=1<<14;

typedef bitset<N> Vecteur;

Vecteur iter[K];
Vecteur nouv[2*K-2];

bool comp (const Vecteur &a,const Vecteur &b){
	return a.count()<b.count();
}

bool lexico(const Vecteur &a,const Vecteur &b){
	int i;
	for (i=0;i<N && a[i]==b[i];i++) ;
	return i<N && b[i];
}

void iteration(){
	for (int i=0;i<K;i++){
		const int debut=i ? K+min(i-2,0) : 0;
		for (int j=i+1;j<K;j++)
			nouv[debut+j-i-1]=iter[i]^iter[j];
		if (i){
/*			for (int j=0;j<debut+K-i-1;j++)
				printf("%s\n",nouv[j].to_string().c_str());
			puts("");*/
			sort(nouv,nouv+debut+K-i-1,lexico);
/*			for (int j=0;j<debut+K-i-1;j++)
				printf("%s\n",nouv[j].to_string().c_str());
			puts("");*/
			Vecteur *fin=unique(nouv,nouv+debut+K-i-1);
/*			for (int j=0;j<fin-nouv;j++)
				printf("%s\n",nouv[j].to_string().c_str());
			puts("");*/
			for (int id=fin-nouv;id<K;id++)
				for (int j=0;j<N;j++)
					nouv[id][j]=rand()&1;
			nth_element(nouv,nouv+K,fin,comp);
/*			for (int j=0;j<K;j++)
				printf("%d %s\n",nouv[j].count(),nouv[j].to_string().c_str());
			puts("---");*/
		}
	}
	for (int i=0;i<K;i++)
		iter[i]=nouv[i];
}

double H(double x){
	return (-x*log(x)-(1-x)*log(1-x))/log(2);
}

int main(){
	for (int i=0;i<K;i++)
		for (int j=0;j<N;j++)
			iter[i][j]=rand()&1;
	int ancSomme;
	for (int i=0;i<40;i++){
		int somme=0;
		for (int j=0;j<K;j++){
			somme+=iter[j].count();
//			printf("%s\n",iter[j].to_string().c_str());
		}
		printf("%d %f %g %f\n",i+1,somme/(K+0.)/N,i ? (ancSomme-somme)/(K+0.)/N : 0.,i ? N*(H((1-pow(1-2*ancSomme/(K+0.)/N,2))/2)-H(somme/(K+0.)/N)) : 0);
		fflush(stdout);
		ancSomme=somme;
		iteration();
	}
	return 0;
}
