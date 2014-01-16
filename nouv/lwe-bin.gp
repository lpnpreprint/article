log2(x)=log(x)/log(2);
for(N=7,11,\
	n=2^N;\
	q=n^2.;\
	Ve=n^3/log2(n)^4.;\
	Vs=1/4.; /* binaire */\
	a=truncate(log2(q^2/Ve/(4/log(2)*Pi^2)*2*log(2)*n/log2(log2(n))),&e);\
	C=log2(q^2*n*Vs/Ve/12)/a;\
	K=n/(2*log(1+1/(C-1)));\
	a=truncate(log2(q^2/Ve/(4/log(2)*Pi^2)*K),&e);\
	C=log2(q^2*n*Vs/Ve/12)/a;\
	K=n/(2*log(1+1/(C-1)));\
	V=2^a*Ve;\
	Vini=V;\
	lTot=0;\
	for(i=0,a-1,\
		l=truncate(1+2*K/(C*a-i),&e);\
		lTot+=l;\
		p=truncate(1+q/2^(K/l),&e);\
		V+=2^(a-i)*(p^2-1)/12.*l*Vs;\
/*		print(l/(2*K/C/a)," ",2^(a-i)*(p^2-1)/12.*l*Vs/Vini);*/\
	);\
	print(K/n," ",lTot/n+0.," ",a," ",K/(2*Pi^2/log(2)*V/q^2)," ",V/Vini," ",1/(C-1));\
);
quit
