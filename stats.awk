#!/bin/awk

NR>1 { 
    if ($4 == $12 && $1 == $9) { 
	value[$4]+=$16/$8; 
	time[$4]+=$15; 
	num[$4]+=1; 
	minv[$4]=(minv[$4] > 0 && minv[$4]<$16/$8)?minv[$4]:$16/$8;
	maxv[$4]=(maxv[$4]>$16/$8)?maxv[$4]:$16/$8;
	mint[$4]=(mint[$4] > 0 && mint[$4]<$15)?mint[$4]:$15;
	maxt[$4]=(maxt[$4]>$15)?maxt[$4]:$15; 
    }
}

END{ 
    print "\nmin value\n"; 
    for (N in value) { 
	printf("(%d,%g)\n", N, minv[N]) 
    } 
    print "\navg value\n";
    for (N in value) { 
	printf("(%d,%g)\n", N, value[N]/num[N])
    } 
    print "\nmax value\n"; 
    for (N in value) { 
	printf("(%d,%g)\n", N, maxv[N]) 
    } 
    print "\nmin time\n"; 
    for (N in time) { 
	printf("(%d,%g)\n", N, mint[N]) 
    } 
    print "\navg time\n"; 
    for (N in time) { 
	printf("(%d,%g)\n", N, time[N]/num[N]) 
    } 
    print "\nmax time\n"; 
    for (N in time) { 
	printf("(%d,%g)\n", N, maxt[N]) 
    } 
}
