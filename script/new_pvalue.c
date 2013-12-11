#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    
    int exp_count;
    double pvalue;
//    exp_count=atoi(argv[1]);
    pvalue=atof(argv[1]);
   
    if (pvalue>=1){

    pvalue=0.999;  
    printf ("%.3e\n", pvalue);

    }

    else

    printf ("%.3e\n", pvalue);

    return 0;
}

