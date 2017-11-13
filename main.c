#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define KMP() {\
	while( ti < dims[max] && wi < dims[min] )\
        if( secv[min][wi] == secv[max][ti] ){\
		    wi++;\
            ti++;\
        } else if( wi != 0){\
            wi = metaData[min][wi-1];\
        } else {\
            ti++;\
        }\
}

int main(){
    FILE * fi = fopen("adn.in","rt");
    unsigned char n, i, j, cn;
    int rd;
    rd = fscanf( fi, "%hhd\n" , &n);
    if( rd == 0 ){
    	perror("S-au citi 0 caractere !\n");
    	return 1;
    }
    cn = n;
    char **secv = (char**)malloc( n * sizeof(char*) ),
        *aux = (char*)malloc(30001 * sizeof(char));
    unsigned short *dims = (unsigned short *)malloc(n*sizeof(unsigned short));
    n--;
    //Reading stage
    for( i  = 0 ; i < n ; i++ ){
        rd = fscanf(fi,"%s\n",aux);
        if(!rd){
        	perror("S-au citi 0 caractere !\n");
        	for( j = 0 ; j < i ; j++ )
        		free(secv[j]);
        	free(secv);
        	free(aux);
    		return 1;
        }
        dims[i] = strlen(aux);
        secv[i] = (char*)malloc((dims[i]+1) * sizeof(char));
        strcpy(secv[i],aux);
    }
    rd = fscanf(fi,"%s\n",aux);
    if( !rd ){
    	perror("S-au citi 0 caractere !\n");
        for( j = 0 ; j < n ; j++ )
       		free(secv[j]);
       	free(secv);
       	free(aux);
   		return 1;
    }
    dims[n] = strlen(aux);
    secv[n] = (char*)malloc((dims[n]+1) * sizeof(char));
	strcpy(secv[n],aux);

    free(aux);
	close(fi);

    n++;
    unsigned short **metaData = (unsigned short **)malloc(n * sizeof(unsigned short*))
        , ti, wi;
    //KMP pre-processing
    for( i = 0 ; i < n ; i++ ){
        metaData[i] = (unsigned short*)malloc(dims[i] * sizeof(unsigned short));
        metaData[i][0] = 0;
        wi=0;
        ti = 1;
        while( ti < dims[i] ){
            if( secv[i][ti] == secv[i][wi] ){
                wi++;
                metaData[i][ti] = wi;
                ti++;
            } else if( wi > 0 )
                wi = metaData[i][wi-1];
            else {
                metaData[i][ti] = 0;
                ti++;
            }
        }
    }

    unsigned short costs[n][n];
    unsigned char deElim[n], max, min, k = 0;
    for( i = 0 ; i < n ; i++ )
        deElim[i] = 0;
    n--;
    //The cost matrix stores for 2 sequences Si and Sj in cost[i][j]
    //the maximum number of characters that are a sufix for Si and
    //prefix for Sj.
    for( i = 0 ; i < n ; i++ )
        for( j = i + 1 ; j <= n && !deElim[i] ; j++ )
            if( !deElim[j] ){
                if( strlen(secv[i]) > strlen(secv[j]) ){
                    max = i;
                    min = j;
                } else {
                    max = j;
                    min = i;
                }
                wi = ti = 0;
                KMP();
                if( wi == dims[min]){
                    deElim[min] = 1;
                    k++;
                }
                else{
                    costs[max][min] = wi;
                    ti = 1;
                    wi = 0;
                    unsigned char swapAux;
                    swapAux = max;
                    max = min;
                    min = swapAux;

					KMP();

                    costs[max][min] = wi;
                }
            }
    n++;

    for( i = 0 ; i < n ; i++ )
        free(metaData[i]);
    free(metaData);
    free(dims);

    //Some input sequences are entirely contained in other input
    //sequences.
    unsigned char indexes[n-k];
    i = j = 0;
    while( i < n ){
        if(!deElim[i]){
            indexes[j] = i;
            j++;
        }
        i++;
    }
    n = j;

    unsigned int mask, maxLen = 1 << n, a;
    int **variante = (int**)malloc(n*sizeof(int*));
    if( !variante ){
        perror("Nu s-a putut aloca memorie pentru variante !");
        for( i = 0 ; i < cn ; i++ )
                free(secv[i]);
            free(secv);
        return 1;
    }
    unsigned char **after = (unsigned char**)malloc(n*sizeof(unsigned char*));

    for( i = 0 ; i < n ; i++ ){
        variante[i] = (int*)malloc( maxLen * sizeof(int) );
        if( !variante[i] ){
            perror("Nu s-a putut aloca memorie pentru variante !");
            for( j = 0 ; j < i ; j++ )
                free(variante[j]);
            free(variante);
            for( i = 0 ; i < cn ; i++ )
                free(secv[i]);
            free(secv);
            return 1;
        }
        after[i] = (unsigned char*)malloc( maxLen * sizeof(unsigned char) );
        if( !after[i] ){
            perror("Nu s-a putut aloca memorie pentru after !");
            for( j = 0 ; j < i ; j++ )
                free(after[j]);
            free(after);
            for( j = 0 ; j < n ; j++ )
                free(variante[j]);
            free(variante);
            for( i = 0 ; i < cn ; i++ )
                free(secv[i]);
            free(secv);
            return 1;
        }
    }

    for( i = 0 ; i < n ; i++ ){
        for( mask = 0 ; mask < maxLen ; mask++ ){
            variante[i][mask] = -1;
        }
    }

    for( i = 0 ; i < n ; i++ ){
        for( j = 0 ; j < i ; j++ ){
            variante[i][(1<<i)+(1<<j)] = costs[indexes[i]][indexes[j]];
            after[i][(1<<i)+(1<<j)] = j;
        }
        for( j = i+1 ; j < n ; j++ ){
            variante[i][(1<<i)+(1<<j)] = costs[indexes[i]][indexes[j]];
            after[i][(1<<i)+(1<<j)] = j;
        }
    }

    //From here on starts a dynamic programming algorithm in a bottom-up fashion.
    for( mask = 7 ; mask < maxLen ; mask ++ )
        for( i = 0 ; i < n ; i++ )
            if( mask & ( 1 << i ) ){

                a = mask ^ ( 1 << i );

                for( j = 0 ; j < i ; j++ ){
                    if( ( mask & ( 1 << j ) )  &&
                            variante[j][a] != -1 &&
                            variante[i][mask] < costs[indexes[i]][indexes[j]] +
                                variante[j][a]
                        ) {
                        variante[i][mask] = costs[indexes[i]][indexes[j]] +
                            variante[j][a];
                        after[i][mask] = j;
                    }
                }

                for( j = i + 1 ; j < n ; j++ ){
                    if( ( mask & ( 1 << j ) ) &&
                            variante[j][a] != -1 &&
                            variante[i][mask] < costs[indexes[i]][indexes[j]] +
                                variante[j][a]
                        ){
                        variante[i][mask] = costs[indexes[i]][indexes[j]] +
                            variante[j][a];
                        after[i][mask] = j;
                    }
                }
            }
    maxLen--;
    mask = variante[0][maxLen];
    j = 0;
    for( i = 1 ; i < n ; i++ )
        if( variante[i][maxLen] > mask ){
            mask = variante[i][maxLen];
            j = i;
        }


    for( i = 0 ; i < n ; i++ )
        free(variante[i]);
    free(variante);
    deElim[1] = n;

    mask = maxLen;
    ti = 0;
    fi = fopen("adn.out","wt");
    while(1){
        fprintf(fi,"%s", secv[indexes[j]] + ti );
        n--;
        if( n == 0 )
            break;
        deElim[0] = j;
        j = after[j][mask];
        ti = costs[indexes[deElim[0]]][indexes[j]];
        mask -= (1 << deElim[0]);
    }
    fclose(fi);
    for( i = 0 ; i < deElim[1] ; i++ )
        free(after[i]);
    free(after);
    for( i = 0 ; i < cn ; i++ )
        free(secv[i]);
    free(secv);

    return 0;
}
