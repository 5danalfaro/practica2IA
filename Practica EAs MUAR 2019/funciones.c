/*============================================================================
| INCULIR ESTAS INSTRUCCIONES AL PRINCIPIO DEL PROGRAMA
============================================================================*/
#include "ga.h" 

#define MAXR 10

 /*--- Forward declarations ---*/  
void chrom2chessboard(Chrom_Ptr c, char tablero[MAXR][MAXR], int n);
int obj_fun();
int cuentaamenazas(char tablero[MAXR][MAXR], int n);

// Variables globales
int tipo;
int size;
char tablero[MAXR][MAXR];


/*----------------------------------------------------------------------------
| FUNCIONES
----------------------------------------------------------------------------*/
   
/*----------------------------------------------------------------------------
| obj_fun() - user specified objective function - FITNESS
----------------------------------------------------------------------------*/
int obj_fun(Chrom_Ptr chrom) 
{
  int i, j; 
  double val = 0.0;
  int amenazas = 0;
  int reinas[MAXR], sobran = 0, combinacion = 0;

  chrom2chessboard(chrom, tablero, size);
  
  amenazas = cuentaamenazas(tablero, size);
  
  if(tipo==DT_BIT)
    {
	for(i=0;i<size;i++)
		for(j=0;j<size;j++){
			if (tablero[i][j] == 1)
				++reinas[i];
	}

      
	for(i=0;i<size;i++)
		if (reinas[i] > 1) ++sobran;

	chrom->fitness = size * 10 - 0.05 * amenazas + 10 * sobran; 
    }
  

  else if(tipo==DT_INT_PERM)
    {
	chrom->fitness = size - amenazas;
	
	/*
	if (chrom->fitness == 10){
		for (j=0; j<size; ++j){
			chrom->gene[j] = j+1;
		}		
		++combinacion;	
	}
 	*/
    }

  else printf("Te eta equivocando chico\n");

  
  return 0;
  
}

/*----------------------------------------------------------------------------
| Funciones proporcionadas
----------------------------------------------------------------------------*/

void chrom2chessboard(Chrom_Ptr c, char tablero[MAXR][MAXR], int size)
{
  int i,j;

  if(tipo==DT_BIT)
    {
      for(i=0;i<size;i++)
	for(j=0;j<size;j++)
	  tablero[i][j] = c->gene[i*size+j];
	
    }
  else
    if(tipo==DT_INT_PERM)
      {
	for(i=0;i<size;i++)
	  for(j=0;j<size;j++)
	    tablero[i][j] = 0;
	  
	for(i=0;i<size;i++)
	  tablero[i][(int)c->gene[i]-1] = 1;
      }

    else printf("WTF???\n");
}


   
int cuentaamenazas(char tablero[MAXR][MAXR], int n)
{
  int i,sum,f,c;
  int amenazas = 0;

  // checking rows
  for(f=0;f<n;f++)
    {
      sum = 0;
      for(c=0;c<n;c++)
	sum += (int)tablero[f][c];
      if(sum>1) amenazas += sum*(sum-1)/2;
    }

  // checking columns
  for(c=0;c<n;c++)
    {
      sum = 0;
      for(f=0;f<n;f++)
	sum += (int)tablero[f][c];
      if(sum>1) amenazas += sum*(sum-1)/2;
    }

  // checking diagonals 
  for(f=0;f<n;f++)
    for(c=0;c<n;c++)
      {
	sum = 0;
	for(i=0;i<n;i++)
	  {
	    if(f+i<n && c+i<n)
	      sum += (int)tablero[f+i][c+i];
	  }
	if(sum>1) amenazas += sum*(sum-1)/2;
	  

	sum = 0;
	for(i=0;i<n;i++)
	  {
	    if(f+i<n && c-i>=0)
	      sum += (int)tablero[f+i][c-i];
	    }
	if(sum>1) amenazas += sum*(sum-1)/2; 


	sum = 0;
	for(i=0;i<n;i++)
	  { if(f-i>=0 && c+i<n)
	      sum += (int)tablero[f-i][c+i];
	   }
	if(sum>1) amenazas += sum*(sum-1)/2;
	
	sum = 0;
	for(i=0;i<n;i++)
	  {if(f-i>=0 && c-i>=0)
	      sum += (int)tablero[f-i][c-i];
	    }
	if(sum>1) amenazas += sum*(sum-1)/2;

      }

  return amenazas;
}



/*----------------------------------------------------------------------------
| MAIN
----------------------------------------------------------------------------*/

int main() 
{
   GA_Info_Ptr ga_info;
   int i, aux;

   printf("Introduce tipo de dato (1: permutacion; otro: bit string)");
   scanf("%d", &aux);

   if (aux == 1){

   /*--- Initialize the genetic algorithm ---*/
   ga_info = GA_config("GAconfig_permut", obj_fun);

   }

   else {

   /*--- Initialize the genetic algorithm ---*/
   ga_info = GA_config("GAconfig_bit", obj_fun);

   }

 
   tipo = ga_info->datatype;
   
   if(tipo==DT_BIT)
      size = (int)sqrt(ga_info->chrom_len);
   else if(tipo==DT_INT_PERM)
      size = (int)ga_info->chrom_len;
   else {printf("Something went wrong...%d\n",tipo);exit(-1);}


   /*--- Run the GA ---*/
   GA_run(ga_info);


   printf("\nBest chrom:  \n");
   for(i=0;i<ga_info->chrom_len;i++)
	//if (i%10 == 0) printf("\n");	     
	printf("%5.0f  ",ga_info->best->gene[i]);
   
   printf("   (fitness: %g)\n\n",ga_info->best->fitness);

   return 0;

}


