/*============================================================================
| INCULIR ESTAS INSTRUCCIONES AL PRINCIPIO DEL PROGRAMA
============================================================================*/
 

#define MAXR 10

 /*--- Forward declarations ---*/  
void chrom2chessboard(Chrom_Ptr c, char tablero[MAXR][MAXR], int n);

// Variables globales
int tipo;
int size;
char tablero[MAXR][MAXR];




/*----------------------------------------------------------------------------
| AÑADIR ESTAS INSTRUCCIONES EN EL MAIN
----------------------------------------------------------------------------*/
   
    
   tipo = ga_info->datatype;
   
   if(tipo==DT_BIT)
     size = (int)sqrt(ga_info->chrom_len);
   else if(tipo==DT_INT_PERM)
      size =  (int)ga_info->chrom_len;
    else {printf("Something went wrong...%d\n",tipo);exit(-1);}
     
  
}




/*----------------------------------------------------------------------------
| FUNCIONES
----------------------------------------------------------------------------*/
   

void chrom2chessboard(Chrom_Ptr c, char tablero[MAXR][MAXR], int n)
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

