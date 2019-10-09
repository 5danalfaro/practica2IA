

//VARIABLES GLOBALES
int **ADJ=NULL;
int NN;





int load_inst(char *fn)
{
 FILE *fp;
 char tmp[100];
 int i,j;
 char tmpc,tmpe;
 
 if(!(fp=fopen(fn,"r")))
   return -1;
 
 while(!feof(fp))
   {
   
     fscanf(fp,"%c",&tmpc);
     switch(tmpc)
       {
       case 'c':
	 fgets(tmp,50,fp);
	 break;
       case 'p':
	 fscanf(fp,"%s %d %d",tmp,&NN,&j);
	  ADJ = (int**)malloc(NN*sizeof(int*));
	  for (i=0; i<NN; i++) 
	    ADJ[i] = (int *)malloc(NN * sizeof(int));
	  for (i=0; i<NN; i++) 
	    for (j=0; j<NN; j++)
	      ADJ[i][j]=0;
	 break;
       case 'e':
	 fscanf(fp,"%d %d",&i,&j);
	   
	 if(ADJ)
	   {ADJ[i-1][j-1]=1; ADJ[j-1][i-1]=1;}
	 break;  
	 default:
	 //	 printf("Unknown label...\n");
	 //return -1;
	   break;
       }

   }
 

 return 1;
}
