


// VARIABLES GLOBALES
int cities[500][2];
double **DISTANCES;
int NN;




// CARGA LAS CIUDADES DEL FICHERO Y CALCULA LA MATRIZ DE DISTANCIAS
int load_inst(char *fn)
{
 FILE *fp;
 char tmp[50];
 int i,j;

  if(!(fp=fopen(fn,"r")))
{return -1;}

//  read out 3 header lines
fgets(tmp,50,fp);
fgets(tmp,50,fp);
fgets(tmp,50,fp);
// read dimension
fscanf(fp,"%s %d",tmp,&NN);
//  read out 3 header lines   ¡¡¡ OJO algunos ficheros solo tienen 2!!!
fgets(tmp,50,fp);
fgets(tmp,50,fp);
fgets(tmp,50,fp);

 for(i=0;i<NN;i++)
     fscanf(fp,"%d %d %d",&j,&cities[i][0],&cities[i][1]);
     
 fclose(fp);

 printf("Read %d cities out of %d\n",i,NN);

 DISTANCES = (double**)malloc(NN*sizeof(double*));

for (i=0; i<NN; i++) 
         DISTANCES[i] = (double *)malloc(NN * sizeof(double)); 
 
 for(i=0;i<NN;i++)
   for(j=0;j<NN;j++)     
        DISTANCES[i][j]=sqrt((cities[i][0]-cities[j][0])*(cities[i][0]-cities[j][0])+(cities[i][1]-cities[j][1])*(cities[i][1]-cities[j][1]));
     

 return 1;
}
