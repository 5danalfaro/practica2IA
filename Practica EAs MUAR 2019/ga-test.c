/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
|
| Genetic Algorithm Test Program 
============================================================================*/
#include "ga.h"

int obj_fun();    /*--- Forward declaration ---*/

/*----------------------------------------------------------------------------
| main()
----------------------------------------------------------------------------*/
main() 
{
   GA_Info_Ptr ga_info;
   int i;

   /*--- Initialize the genetic algorithm ---*/
   ga_info = GA_config("GAconfig_ejemplo", obj_fun);

   /*--- Run the GA ---*/
   GA_run(ga_info);

   printf("\nBest chrom:  ");
   for(i=0;i<ga_info->chrom_len;i++)
	     printf("%5.4f  ",ga_info->best->gene[i]);
   
   printf("   (fitness: %g)\n\n",ga_info->best->fitness);

}

/*----------------------------------------------------------------------------
| obj_fun() - user specified objective function
----------------------------------------------------------------------------*/
int obj_fun(Chrom_Ptr chrom) 
{
  int i; 
  double val = 0.0;
  
  for(i = 0; i < chrom->length; i++)
    {        
      val +=  chrom->gene[i];
    }
  
  chrom->fitness = val;
  
  return 0;
  
}
