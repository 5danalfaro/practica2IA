/*============================================================================
| (c) Copyright Arthur L. Corcoran, 1992, 1993.  All rights reserved.
|
| Mutation operators 
|
| Bit String Representations
|    MU_simple_invert() - random bit inversion based on mutation rate
|    MU_simple_random() - random bit value based on mutation rate
|
| Any Representation
|    MU_swap()          - random element swap based on mutation rate
|
| Interface
|    MU_table[]   - used in selection of mutation method
|    MU_set_fun() - set and select user defined mutation function
|    MU_select()  - select mutation function by name
|    MU_name()    - get name of current mutation function
|    MU_fun()     - setup and perform current mutation operator
============================================================================*/
#include "ga.h"

int MU_simple_invert(), MU_simple_random(), MU_swap();
 /* rnd float in [0..1] -- introduced by claudio 10/02/2004 */
int MU_float_random(), MU_float_rnd_pert(), MU_float_LS(), MU_float_gauss_pert();

double gaussian_random(void);

/*============================================================================
|                     Mutation interface
============================================================================*/
/*----------------------------------------------------------------------------
| Mutation table
----------------------------------------------------------------------------*/
FN_Table_Type MU_table[] = {
   {  NULL,           NULL             }, /* user defined function */
   { "simple_invert", MU_simple_invert },
   { "simple_random", MU_simple_random },
   { "swap",          MU_swap          },
   { "float_random",  MU_float_random  },
   { "float_rnd_pert",MU_float_rnd_pert},
   { "float_LS",      MU_float_LS      },
   { "float_gauss_pert",MU_float_gauss_pert},
   { NULL,            NULL             }
};

/*----------------------------------------------------------------------------
| Select user mutation method
----------------------------------------------------------------------------*/
MU_set_fun(ga_info, fn_name, fn_ptr)
   GA_Info_Ptr  ga_info;
   char         *fn_name;
   FN_Ptr       fn_ptr;
{
   return FN_set_fun(ga_info, MU_table, fn_name, fn_ptr, &ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Select mutation method
----------------------------------------------------------------------------*/
MU_select(ga_info, fn_name)
   GA_Info_Ptr  ga_info;
   char         *fn_name;
{
   return FN_select(ga_info, MU_table, fn_name, &ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Mutation name
----------------------------------------------------------------------------*/
char *MU_name(ga_info)
   GA_Info_Ptr  ga_info;
{
   return FN_name(ga_info, MU_table, ga_info->MU_fun);
}

/*----------------------------------------------------------------------------
| Mutation interface
----------------------------------------------------------------------------*/
MU_fun(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr   chrom;
{
   /*--- Random chance to mutate ---*/
   if(RAND_FRAC() <= ga_info->mu_rate && ga_info->MU_fun != NULL) {
      ga_info->MU_fun(ga_info, chrom);
      ga_info->num_mut++;
      ga_info->tot_mut++;
   }
}

/*============================================================================
|                             Mutation operators
============================================================================*/
/*----------------------------------------------------------------------------
| Simple mutation, invert bit
----------------------------------------------------------------------------*/
MU_simple_invert(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   int idx;

   /*--- Select bit at random ---*/
   idx = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Invert selected bit ---*/
   chrom->gene[idx] = chrom->gene[idx] ? 0 : 1;
}

/*----------------------------------------------------------------------------
| Simple mutation, random bit
----------------------------------------------------------------------------*/
MU_simple_random(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   int idx;

   /*--- Select bit at random ---*/
   idx = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Assign random value to bit ---*/
   chrom->gene[idx] = RAND_BIT();
}

/*----------------------------------------------------------------------------
| Swap random elements
----------------------------------------------------------------------------*/
MU_swap(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   Gene_Type tmp;
   int       i, j;

   /*--- Select two bits at random (can be same) ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);
   j = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Swap the elements ---*/
   tmp            = chrom->gene[i];
   chrom->gene[i] = chrom->gene[j];
   chrom->gene[j] = tmp;
}


/*----------------------------------------------------------------------------
|   rnd float perturbation in [0..1] -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_rnd_pert(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   int       i;

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);

   if(i==6)
     {
       /*--- Generate random element ---*/
       chrom->gene[i] = RAND_FRAC();
     }
   else
     {
       /*--- Generate randomly perturbed element ---*/
       chrom->gene[i] += ga_info->pert_range*(1.0 - 2.0*RAND_FRAC());
     }
   //   printf("gene %d, bias %g\n",i,ga_info->mut_bias[i]);
   
   if( chrom->gene[i]>1)
     chrom->gene[i]=1;
   if( chrom->gene[i]<0)
     chrom->gene[i]=0;
}


/*----------------------------------------------------------------------------
|   rnd float in [0..1] -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_random(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   int       i;

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Generate random element ---*/
   chrom->gene[i] = RAND_FRAC();
   if( chrom->gene[i]>1)
     chrom->gene[i]=1;
   if( chrom->gene[i]<0)
     chrom->gene[i]=0;

   chrom->gene[i] = chrom->gene[i]*2-1;

}

/*----------------------------------------------------------------------------
|   Local Search ***experimental*** -- introduced by claudio 10/02/2004 
----------------------------------------------------------------------------*/
MU_float_LS(ga_info, chrom)
     GA_Info_Ptr ga_info;
     Chrom_Ptr chrom;
{

  int   i;
  float old_fit, new_fit, prev_fit, prev_val,tmp;
  

  while(1)
    {
      old_fit=chrom->fitness;
      
      /*--- loop thru genes ---*/
      for(i=chrom->idx_min;i<chrom->length;i++)
	{
	  prev_fit=chrom->fitness;
	  prev_val= chrom->gene[i];
	  
	  chrom->gene[i] += 0.1*(1.0 - 2.0*RAND_FRAC());
	  if( chrom->gene[i]>1)
	    chrom->gene[i]=1;
	  if( chrom->gene[i]<0)
	    chrom->gene[i]=0;
	  
	  obj_fun(chrom);
	  new_fit=chrom->fitness;
	  
	  if(new_fit>prev_fit)
	    {
	      chrom->gene[i]=prev_val;
	      chrom->fitness=prev_fit;
	    }
	}
      if(chrom->fitness>=old_fit)
	break;
      

    }  

}



/*----------------------------------------------------------------------------
|   gaussian perturbation in [0..1] -- introduced by claudio 28/10/2004 
----------------------------------------------------------------------------*/
MU_float_gauss_pert(ga_info, chrom)
   GA_Info_Ptr ga_info;
   Chrom_Ptr chrom;
{
   int       i;

   //int q=15;
   int c1,c2;
   float c3, pert;

   
   c1 = 32767;  // c1=(1 << q)-1;
   c2 = 10922;  // c2=(c1 / 3);
   c3 = 3.05185e-05;  // c3=1.0/c1;
   

   // NB pert in [-1,1], gaussian: mean=0, var=1
   pert=(2*(RAND_DOM(0,c2)+RAND_DOM(0,c2)+RAND_DOM(0,c2))-3*(c2))*c3;

   //pert=gaussian_random();

   /*--- Select one element at random ---*/
   i = RAND_DOM(chrom->idx_min, chrom->length-1);

   /*--- Generate randomly perturbed element ---*/
   chrom->gene[i] += ga_info->pert_range*pert;

   //   printf("gene %d, bias %g\n",i,ga_info->mut_bias[i]);
   

}



double gaussian_random(void)
{
  static int next_gaussian = 0;
  static double saved_gaussian_value;

  double fac, rsq, v1, v2;

  if (next_gaussian == 0) {
    do {
      v1 = 2.0*RAND_FRAC()-1.0;
      v2 = 2.0*RAND_FRAC()-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    saved_gaussian_value=v1*fac;
    next_gaussian=1;
    return v2*fac;
  } else {
    next_gaussian=0;
    return saved_gaussian_value;
  }
}





