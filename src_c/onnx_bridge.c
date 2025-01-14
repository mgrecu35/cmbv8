#include "stdlib.h"
#include "math.h"

#include <npy_array_list.h>
#include <stdlib.h>
#include <assert.h>

#pragma GCC diagnostic ignored "-Wmissing-field-initializers"

void get_dprdata_c_loc_(float *zku,float *zka,int *rain_type,
				    int *node, int *bin_sfc, float *freez_h);

void onnx_bridge(int ichunk, int orbnumber)
{
  float *zku, *zka;
  int *rain_type,*node, *bin_sfc;
  float *freez_h;
  char *fname;
  printf("ichunk=%i iorb=%i\n",ichunk,orbnumber);
  zku=(float *) malloc(sizeof(float)*300*49*88);
  zka=(float *) malloc(sizeof(float)*300*49*88);
  rain_type=(int *)malloc(sizeof(int)*300*49);
  node=(int *)malloc(sizeof(int)*300*49*5);
  bin_sfc=(int *)malloc(sizeof(int)*300*49);
  freez_h=(float *)malloc(sizeof(float)*300*49);
  
  get_dprdata_c_loc_(zku,zka,rain_type,node,bin_sfc,freez_h);
  fname=(char *)malloc(sizeof(char)*100);
  sprintf(fname, "npy_tmp_dir/cmb_npy_%2.2i_%6.6i.npz", ichunk, orbnumber);
  printf("Filename: %s\n", fname);

  npy_array_list_t* list_var=NULL;
  list_var= npy_array_list_append(list_var, 
				  NPY_ARRAY_BUILDER_COPY( rain_type, SHAPE( 300, 49 ), NPY_DTYPE_INT32 ), "rain_type" );
  list_var= npy_array_list_append(list_var,
				  NPY_ARRAY_BUILDER_DEEPCOPY( zku, SHAPE( 300, 49, 88 ), NPY_DTYPE_FLOAT32 ), "zku" );

  npy_array_list_save_compressed( fname, list_var, ZIP_CM_DEFAULT, 0 );
  npy_array_list_free(list_var);
  
  free(zku);
  free(zka);
  free(rain_type);
  free(node);
  free(bin_sfc);
  free(freez_h);
}
