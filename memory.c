#include <stdlib.h>
#include <string.h>
#include <iscl/log/log.h>
#include "cvm.h"

//============================================================================//
/*!
 * @brief Frees the parameter structure 
 *
 * @param[in,out] parms     parameter structure to be freed
 *
 * @author Ben Baker, ISTI
 *
 */
void cvm_memory_free__parms(struct cvm_parms_struct *parms)
{
    if (parms->dxl_cvm     != NULL){free(parms->dxl_cvm);}
    if (parms->dyl_cvm     != NULL){free(parms->dyl_cvm);}
    if (parms->dzl_cvm     != NULL){free(parms->dzl_cvm);}
    if (parms->z0_cvm      != NULL){free(parms->z0_cvm);}
    if (parms->zcoarsen    != NULL){free(parms->zcoarsen);}
    if (parms->nxl_cvm     != NULL){free(parms->nxl_cvm);}
    if (parms->nyl_cvm     != NULL){free(parms->nyl_cvm);}
    if (parms->nzl_cvm     != NULL){free(parms->nzl_cvm);}
    if (parms->coarsen     != NULL){free(parms->coarsen);}
    memset(parms, 0, sizeof(struct cvm_parms_struct));
    return;
}
//============================================================================//
/*!
 * @brief Frees the CVM model structure
 *
 * @param[in] nlay            number of layers in model
 *
 * @param[in,out] cvm_model   CVM model to free
 *
 * @author Ben Baker, ISTI
 */
void cvm_memory_free__model(int nlay, struct cvm_model_struct *cvm_model)
{
    int k;
    if (cvm_model == NULL){return;}
    for (k=0; k<nlay; k++){
        if (cvm_model[k].vp    != NULL){free(cvm_model[k].vp);}
        if (cvm_model[k].vs    != NULL){free(cvm_model[k].vs);}
        if (cvm_model[k].dens  != NULL){free(cvm_model[k].dens);}
        if (cvm_model[k].Qp    != NULL){free(cvm_model[k].Qp);}
        if (cvm_model[k].Qs    != NULL){free(cvm_model[k].Qs);}
        if (cvm_model[k].xlocs != NULL){free(cvm_model[k].xlocs);}
        if (cvm_model[k].ylocs != NULL){free(cvm_model[k].ylocs);}
        if (cvm_model[k].zlocs != NULL){free(cvm_model[k].zlocs);}
        memset(&cvm_model[k], 0, sizeof(struct cvm_model_struct));
    }
    free(cvm_model);
    cvm_model = NULL;
    return;
}
