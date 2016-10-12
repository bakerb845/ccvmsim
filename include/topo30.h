#include "cvm_struct.h"
#ifndef _topo30_h__
#define _topo30_h__ 1
#ifdef __cplusplus
extern "C" 
{
#endif

int topo30_deformMesh(const char *topofl, int utm_zone,
                      double ztopo_min,
                      struct mesh_struct *mesh);

#ifdef __cplusplus
}
#endif
#endif /* __TOPO30_H__ */
