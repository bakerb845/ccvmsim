#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "mesh.h"

int layeredMesh_driver(struct cvm_parms_struct *parms,
                       struct mesh_struct *mesh)
{
    memset(mesh, 0, sizeof(struct mesh_struct));
    // First fix the number of x and y grid points in the top layer
    return 0;
}
