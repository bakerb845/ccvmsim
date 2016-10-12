#ifndef _h5_cinter_h__
#define _h5_cinter_h__ 1

#include <stdbool.h>
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif
#include <hdf5.h>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef __cplusplus
extern "C" 
{
#endif

hid_t h5_open_rdonly(const char *flname);
hid_t h5_open_rdwt(const char *flname);
int h5_close(const hid_t file_id);
int h5_write_array__float(const char *dset_name, const hid_t file_id,
                          const int n, const float *x);
int h5_write_array__double(const char *dset_name, const hid_t file_id,
                           const int n, const double *x);
int h5_write_array__int(const char *dset_name, const hid_t file_id,
                        const int n, const int *x);
int h5_write_array__chars(const char *citem_chr, const hid_t file_id,
                          const int n, const char **c);
int h5_read_array__double(const char *dset_name, const hid_t file_id,
                          const int nref, double *x);
int h5_read_array__float(const char *dset_name, const hid_t file_id,
                         const int nref, float *x);
int h5_read_array__int(const char *dset_name, const hid_t file_id,
                       const int nref, int *x);
int h5_write_attribute__double(const char *citem, const hid_t hdf5_id,
                               const int n, const double *attr_data);
int h5_write_attribute__int(const char *citem, const hid_t hdf5_id,
                            const int n, const int *attr_data);
int h5_write_attribute__char(const char *citem, const hid_t hdf5_id,
                             const int n, const char **cattr);
int h5_n_group_members(const char *group_name, const hid_t file_id);
int h5_get_array_size(const hid_t file_id, const char *citem);
bool h5_item_exists(const hid_t file_id, const char *citem_in);
herr_t h5_create_group(const hid_t file_id, const char *cgroup);

#ifdef __cplusplus
}
#endif
#endif /* __H5_C_INTER_H__ */
