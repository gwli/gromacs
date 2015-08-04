/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 *  \brief Define functions for detection and initialization for OpenCL devices.
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 */

#include "gmxpre.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <memory.h>

#include "gromacs/gmxlib/gpu_utils/gpu_utils.h"
#include "gromacs/gmxlib/gpu_utils/ocl_compiler.h"
#include "gromacs/gmxlib/ocl_tools/oclutils.h"
#include "gromacs/legacyheaders/types/enums.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Helper macro for error handling */
#define CALLOCLFUNC_LOGERROR(func, err_str, retval) { \
        cl_int opencl_ret = func; \
        if (CL_SUCCESS != opencl_ret) \
        { \
            sprintf(err_str, "OpenCL error %d", opencl_ret); \
            retval = -1; \
        } \
        else{ \
            retval = 0; } \
}


/*! \brief Helper function that checks whether a given GPU status indicates compatible GPU.
 *
 * \param[in] stat  GPU status.
 * \returns         true if the provided status is egpuCompatible, otherwise false.
 */
static bool is_compatible_gpu(int stat)
{
    return (stat == egpuCompatible);
}

/*! \brief Returns true if the gpu characterized by the device properties is
 *  supported by the native gpu acceleration.
 * \returns             true if the GPU properties passed indicate a compatible
 *                      GPU, otherwise false.
 */
static int is_gmx_supported_gpu_id(struct gmx_device_info_t *ocl_gpu_device)
{
    /* Only AMD and NVIDIA GPUs are supported for now */
    if ((OCL_VENDOR_NVIDIA == ocl_gpu_device->vendor_e) ||
        (OCL_VENDOR_AMD == ocl_gpu_device->vendor_e))
    {
        return egpuCompatible;
    }

    return egpuIncompatible;
}

/*! \brief Returns an ocl_vendor_id_t value corresponding to the input OpenCL vendor name.
 *
 *  \param[in] vendor_name String with OpenCL vendor name.
 *  \returns               ocl_vendor_id_t value for the input vendor_name
 */
ocl_vendor_id_t get_vendor_id(char *vendor_name)
{
    if (vendor_name)
    {
        if (strstr(vendor_name, "NVIDIA"))
        {
            return OCL_VENDOR_NVIDIA;
        }
        else
        if (strstr(vendor_name, "AMD") ||
            strstr(vendor_name, "Advanced Micro Devices"))
        {
            return OCL_VENDOR_AMD;
        }
        else
        if (strstr(vendor_name, "Intel"))
        {
            return OCL_VENDOR_INTEL;
        }
    }
    return OCL_VENDOR_UNKNOWN;
}


//! This function is documented in the header file
int detect_gpus(gmx_gpu_info_t *gpu_info, char *err_str)
{
    int             retval;
    cl_uint         ocl_platform_count;
    cl_platform_id *ocl_platform_ids;
    cl_device_type  req_dev_type = CL_DEVICE_TYPE_GPU;

    retval           = 0;
    ocl_platform_ids = NULL;

    if (getenv("GMX_OCL_FORCE_CPU") != NULL)
    {
        req_dev_type = CL_DEVICE_TYPE_CPU;
    }

    while (1)
    {
        CALLOCLFUNC_LOGERROR(clGetPlatformIDs(0, NULL, &ocl_platform_count), err_str, retval)
        if (0 != retval)
        {
            break;
        }

        if (1 > ocl_platform_count)
        {
            break;
        }

        snew(ocl_platform_ids, ocl_platform_count);

        CALLOCLFUNC_LOGERROR(clGetPlatformIDs(ocl_platform_count, ocl_platform_ids, NULL), err_str, retval)
        if (0 != retval)
        {
            break;
        }

        for (unsigned int i = 0; i < ocl_platform_count; i++)
        {
            cl_uint ocl_device_count;

            /* If requesting req_dev_type devices fails, just go to the next platform */
            if (CL_SUCCESS != clGetDeviceIDs(ocl_platform_ids[i], req_dev_type, 0, NULL, &ocl_device_count))
            {
                continue;
            }

            if (1 <= ocl_device_count)
            {
                gpu_info->n_dev += ocl_device_count;
            }
        }

        if (1 > gpu_info->n_dev)
        {
            break;
        }

        snew(gpu_info->gpu_dev, gpu_info->n_dev);

        {
            int           device_index;
            cl_device_id *ocl_device_ids;

            snew(ocl_device_ids, gpu_info->n_dev);
            device_index = 0;

            for (unsigned int i = 0; i < ocl_platform_count; i++)
            {
                cl_uint ocl_device_count;

                /* If requesting req_dev_type devices fails, just go to the next platform */
                if (CL_SUCCESS != clGetDeviceIDs(ocl_platform_ids[i], req_dev_type, gpu_info->n_dev, ocl_device_ids, &ocl_device_count))
                {
                    continue;
                }

                if (1 > ocl_device_count)
                {
                    break;
                }

                for (unsigned int j = 0; j < ocl_device_count; j++)
                {
                    gpu_info->gpu_dev[device_index].ocl_gpu_id.ocl_platform_id = ocl_platform_ids[i];
                    gpu_info->gpu_dev[device_index].ocl_gpu_id.ocl_device_id   = ocl_device_ids[j];

                    gpu_info->gpu_dev[device_index].device_name[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_NAME, sizeof(gpu_info->gpu_dev[device_index].device_name), gpu_info->gpu_dev[device_index].device_name, NULL);

                    gpu_info->gpu_dev[device_index].device_version[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_VERSION, sizeof(gpu_info->gpu_dev[device_index].device_version), gpu_info->gpu_dev[device_index].device_version, NULL);

                    gpu_info->gpu_dev[device_index].device_vendor[0] = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_VENDOR, sizeof(gpu_info->gpu_dev[device_index].device_vendor), gpu_info->gpu_dev[device_index].device_vendor, NULL);

                    gpu_info->gpu_dev[device_index].compute_units = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(gpu_info->gpu_dev[device_index].compute_units), &(gpu_info->gpu_dev[device_index].compute_units), NULL);

                    gpu_info->gpu_dev[device_index].adress_bits = 0;
                    clGetDeviceInfo(ocl_device_ids[j], CL_DEVICE_ADDRESS_BITS, sizeof(gpu_info->gpu_dev[device_index].adress_bits), &(gpu_info->gpu_dev[device_index].adress_bits), NULL);

                    gpu_info->gpu_dev[device_index].vendor_e = get_vendor_id(gpu_info->gpu_dev[device_index].device_vendor);

                    gpu_info->gpu_dev[device_index].stat = is_gmx_supported_gpu_id(gpu_info->gpu_dev + device_index);

                    if (egpuCompatible == gpu_info->gpu_dev[device_index].stat)
                    {
                        gpu_info->n_dev_compatible++;
                    }

                    device_index++;
                }
            }

            gpu_info->n_dev = device_index;

            /* Dummy sort of devices -  AMD first, then NVIDIA, then Intel */
            // TODO: Sort devices based on performance.
            if (0 < gpu_info->n_dev)
            {
                int last = -1;
                for (int i = 0; i < gpu_info->n_dev; i++)
                {
                    if (OCL_VENDOR_AMD == gpu_info->gpu_dev[i].vendor_e)
                    {
                        last++;

                        if (last < i)
                        {
                            gmx_device_info_t ocl_gpu_info;
                            ocl_gpu_info            = gpu_info->gpu_dev[i];
                            gpu_info->gpu_dev[i]    = gpu_info->gpu_dev[last];
                            gpu_info->gpu_dev[last] = ocl_gpu_info;
                        }
                    }
                }

                /* if more than 1 device left to be sorted */
                if ((gpu_info->n_dev - 1 - last) > 1)
                {
                    for (int i = 0; i < gpu_info->n_dev; i++)
                    {
                        if (OCL_VENDOR_NVIDIA == gpu_info->gpu_dev[i].vendor_e)
                        {
                            last++;

                            if (last < i)
                            {
                                gmx_device_info_t ocl_gpu_info;
                                ocl_gpu_info            = gpu_info->gpu_dev[i];
                                gpu_info->gpu_dev[i]    = gpu_info->gpu_dev[last];
                                gpu_info->gpu_dev[last] = ocl_gpu_info;
                            }
                        }
                    }
                }
            }

            sfree(ocl_device_ids);
        }

        break;
    }

    sfree(ocl_platform_ids);

    return retval;
}

//! This function is documented in the header file
void free_gpu_info(const gmx_gpu_info_t gmx_unused *gpu_info)
{
    if (gpu_info)
    {
        for (int i = 0; i < gpu_info->n_dev; i++)
        {
            cl_int gmx_unused cl_error;

            if (gpu_info->gpu_dev[i].context)
            {
                cl_error                     = clReleaseContext(gpu_info->gpu_dev[i].context);
                gpu_info->gpu_dev[i].context = NULL;
                assert(CL_SUCCESS == cl_error);
            }

            if (gpu_info->gpu_dev[i].program)
            {
                cl_error                     = clReleaseProgram(gpu_info->gpu_dev[i].program);
                gpu_info->gpu_dev[i].program = NULL;
                assert(CL_SUCCESS == cl_error);
            }
        }

        sfree(gpu_info->gpu_dev);
    }
}

//! This function is documented in the header file
void pick_compatible_gpus(const gmx_gpu_info_t *gpu_info,
                          gmx_gpu_opt_t        *gpu_opt)
{
    int  i, ncompat;
    int *compat;

    assert(gpu_info);
    /* gpu_dev/n_dev have to be either NULL/0 or not (NULL/0) */
    assert((gpu_info->n_dev != 0 ? 0 : 1) ^ (gpu_info->gpu_dev == NULL ? 0 : 1));

    snew(compat, gpu_info->n_dev);
    ncompat = 0;
    for (i = 0; i < gpu_info->n_dev; i++)
    {
        if (is_compatible_gpu(gpu_info->gpu_dev[i].stat))
        {
            ncompat++;
            compat[ncompat - 1] = i;
        }
    }

    gpu_opt->n_dev_compatible = ncompat;
    snew(gpu_opt->dev_compatible, ncompat);
    memcpy(gpu_opt->dev_compatible, compat, ncompat*sizeof(*compat));
    sfree(compat);
}

//! This function is documented in the header file
gmx_bool check_selected_gpus(int                  *checkres,
                             const gmx_gpu_info_t *gpu_info,
                             gmx_gpu_opt_t        *gpu_opt)
{
    int  i, id;
    bool bAllOk;

    assert(checkres);
    assert(gpu_info);
    assert(gpu_opt->n_dev_use >= 0);

    if (gpu_opt->n_dev_use == 0)
    {
        return TRUE;
    }

    assert(gpu_opt->dev_use);

    /* we will assume that all GPUs requested are valid IDs,
       otherwise we'll bail anyways */

    bAllOk = true;
    for (i = 0; i < gpu_opt->n_dev_use; i++)
    {
        id = gpu_opt->dev_use[i];

        /* devices are stored in increasing order of IDs in gpu_dev */
        gpu_opt->dev_use[i] = id;

        checkres[i] = (id >= gpu_info->n_dev) ?
            egpuNonexistent : gpu_info->gpu_dev[id].stat;

        bAllOk = bAllOk && is_compatible_gpu(checkres[i]);
    }

    return bAllOk;
}

//! This function is documented in the header file
void get_gpu_device_info_string(char gmx_unused *s, const gmx_gpu_info_t gmx_unused *gpu_info, int gmx_unused index)
{
    assert(s);
    assert(gpu_info);

    if (index < 0 && index >= gpu_info->n_dev)
    {
        return;
    }

    gmx_device_info_t  *dinfo = &gpu_info->gpu_dev[index];

    bool                bGpuExists =
        dinfo->stat == egpuCompatible ||
        dinfo->stat == egpuIncompatible;

    if (!bGpuExists)
    {
        sprintf(s, "#%d: %s, stat: %s",
                index, "N/A",
                gpu_detect_res_str[dinfo->stat]);
    }
    else
    {
        sprintf(s, "#%d: name: %s, vendor: %s, device version: %s, stat: %s",
                index, dinfo->device_name, dinfo->device_vendor,
                dinfo->device_version,
                gpu_detect_res_str[dinfo->stat]);
    }
}

//! This function is documented in the header file
gmx_bool init_gpu(FILE gmx_unused                 *fplog,
                  int                              mygpu,
                  char                            *result_str,
                  const gmx_gpu_info_t gmx_unused *gpu_info,
                  const gmx_gpu_opt_t             *gpu_opt
                  )
{
    assert(result_str);

    result_str[0] = 0;

    if (mygpu < 0 || mygpu >= gpu_opt->n_dev_use)
    {
        char        sbuf[STRLEN];
        sprintf(sbuf, "Trying to initialize an inexistent GPU: "
                "there are %d %s-selected GPU(s), but #%d was requested.",
                gpu_opt->n_dev_use, gpu_opt->bUserSet ? "user" : "auto", mygpu);
        gmx_incons(sbuf);
    }

    return TRUE;
}

//! This function is documented in the header file
int get_gpu_device_id(const gmx_gpu_info_t  *,
                      const gmx_gpu_opt_t  *gpu_opt,
                      int                   idx)
{
    assert(gpu_opt);
    assert(idx >= 0 && idx < gpu_opt->n_dev_use);

    return gpu_opt->dev_use[idx];
}

//! This function is documented in the header file
char* get_ocl_gpu_device_name(const gmx_gpu_info_t *gpu_info,
                              const gmx_gpu_opt_t  *gpu_opt,
                              int                   idx)
{
    assert(gpu_info);
    assert(gpu_opt);
    assert(idx >= 0 && idx < gpu_opt->n_dev_use);

    return gpu_info->gpu_dev[gpu_opt->dev_use[idx]].device_name;
}

//! This function is documented in the header file
size_t sizeof_gpu_dev_info(void)
{
    return sizeof(gmx_device_info_t);
}

/*! \brief Prints the name of a kernel function pointer.
 *
 * \param[in]    kernel   OpenCL kernel
 * \returns               CL_SUCCESS if the operation was successful, an OpenCL error otherwise.
 */
cl_int dbg_ocl_kernel_name(const cl_kernel kernel)
{
    cl_int cl_error;
    char   kernel_name[256];
    cl_error = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME,
                               sizeof(kernel_name), &kernel_name, NULL);
    if (cl_error)
    {
        printf("No kernel found!\n");
    }
    else
    {
        printf("%s\n", kernel_name);
    }
    return cl_error;
}

/*! \brief Prints the name of a kernel function pointer.
 *
 * \param[in]    kernel   OpenCL kernel
 * \returns               CL_SUCCESS if the operation was successful, an OpenCL error otherwise.
 */
cl_int dbg_ocl_kernel_name_address(void* kernel)
{
    cl_int cl_error;
    char   kernel_name[256];
    cl_error = clGetKernelInfo((cl_kernel)kernel, CL_KERNEL_FUNCTION_NAME,
                               sizeof(kernel_name), &kernel_name, NULL);
    if (cl_error)
    {
        printf("No kernel found!\n");
    }
    else
    {
        printf("%s\n", kernel_name);
    }
    return cl_error;
}

void gpu_set_host_malloc_and_free(bool               bUseGpuKernels,
                                  gmx_host_alloc_t **nb_alloc,
                                  gmx_host_free_t  **nb_free)
{
    if (bUseGpuKernels)
    {
        *nb_alloc = &ocl_pmalloc;
        *nb_free  = &ocl_pfree;
    }
    else
    {
        *nb_alloc = NULL;
        *nb_free  = NULL;
    }
}
