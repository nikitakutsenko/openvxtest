#include "../ref.h"
#include <stdio.h>

#define min_eps 1e-8
#define CV_PI   3.14159265

vx_uint32 num_of_solve = 0;

vx_status  transpose  (vx_matrix in,  vx_matrix out);
vx_status  mulm       (vx_matrix A,   vx_matrix B,   vx_matrix C);
vx_status  mulv       (vx_matrix A,   vx_array  B,   vx_array C);
vx_status  solve      (vx_matrix A,   vx_array  b,   vx_array x);


vx_status ref_FitEllipse(vx_array Xs, vx_array Ys, vx_array Out)
{
    if (Xs -> size != Ys -> size)
        return VX_ERROR_INVALID_PARAMETERS;

    vx_uint32 n = Xs -> size;

    if( n < 5 )
        return VX_ERROR_INVALID_PARAMETERS;

    vx_float32 x0 = 0, y0 = 0;
    vx_float32 *xs = Xs -> data, *ys = Ys -> data;

    for(vx_uint32 i = 0; i < n; i++ )
    {
        x0 += xs[i];
        y0 += ys[i];
    }

    x0 /= n;
    y0 /= n;

    struct _vx_matrix A = {
        calloc(n * 5, sizeof(vx_float64)),  // memory
        n,                                  // rows
        5,                                  // cols
        VX_TYPE_FLOAT64                     // type
    };

    struct _vx_array b = {
        calloc(n, sizeof(vx_float64)),      // memory
        n,                                  // size
        VX_TYPE_FLOAT64                     // type
    };

    struct _vx_array x = {
        calloc(5, sizeof(vx_float64)),      // memory
        5,                                  // size
        VX_TYPE_FLOAT64                     // type
    };

    vx_float64 *Ad = A.data, *bd = b.data, *xd = x.data;

    for (vx_uint32 i = 0; i < 5; i++)
        xd[i] = 0;

    for (vx_uint32 i = 0; i < n; i++) {
        vx_float64 xc = xs[i] - x0, yc = ys[i] - y0; 

        bd[i] = 10000.0; // 1.0?
        Ad[i*5]     = - xc * xc;
        Ad[i*5 + 1] = - yc * yc;
        Ad[i*5 + 2] = - xc * yc;
        Ad[i*5 + 3] = xc;
        Ad[i*5 + 4] = yc;
    }


    if (solve(&A, &b, &x) != VX_SUCCESS)
        return VX_ERROR_INVALID_PARAMETERS;


    A.height = 2;
    A.width  = 2;
    b.size   = 2;
    x.size   = 2;

    Ad[0] = 2 * xd[0];
    Ad[1] = Ad[2] = xd[2];
    Ad[3] = 2 * xd[1];
    bd[0] = xd[3];
    bd[1] = xd[4];

    if (solve(&A, &b, &x) != VX_SUCCESS)
        return VX_ERROR_INVALID_PARAMETERS;

    A.height = n;
    A.width  = 3;
    b.size   = n;
    x.size   = 3;

    for(vx_uint32 i = 0; i < n; i++ )
    {
        vx_float64 xc = xs[i] - x0, yc = ys[i] - y0; 

        bd[i] = 1.0;
        Ad[i * 3    ] = (xc - xd[0]) * (xc - xd[0]);
        Ad[i * 3 + 1] = (yc - xd[1]) * (yc - xd[1]);
        Ad[i * 3 + 2] = (xc - xd[0]) * (yc - xd[1]);
    }
    
    vx_float64 * xd_ = calloc(5, sizeof(vx_float64));
    vx_float32 * Outd = Out -> data;

    xd_[0] = xd[0];
    xd_[1] = xd[1];
    

    if (solve(&A, &b, &x) != VX_SUCCESS)
        return VX_ERROR_INVALID_PARAMETERS;

    vx_float64 t;
    xd_[4] = -0.5 * atan2(xd[2], xd[1] - xd[0]);

    if( fabs(xd[2]) > min_eps )
        t = xd[2] / sin(-2.0 * xd_[4]);
    else
        t = xd[1] - xd[0];

    xd_[2] = fabs(xd[0] + xd[1] - t);

    if( xd_[2] > min_eps )
        xd_[2] = sqrt(2.0 / xd_[2]);
    
    xd_[3] = fabs(xd[0] + xd[1] + t);

    if( xd_[3] > min_eps )
        xd_[3] = sqrt(2.0 / xd_[3]);

    Outd[0] = (vx_float32) xd_[0] + x0;
    Outd[1] = (vx_float32) xd_[1] + y0;

    Outd[2] = (vx_float32) (xd_[2]*2);
    Outd[3] = (vx_float32) (xd_[3]*2);

    if( Outd[2] > Outd[3] ) 
    {
        vx_float32 tmp = Outd[3];
        Outd[3] = Outd[2];
        Outd[2] = tmp;
        Outd[4] = (vx_float32)(90 + xd_[4]*180 / CV_PI);
    }

    if( Outd[4] < -180 )
        Outd[4] += 360;
    if( Outd[4] > 360 )
        Outd[4] -= 360;

    for (vx_uint32 i = 0; i < 5; i ++)
        if (Outd[i] != Outd[i]) 
            return 1;
    return VX_SUCCESS;
}


vx_status solve(vx_matrix A, vx_array b, vx_array x) {
    vx_uint32 n = A -> height; // колво строк
    vx_uint32 m = A ->  width; // колво столбцов

    vx_float64 * Ad = A -> data;
    vx_float64 * bd = b -> data;
    vx_float64 * xd = x -> data;

    if (n != b -> size || m != x -> size)
        return VX_ERROR_INVALID_PARAMETERS;

    if(n == m) {
        for (vx_uint32 i = 0; i < n; i++)
        {
            vx_float64 tmp = Ad[i * n + i];

            for (vx_int32 j = n - 1; j >= (vx_int32) i; j--){
                Ad[i * n + j] = Ad[i * n + j] / tmp;
            }

            bd[i] = bd[i] / tmp;

            for (vx_int32 j = i + 1; j < (vx_int32) n; j++)
            {
                tmp = Ad[j * n + i];

                for (vx_int32 k = n - 1; k >= (vx_int32)i; k--)
                    Ad[j * n + k] -= tmp * Ad[i * n + k];
                
                bd[j] -= tmp * bd[i];
            }
        }

        xd[n-1] = bd[n-1];
        for (vx_int32 i = n - 2; i >= 0; i--)
        {
            xd[i] = bd[i];
            for (vx_int32 j = i + 1; j < (vx_int32) n; j++) 
                xd[i] -= Ad[i * n + j] * xd[j];
        }

        return VX_SUCCESS;
    }
    
    if (n != m) {
        struct _vx_matrix _A = {
            calloc (n*m, sizeof(vx_float64)),
            m, n, VX_TYPE_FLOAT64 // m -строк n -столбцов
        };

        transpose(A, &_A); // транспонированная A

        vx_float64 * AD = _A.data;

        struct _vx_matrix ATA = {
            calloc (m*m, sizeof(vx_float64)),
            m, m, VX_TYPE_FLOAT64  // m x m
        };

        struct _vx_array  ATb = {
            calloc (m, sizeof(vx_float64)),
            m, VX_TYPE_FLOAT64     // m
        };

        mulm(&_A, A, &ATA);         // ATA = A' * A
        mulv(&_A, b, &ATb);         // ATb = A' * b

        solve(&ATA, &ATb, x);       // A * x = b  =>  A' * A * x = A' * b
    }

    return VX_SUCCESS;
}

vx_status transpose (vx_matrix in, vx_matrix out) {
    vx_float64 * outd = out -> data;
    vx_float64 * ind  =  in -> data;

    vx_uint32 n = in -> height; // колво строк или длинна столбца
    vx_uint32 m = in -> width;  // колво столбцов или длинна строки
    
    for (vx_uint32 i = 0; i < n; i++)
        for (vx_uint32 j = 0; j < m; j++)
            outd[j * n + i] = ind[i * m + j];

    return VX_SUCCESS;
}

vx_status mulm(vx_matrix A, vx_matrix B, vx_matrix C) {
    vx_float64 * Ad = A -> data;
    vx_float64 * Bd = B -> data;
    vx_float64 * Cd = C -> data;

    /*
        A -> n x m
        B -> m x l
        C -> n x l
    */

    vx_uint32 
        n = A -> height,
        m = B -> height,
        l = B -> width;

    for (vx_uint32 i = 0; i < n * l; i++)
        Cd[i] = 0;

    for (vx_uint32 i = 0; i < n; i++){
        for (vx_uint32 j = 0; j < l; j++){
            for (vx_uint32 k = 0; k < m; k++)
                Cd[i * l + j] += Ad[i * m + k] * Bd[k * l + j];
        }
    }


    return VX_SUCCESS;
}

vx_status mulv(vx_matrix A, vx_array B, vx_array C) {
    vx_float64 * Ad = A -> data;
    vx_float64 * Bd = B -> data;
    vx_float64 * Cd = C -> data;

    /*
        A -> n x m
        B -> m * 1
        C -> n * 1
    */

    vx_uint32
        n = A -> height, 
        m = B -> size;

    for (vx_uint32 i = 0; i < n; i++)
        Cd[i] = 0;

    
    for (vx_uint32 i = 0; i < n; i++){
            for (vx_uint32 k = 0; k < m; k++)
                Cd[i] += Ad[i * m + k] * Bd[k];
    }

    return VX_SUCCESS;
}
