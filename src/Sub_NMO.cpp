//----------------------------------------------------------------------------//
//程序目的： 1、X-T NMO
//           2、THeta-T NMO
//
//程序原理：
//
//
//程序参数说明：
//          t_xt_norm ：按照时距关系公式计算的每一点的时间
//          t_xt_nmo  ：动校正量
//          seis_xt_after_nmo：动校正后的剖面
//          seis_xt[][] : X_T域数据
//          nx：横向样点数（个）
//          dx：横向取样间隔（米）
//          nt：时间样点数（个）
//          dt：时间采样（秒）
//          v_rms[][] ：均方根速度（米/秒）
//          t_thetat_norm ：按照时距关系公式计算的每一点的时间
//          t_thetat_nmo  ：动校正量
//          seis_thetat_after_nmo：动校正后的剖面
//
//程序依赖关系说明：
//          1、需要armadillo库
//          2、-std=c++11 或者 更高
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 08 01
//----------------------------------------------------------------------------//

#include <iostream>
#include <armadillo>
#include "../lib/TJWPI_ShengShen.h"

using namespace std;
using namespace arma;


int xt_nmo_func( arma::Mat<float> &seis_xt,arma::Mat<float> &seis_xt_after_nmo,
        arma::Col<float> &v_rms, int nx, int nt, float dx, float dt)
{
    //时间-空间域（CMP）动校正
    //参数说明：
    //      t_xt_norm ：按照时距关系公式计算的每一点的时间
    //      t_xt_nmo  ：动校正量
    //      seis_xt_after_nmo：动校正后的剖面
    //      seis_xt_op2：动校正后域零偏移距剖面做差结果
    arma::Mat<float> t_xt_norm(nt, nx, fill::zeros);
    arma::Mat<float> t_xt_nmo(nt, nx, fill::zeros);

    for (int ix = 0; ix < nx; ix++)
    {
        for (int it = 0; it < nt; it++)
        {
            float t0_temp = it * dt;
            float x_temp = abs((ix - (nx - 1) / 2) * dx);
            float v=v_rms(it);
            float t = sqrt(t0_temp * t0_temp + 4.0 * x_temp * x_temp * 1.0 / v * 1.0 / v);
            t_xt_norm(it, ix) = t;
            t_xt_nmo(it, ix) = t_xt_norm(it, ix) - t0_temp;
        }
    }

    for (int ix = 0; ix < nx; ix++)
    {
        for (int it = 0; it < nt; it++)
        {
            float t0_temp = it * dt;
            float t = t0_temp + t_xt_nmo(it, ix);
            float max_sup = t * 1.0 / dt;
            float a = (int(max_sup * 1000000) % 1000000) * 1.0 / 1000000;
            float b = 1 - a;

            int it_ori = int(max_sup);
            if (it_ori + 1 < nt)
            {
                seis_xt_after_nmo(it, ix) = b * seis_xt(it_ori, ix) +
                    a * seis_xt(it_ori + 1, ix);
            }
        }
    }

    return 0;
}





int thetat_nmo_func( arma::Mat<float> &seis_thetat,arma::Mat<float> &seis_thetat_after_nmo,
        arma::Col<float> &v_rms,
        int nx, int n_theta, int nt,
        float dx, float d_theta, float dt)
{

    arma::Mat<float> t_thetat_norm(nt, n_theta, fill::zeros);
    arma::Mat<float> t_thetat_nmo(nt, n_theta, fill::zeros);

    for (int ix = 0; ix < n_theta; ix++)
    {
        for (int it = 1; it < nt; it++)
        {
            float theta = (ix - (nx - 1) / 2) * d_theta;
            float v=v_rms(it);
            float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / (it * dt * v / 2));
            if (abs(theta) < theta_max_here)
            {
                float h_temp = v * it * dt / 2;
                float t0_temp = it * dt;
                float t = sqrt(t0_temp * t0_temp + 4.0 * h_temp * h_temp * sin(theta) * sin(theta) * 1.0 / v * 1.0 / v * 1.0 / cos(theta) * 1.0 / cos(theta));

                if (int(t * 1.0 / dt) < nt)
                {
                    t_thetat_norm(it, ix) = t;
                    t_thetat_nmo(it, ix) = t - t0_temp;
                }
            }
        }
    }

    for (int ix = 0; ix < n_theta; ix++)
    {
        for (int it = 0; it < nt; it++)
        {
            float v=v_rms(it);
            float theta = (ix - (nx - 1) / 2) * d_theta;
            float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / (it * dt * v / 2));
            if (abs(theta) < theta_max_here)
            {
                float t0_temp = it * dt;
                float t = t0_temp + t_thetat_nmo(it, ix);
                float max_sup = t * 1.0 / dt;
                float a = (int(max_sup * 1000000) % 1000000) * 1.0 / 1000000;
                float b = 1 - a;

                int it_ori = int(max_sup);
                //cout<<it_ori<<endl;
                if (it_ori + 1 < nt)
                {
                    seis_thetat_after_nmo(it, ix) = b * seis_thetat(it_ori, ix) + a * seis_thetat(it_ori + 1, ix);
                }
            }
        }
    }

    return 0;

}





