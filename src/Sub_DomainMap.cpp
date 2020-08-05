//----------------------------------------------------------------------------//
//程序目的： 1、X_T域数据 映射到 Theta_T域
//           2、Theta_T域数据 映射到 X_T域
//           3、Theta_T NMO之后数据 映射到 X_T域
//
//程序原理：
//
//
//程序参数说明：
//          seis_xt[][] : X_T域数据
//          seis_thetat[][]: Theta_T域数据
//          nx：横向样点数（个）
//          dx：横向取样间隔（米）
//          nt：时间样点数（个）
//          dt：时间采样（秒）
//          v_rms[][] ：均方根速度（米/秒）
//
//程序依赖关系说明：
//          1、需要armadillo库
//          2、-std=c++11 或者 更高
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 08 05
//----------------------------------------------------------------------------//

#include <iostream>
#include <armadillo>
#include "../lib/TJWPI_ShengShen.h"

using std::cout;
using namespace arma;

int xt_to_thetat_func(arma::Mat<float> &seis_xt, arma::Mat<float> &seis_thetat,
        arma::Col<float> &v_rms,
        int nx,  int nt, int dx, int dz, float dt)
{

    // CMP Domain-> Theta-T Domain
    //时间-角度域剖面（T-Theta）合成（基于褶积运算）
    //参数说明：
    //      ref_cof_thetat :反射稀疏矩阵
    //      seis_thetat    :合成地震记录（CMP）
    int n_theta = nx;
    float theta_max = atan(tan(80*3.1415926/180));
    float d_theta = theta_max * 1.0 / ((n_theta - 1) / 2);
    cout << "         Theta max is      : " << theta_max << " radians";
    cout << " : " << theta_max * 180 * 1.0 / 3.1415926 << " degrees" << endl;
    cout << "         D_Theta  is       : " << d_theta << " radinans";
    cout << " : " << d_theta * 180 * 1.0 / 3.1415926 << " degrees" << endl;

    //CMP道集加密
    cout << "     CMP Matrix Interpolation" << endl;
    arma::Mat<float> seis_xt_interp;
    {
        float x_interp_times = 50;
        float t_interp_times = 1;
        arma::fvec x = regspace<fvec>(0, nx - 1);
        arma::fvec y = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / x_interp_times, nx - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(x, y, seis_xt, xi, yi, seis_xt_interp);
    }
    //seis_xt_interp.save("../file_DataBased/seis_x_t_interp.dat",raw_binary);

    int nx_interp = seis_xt_interp.n_cols;
    int nt_interp = seis_xt_interp.n_rows;
    float dx_interp = nx * 1.0 / nx_interp * dx;
    float dt_interp = nt * 1.0 / nt_interp * dt;

    cout << "         The NX of Interpolated CMP is :" << nx_interp << endl;
    cout << "         The dx of Interpolated CMP is :" << dx_interp << " m" << endl;
    cout << "         The NT of Interpolated CMP is :" << nt_interp << endl;
    cout << "         The dt of Interpolated CMP is :" << dt_interp << " s" << endl;

    cout << "     Here We Transform CMP Gathers to Theta-T Gathers" << endl;
    for (int it = 0; it < nt; it++)
    {
        for (int i_theta = 0; i_theta < n_theta; i_theta++)
        {
            float theta = (i_theta - n_theta / 2) * d_theta;
            float theta_max_here = atan((nx / 2) * dx * 1.0 / (it * dt * v_rms(it) / 2));
            if (abs(theta) < theta_max_here)
            {
                //给定 theta t0，确定x
                float t0_temp = it * dt;
                float x = v_rms(it)* 1.0 / 2 * t0_temp * tan(theta);
                float x_coordinate_float = x * 1.0 / dx_interp;
                int x_coordinate = int(x_coordinate_float) + nx_interp / 2;

                //给定 Theta t0，确定t
                float t = t0_temp * 1.0 / cos(theta);
                float t_coordinate_float_interp = t * 1.0 / dt_interp;
                int t_coordinate_interp = int(t_coordinate_float_interp);

                float t_coordinate_float = t * 1.0 / dt;
                int t_coordinate = int(t_coordinate_float);
                if (x_coordinate + 1 < nx_interp && x_coordinate >= 0 &&
                        t_coordinate_interp + 1 < nt_interp && t_coordinate_interp >= 0 &&
                        t_coordinate + 1 < nt && t_coordinate >= 0)
                {
                    float a_x = (int(x_coordinate_float * 1000000) % 1000000) * 1.0 / 1000000;
                    float b_x = 1 - a_x;

                    float a_t = (int(t_coordinate_float_interp * 1000000) % 1000000) * 1.0 / 1000000;
                    float b_t = 1 - a_t;

                    seis_thetat(t_coordinate, i_theta) =
                        a_t *
                        (a_x * seis_xt_interp(t_coordinate_interp, x_coordinate) +
                         b_x * seis_xt_interp(t_coordinate_interp, x_coordinate + 1)) +
                        b_t *
                        (a_x * seis_xt_interp(t_coordinate_interp + 1, x_coordinate) +
                         b_x * seis_xt_interp(t_coordinate_interp + 1, x_coordinate + 1));
                }
            }
        }
    }

}



int thetat_to_xt_func(arma::Mat<float> &seis_thetat, arma::Mat<float> &seis_xt,
        arma::Col<float> &v_rms,
        int nx, int n_theta, int nt, int dx, float d_theta, int dz, float dt)
{
    //Theta-T域数据转换为X-T域

    //Theta-T道集加密
    cout << "       Theta-T Matrix Interpolation" << endl;
    arma::Mat<float> seis_thetat_interp;
    {
        float theta_interp_times = 50;
        float t_interp_times = 1;
        arma::fvec x = regspace<fvec>(0, n_theta - 1);
        arma::fvec y = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / theta_interp_times, n_theta - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(x, y, seis_thetat, xi, yi, seis_thetat_interp);
    }

    int ntheta_interp = seis_thetat_interp.n_cols;
    int nt_interp = seis_thetat_interp.n_rows;
    float dtheta_interp = nx * 1.0 / ntheta_interp * d_theta;
    float dt_interp = nt * 1.0 / nt_interp * dt;

    cout << "         The Nthetaof Interpolated Theta-T is :" << ntheta_interp << endl;
    cout << "         The dtheta of Interpolated  Theta-T is :" << dtheta_interp << " m" << endl;
    cout << "         The NT of Interpolated Theta-T is :" << nt_interp << endl;
    cout << "         The dt of Interpolated Theta-T is :" << dt_interp << " s" << endl;

    cout << "     Here We Transform Theta-T Gathers to CMP Gathers" << endl;

    for (int it = 0; it < nt; it++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            float x_temp = (ix - nx / 2) * dx;
            if (1)
            {
                //给定 x t0，确定theta
                float t0_temp = it * dt;
                float v=v_rms(it);
                float theta_temp = atan(x_temp * 1.0 / (v * t0_temp/2));
                float theta_coordinate_float = theta_temp * 1.0 / dtheta_interp;
                int theta_coordinate = int(theta_coordinate_float + ntheta_interp / 2);

                //给定 x t0，确定t
                // float t = sqrt(t0_temp * t0_temp + 4.0 * x_temp * x_temp * 1.0 / v * 1.0 / v);
                float t = t0_temp * 1.0 / cos(theta_temp);
                float t_coordinate_float_interp = t * 1.0 / dt_interp;
                int t_coordinate_interp = int(t_coordinate_float_interp);

                float t_coordinate_float = t * 1.0 / dt;
                int t_coordinate = int(t_coordinate_float);

                if (theta_coordinate + 1 < ntheta_interp && theta_coordinate >= 0 &&
                        t_coordinate_interp + 1 < nt_interp && t_coordinate_interp >= 0 &&
                        t_coordinate + 1 < nt && t_coordinate >= 0)
                {
                    float a_x = (int(theta_coordinate_float * 1000000) % 1000000) * 1.0 / 1000000;
                    float b_x = 1 - a_x;

                    float a_t = (int(t_coordinate_float_interp * 1000000) % 1000000) * 1.0 / 1000000;
                    float b_t = 1 - a_t;

                    seis_xt(t_coordinate, ix) =
                        a_t *
                        (a_x * seis_thetat_interp(t_coordinate_interp, theta_coordinate) +
                         b_x * seis_thetat_interp(t_coordinate_interp, theta_coordinate + 1)) +
                        b_t *
                        (a_x * seis_thetat_interp(t_coordinate_interp + 1, theta_coordinate) +
                         b_x * seis_thetat_interp(t_coordinate_interp + 1, theta_coordinate + 1));
                }
            }
        }
    }


    return 0;


}




int thetat_afternmo_to_xt_func(arma::Mat<float> &seis_thetat_afternmo, arma::Mat<float> &seis_xt,
        arma::Col<float> &v_rms,
        int nx, int n_theta, int nt, int dx, float d_theta, int dz, float dt)
{
    //Theta-T域数据转换为X-T域

    //Theta-T道集加密
    cout << "       Theta-T Matrix Interpolation" << endl;
    arma::Mat<float> seis_thetat_interp;
    {
        float theta_interp_times = 50;
        float t_interp_times = 1;
        arma::fvec x = regspace<fvec>(0, n_theta - 1);
        arma::fvec y = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / theta_interp_times, n_theta - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(x, y, seis_thetat_afternmo, xi, yi, seis_thetat_interp);
    }

    int ntheta_interp = seis_thetat_interp.n_cols;
    int nt_interp = seis_thetat_interp.n_rows;
    float dtheta_interp = nx * 1.0 / ntheta_interp * d_theta;
    float dt_interp = nt * 1.0 / nt_interp * dt;

    cout << "         The Nthetaof Interpolated Theta-T is :" << ntheta_interp << endl;
    cout << "         The dtheta of Interpolated  Theta-T is :" << dtheta_interp << " m" << endl;
    cout << "         The NT of Interpolated Theta-T is :" << nt_interp << endl;
    cout << "         The dt of Interpolated Theta-T is :" << dt_interp << " s" << endl;

    cout << "     Here We Transform Theta-T Gathers to CMP Gathers" << endl;

    for (int it = 0; it < nt; it++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            float x_temp = (ix - nx / 2) * dx;
            //给定 x t0，确定theta
            float t0_temp = it * dt;
            float v=v_rms(it);
            float theta_temp = atan(x_temp * 1.0 / (v * t0_temp/2));
            float theta_coordinate_float = theta_temp * 1.0 / dtheta_interp;
            int theta_coordinate = int(theta_coordinate_float + ntheta_interp / 2);

            //给定 x t0，确定t
            float t = t0_temp * 1.0 / cos(theta_temp);
            float t_coordinate_float_interp = t * 1.0 / dt_interp;
            int t_coordinate_interp = int(t_coordinate_float_interp);

            float t_coordinate_float = t * 1.0 / dt;
            int t_coordinate = int(t_coordinate_float);

            if (theta_coordinate + 1 < ntheta_interp && theta_coordinate >= 0 &&
                    t_coordinate_interp + 1 < nt_interp && t_coordinate_interp >= 0 &&
                    t_coordinate + 1 < nt && t_coordinate >= 0)
            {
                float a_x = (int(theta_coordinate_float * 1000000) % 1000000) * 1.0 / 1000000;
                float b_x = 1 - a_x;

                float a_t = (int(t_coordinate_float_interp * 1000000) % 1000000) * 1.0 / 1000000;
                float b_t = 1 - a_t;

                seis_xt(t_coordinate, ix) =
                    a_t *
                    (a_x * seis_thetat_interp(it, theta_coordinate) +
                     b_x * seis_thetat_interp(it, theta_coordinate + 1)) +
                    b_t *
                    (a_x * seis_thetat_interp(it + 1, theta_coordinate) +
                     b_x * seis_thetat_interp(it + 1, theta_coordinate + 1));
            }
        }
    }


    return 0;


}




