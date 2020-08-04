//----------------------------------------------------------------------------//
//程序目的： 1、X_T域数据 映射到 Theta_T域
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
//Time    ：2020 08 01
//----------------------------------------------------------------------------//

#include <iostream>
#include <armadillo>

using std::cout;
using namespace arma;

void xt_to_theta_func(arma::Mat<float> &seis_xt, arma::Mat<float> &seis_thetat,
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
