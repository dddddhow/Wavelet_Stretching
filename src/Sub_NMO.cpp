//----------------------------------------------------------------------------//
//程序目的： 1、X-T NMO
//           2、THeta-T NMO
//           3、X_T域数据在Theta_T域进行NMO并返回到XT域
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
    //Theta-T域NMO
    //参数说明：
    //      t_thetat_norm ：按照时距关系公式计算的每一点的时间
    //      t_thetat_nmo  ：动校正量
    //      seis_thetat_after_nmo：动校正后的剖面

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
                float t = sqrt(t0_temp * t0_temp + 4.0 * h_temp * h_temp
                        * sin(theta) * sin(theta) * 1.0 / v * 1.0 / v
                        * 1.0 / cos(theta) * 1.0 / cos(theta));

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
            float theta = (ix - (n_theta - 1) / 2) * d_theta;
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
                    seis_thetat_after_nmo(it, ix) = b * seis_thetat(it_ori, ix)
                        + a * seis_thetat(it_ori + 1, ix);
                }
            }
        }
    }

    return 0;

}




int xt_nmo_thetatbased_func( arma::Mat<float> &seis_xt,arma::Mat<float> &seis_xt_after_nmo,
        arma::Col<float> &v_rms, int nx, int nt, float dx, float dt)
{
    //基于Theta-T关系式的时间-空间域（CMP）动校正
    //参数说明：
    //      seis_xt           : 地震数据
    //      seis_xt_interp    : 插值后的地震数据
    //      seis_xt_after_nmo : 动校正后的剖面
    arma::Mat<float> t_xt_norm(nt, nx, fill::zeros);
    arma::Mat<float> t_xt_nmo(nt, nx, fill::zeros);
    arma::Mat<float> seis_xt_interp;

    //X_T域道集加密
    cout << "     CMP Matrix Interpolation" << endl;
    {
        float x_interp_times = 50;
        float t_interp_times = 1;
        arma::fvec x = regspace<fvec>(0, nx - 1);
        arma::fvec y = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / x_interp_times, nx - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(x, y, seis_xt, xi, yi, seis_xt_interp);
    }

    int nx_interp = seis_xt_interp.n_cols;
    int nt_interp = seis_xt_interp.n_rows;
    float dx_interp = nx * 1.0 / nx_interp * dx;
    float dt_interp = nt * 1.0 / nt_interp * dt;


    for (int ix = 0; ix < nx; ix++)
    {
        for (int it = 0; it < nt; it++)
        {
            //X-T域    ：  A1（t，x）
            //Theta-T域：  B1（theta，t’）    ：A1（t，x）在ThetaT域对应的点
            //Theta-T域：  B2（theta，t’’）   ：B1（theta，t’）在Theta-T域对应的NMO点
            //X-T域    ：  A2（x’，t’’’）     ：B2（theta，t’’）在XT域对应的点


            float v=v_rms(it);

            //cout << "       Step1: choose A1(t,x)" << endl;
            //Step1: 选择 A1（t，x）
            float a1_x  = (ix - nx / 2) * dx;
            float a1_t  = it * dt;
            float a1_t0 = it*dt;;

            //cout << "       Step2: Change to B1(theta,t')" << endl;
            //Step2：变换到B1（theta，t’)
            float b1_theta=atan(a1_x*1.0/(v*a1_t0*1.0/2));
            float b1_t=a1_t0*1.0/cos(b1_theta);

            //cout << "       Step3: Calculate B2(theta,t'')" << endl;
            //Step3：寻找 B2（theta，t’’）
            float b2_theta=b1_theta;
            float b2_t=b1_t*1.0/cos(b2_theta);
            //B2（theta，t’)对应的b2_t0
            float b2_t0=b2_t*cos(b2_theta);

            //cout << "       Step4: Reverse to A2(x',t''')" << endl;
            //Step4：反变换到A2（x’，t’’’）
            if(int(b2_t0*1.0/dt)<nt && int(b2_t0*1.0/dt)>=0)
            {
                v=v_rms(int(b2_t0*1.0/dt));
            }
            float a2_x=tan(b2_theta)*v*b2_t0*1.0/2;
            float a2_t=sqrt(b2_t0*b2_t0+4.0*a2_x*a2_x*1.0/v*1.0/v);

            //cout << "       Step5: Calculate ix & it of A2(x',t''')" << endl;
            //Step5：插值求出A2（x’，t’’’）的坐标(a2_it, a2_ix)
            float max_sup_t = a2_t * 1.0 / dt_interp;
            float p_t = (int(max_sup_t * 1000000) % 1000000) * 1.0 / 1000000;
            float q_t = 1 - p_t;

            float max_sup_x = a2_x * 1.0 / dx_interp +nx_interp/2;
            float p_x = (int(max_sup_x * 1000000) % 1000000) * 1.0 / 1000000;
            float q_x = 1 - p_x;

            int a2_it_interp = int(max_sup_t);
            int a2_ix_interp = int(max_sup_x);
            if (a2_it_interp + 1 < nt_interp && a2_it_interp>=0 &&
                    a2_ix_interp +1 <nx_interp && a2_ix_interp >= 0)
            {
                seis_xt_after_nmo(it, ix) =
                    p_t *
                    (p_x * seis_xt_interp(a2_it_interp,a2_ix_interp) +
                     q_x * seis_xt_interp(a2_it_interp, a2_ix_interp + 1)) +
                    q_t *
                    (p_x * seis_xt_interp(a2_it_interp + 1, a2_ix_interp) +
                     q_x * seis_xt_interp(a2_it_interp + 1, a2_ix_interp + 1));
            }
        }
    }

    return 0;
}



int thetat_nmo_xtbased_func( arma::Mat<float> &seis_xt,arma::Mat<float> &seis_xt_after_nmo,
        arma::Col<float> &v_rms, int nx, int n_theta, int nt, float dx, float dt, float d_theta)
{
    //基于Theta-T关系式的时间-空间域（CMP）动校正
    //参数说明：
    //      seis_xt           : 地震数据
    //      seis_xt_interp    : 插值后的地震数据
    //      seis_xt_after_nmo : 动校正后的剖面
    arma::Mat<float> t_xt_norm(nt, nx, fill::zeros);
    arma::Mat<float> t_xt_nmo(nt, nx, fill::zeros);
    arma::Mat<float> seis_xt_interp;
    arma::Mat<float> seis_xt_after_nmo_interp;

    //X_T域道集加密
    cout << "     CMP Matrix Interpolation" << endl;
    {
        float x_interp_times = 50;
        float t_interp_times = 10;

        arma::fvec x  = regspace<fvec>(0, nx - 1);
        arma::fvec y  = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / x_interp_times, nx - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(x, y, seis_xt, xi, yi, seis_xt_interp);
        interp2(x, y, seis_xt_after_nmo, xi, yi, seis_xt_after_nmo_interp);
    }

    int nx_interp   = seis_xt_interp.n_cols;
    int nt_interp   = seis_xt_interp.n_rows;
    float dx_interp = nx * 1.0 / nx_interp * dx;
    float dt_interp = nt * 1.0 / nt_interp * dt;

    for (int i_theta = 0; i_theta < n_theta; i_theta++)
    {
        for (int it = 0; it < nt; it++)
        {
            //Theta-T域：  B1（t,theta）
            //Theta-T域：  B2（t',theta）    ：B1（theta，t）在ThetaT域NMO对应的点
            //X-T域    ：  A1（x,t'）        ：B1（t，theta）在XT域对应的点
            //X-T域    ：  A2（x',t''）      ：B2（t',theta）在XT域对应的点

            float v_a1 = v_rms(it);
            float v_a2 = v_rms(it);

            //Step1: 选择B1（t,theta）
            //cout << "       Step1: choose B1（t,theta）" << endl;
            float b1_theta = (i_theta - n_theta / 2) * d_theta;
            float b1_t0    = it*dt;
            float b1_t     = b1_t0*1.0/cos(b1_theta);

            //Step2: 计算B2（t',theta）
             //cout << "       Step2: Calculate B2（t',theta）" << endl;
             float b2_theta=b1_theta;
             float b2_t=b1_t;
             float b2_t0=b1_t0;
             if(int(b2_t/dt)<nt && int(b2_t/dt)>=0)
             {
             v_a2=v_rms(int(b2_t0*1.0/dt));
             }

            //Step3: B1（t,theta）变换到 A1（x,t'）
            //cout << "       Step3:B1（t,theta） Change to  A1（x,t'）" << endl;
            float a1_x=v_a1*b1_t0*1.0/2*tan(b1_theta);
            float a1_t=b1_t0;

            //Step4:  B2（t',theta）变换到  A2（x',t''）
            //cout << "       Step4: B2（t',theta）Change to  A2（x',t''）" << endl;
            float a2_x= v_a2*b2_t0*1.0/2*tan(b2_theta);
            float a2_t=b2_t;


            //Step5：插值求出 A1（x,t'） A2（x',t''）的坐标
            //cout << "       Step5: Calculate ix & it of A1（x,t'）& A2（x',t''）" << endl;
            float a1_max_sup_t = a1_t * 1.0 / dt_interp;
            float a1_p_t       = (int(a1_max_sup_t * 1000000) % 1000000) * 1.0 / 1000000;
            float a1_q_t       = 1 - a1_p_t;
            float a1_max_sup_x = a1_x * 1.0 / dx_interp +nx_interp/2;
            float a1_p_x       = (int(a1_max_sup_x * 1000000) % 1000000) * 1.0 / 1000000;
            float a1_q_x       = 1 - a1_p_x;
            float a2_max_sup_t = a2_t * 1.0 / dt_interp;
            float a2_p_t       = (int(a2_max_sup_t * 1000000) % 1000000) * 1.0 / 1000000;
            float a2_q_t       = 1 - a2_p_t;
            float a2_max_sup_x = a2_x * 1.0 / dx_interp +nx_interp/2;
            float a2_p_x       = (int(a2_max_sup_x * 1000000) % 1000000) * 1.0 / 1000000;
            float a2_q_x       = 1 - a2_p_x;

            int a1_it_interp = int(a1_max_sup_t);
            int a1_ix_interp = int(a1_max_sup_x);
            int a2_it_interp = int(a2_max_sup_t);
            int a2_ix_interp = int(a2_max_sup_x);

            if (a2_it_interp + 1 < nt_interp && a2_it_interp>=0 &&
                    a2_ix_interp +1 <nx_interp && a2_ix_interp >= 0 &&
               a1_it_interp + 1 < nt_interp && a1_it_interp>=0 &&
                    a1_ix_interp +1 <nx_interp && a1_ix_interp >= 0)
            {
                seis_xt_after_nmo_interp(a1_it_interp, a1_ix_interp) =
                    a2_p_t *
                    (a2_p_x * seis_xt_interp(a2_it_interp,a2_ix_interp) +
                     a2_q_x * seis_xt_interp(a2_it_interp, a2_ix_interp + 1)) +
                    a2_q_t *
                    (a2_p_x * seis_xt_interp(a2_it_interp + 1, a2_ix_interp) +
                     a2_q_x * seis_xt_interp(a2_it_interp + 1, a2_ix_interp + 1));
            }
        }
    }


    cout << "     CMP Matrix Inv-Interpolation" << endl;
    {
        float x_interp_times = 50;
        float t_interp_times = 10;

        arma::fvec x  = regspace<fvec>(0, nx - 1);
        arma::fvec y  = regspace<fvec>(0, nt - 1);

        arma::fvec xi = regspace<fvec>(0, 1.0 / x_interp_times, nx - 1);
        arma::fvec yi = regspace<fvec>(0, 1.0 / t_interp_times, nt - 1);

        interp2(xi, yi, seis_xt_after_nmo_interp, x, y, seis_xt_after_nmo);
    }


    return 0;
}



int xt_nmo_pilot_trace_based_func(arma::Mat<float> &seis_xt, arma::Mat<float> &seis_xt_after_nmo,
        arma::Col<float> &v_rms, int nx, int nt, float dx, float dt, int pilottrace)
{
    //时间-空间域（CMP）动校正
    //参数说明：
    //      t_xt_norm ：按照时距关系公式计算的每一点的时间
    //      t_xt_nmo  ：动校正量
    //      seis_xt_after_nmo：动校正后的剖面
    //      seis_xt_op2：动校正后域零偏移距剖面做差结果
    arma::Mat<float> t_xt_norm(nt, nx, fill::zeros);
    arma::Mat<float> t_xt_nmo(nt, nx, fill::zeros);

    for (int it = 0; it < nt; it++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            float t0_temp=it*dt;
            float v=v_rms(t0_temp*1.0/dt);

            float a1_x=(ix-nx/2)*dx;
            float a1_t=sqrt(t0_temp*t0_temp+4.0*a1_x*a1_x*1.0/v*1.0/v);

            float x_pilottrace=(pilottrace-nx/2)*dx;
            float t_pilottrace=sqrt(t0_temp*t0_temp+4.0*x_pilottrace*x_pilottrace*1.0/v*1.0/v);

            float max_sup = a1_t * 1.0 / dt;
            float a = (int(max_sup * 1000000) % 1000000) * 1.0 / 1000000;
            float b = 1 - a;

            int it_ori = int(max_sup);
            int it_pilottrace=int(t_pilottrace*1.0/dt);

            if (it_ori + 1 < nt && it_ori>=0
                    && it_pilottrace <nt && it_pilottrace >=0)
            {
                seis_xt_after_nmo(it_pilottrace, ix) = b * seis_xt(it_ori, ix) +
                    a * seis_xt(it_ori + 1, ix);
            }

        }
    }

    return 0;
}
