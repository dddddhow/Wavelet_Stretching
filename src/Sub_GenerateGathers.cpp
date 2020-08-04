//----------------------------------------------------------------------------//
//程序目的： 1、生成CMP道集
//           2、生成CMP零偏移距道集
//           3、生成Theta-T道集
//           4、生成Theta-T零角度道集
//
//程序原理：
//
//程序参数说明：
//          nx：横向样点数（个）
//          dx：横向取样间隔（米）
//          nt：时间样点数（个）
//          dt：时间采样（秒）
//          nw：子波长度 （个）
//          v_rms[][] ：均方根速度（米/秒）
//          t0：分界面零偏移距时间（秒）
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 08 01
//----------------------------------------------------------------------------//

#include <armadillo>
#include <iostream>
#include "../lib/TJWPI_ShengShen.h"

using namespace std;
using namespace arma;


int generate_cmpgathers_func(arma::Mat<float> &seis_xt, arma::Col<float> &v_rms,
        arma::Col<float> &wavelet, int nw,
        int nx, int nt, float dx, float dt , float t0)
{
    //时间-空间域剖面（CMP）合成（基于褶积运算）
    //参数说明：
    //      ref_cof_xt :反射系数矩阵
    //      seis_xt    :合成地震记录（CMP）
    arma::Mat<float> ref_cof_xt(nt, nx, fill::zeros);

    for (int ix = 0; ix < nx; ix++)
    {
        float x_temp = (ix - nx / 2) * dx;
        float t0_temp = t0;
        float v=v_rms(0);
        float t = sqrt(t0_temp * t0_temp + 4.0 * x_temp * x_temp * 1.0 / v * 1.0 / v);
        int it = int(t / dt) + nw / 2;
        if (it < nt)
        {
            ref_cof_xt(it, ix) = 1;
        }
    }

    for (int ix = 0; ix < nx; ix++)
    {
        arma::Col<float> vec_temp = conv(ref_cof_xt.col(ix), wavelet);
        seis_xt.col(ix) = vec_temp(span(nw - 1, nt + nw - 2));
    }

    return 0;
}




int generate_thetatgathers_func(arma::Mat<float> &seis_thetat, arma::Col<float> &v_rms,
        arma::Col<float> &wavelet, int nw,
        int nx, int nt, float dx, float dt , float t0, float h)
{
    //时间-角度域剖面（T-Theta）合成（基于褶积运算）
    //参数说明：
    //      ref_cof_thetat :反射稀疏矩阵
    //      seis_thetat    :合成地震记录（CMP）
    float theta_max=atan(tan(80.0/180*3.1435926));
    float d_theta=theta_max*1.0/((nx-1)/2);

    int n_theta=nx;
    arma::Mat<float> ref_cof_thetat(nt,n_theta,fill::zeros);

    for(int ix=0; ix<n_theta; ix++)
    {
        float v=v_rms(0);
        float theta_max_here=atan((nx-1)*dx*1.0/2*1.0/h);
        float theta=(ix-n_theta/2)*d_theta;
        if(abs(theta)<theta_max_here)
        {
            //float t=sqrt(t0*t0+4.0*h*h*sin(theta)*sin(theta)
            //*1.0/v*1.0/v*1.0/cos(theta)*1.0/cos(theta));
            float t=t0*1.0/cos(theta);
            int it=int(t/dt)+nw/2;
            if(it<nt)
            {
                ref_cof_thetat(it,ix)=1;
            }
        }
    }

    for(int i_theta=0; i_theta<nx; i_theta++)
    {
        arma::Col<float> vec_temp=conv(ref_cof_thetat.col(i_theta),wavelet);
        seis_thetat.col(i_theta)=vec_temp(span(nw-1,nt-2+nw));
    }

return 0;
}




int generate_cmpgathers_modelbased_func(arma::Mat<float> &seis_xt, arma::Col<float> &v_rms,
        arma::Col<float> &model,
        arma::Col<float> &wavelet, int nw,
        float dz, int nvz,
        int nx, int nt, float dx, float dt)
{
        //时间-空间域剖面（CMP）合成（基于褶积运算）
        //参数说明：
        //      ref_cof_xt :反射系数矩阵
        //      seis_xt    :合成地震记录（CMP）

        arma::Mat<float> ref_cof_xt(nt, nx, fill::zeros);
        arma::Col<float> ref_location(nvz,fill::zeros);
        arma::Col<float> ref_location_t0(nvz,fill::zeros);
        int n_ref=0;

        for(int ivz=0; ivz<nvz-1; ivz++)
        {
            if(model(ivz+1) -model(ivz) > 100)
            {
                ref_location(n_ref)=(ivz+1)*dz;
                n_ref++;
            }
        }
        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ivz = 0; ivz < ref_location(ir) / dz; ivz++)
            {
                ref_location_t0(ir) = ref_location_t0(ir) + dz * 1.0 / model(ivz);
            }
        }

        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                float x_temp = (ix - nx / 2) * dx;
                float t0_temp = ref_location_t0(ir);
                int iz = ref_location(ir) / dz;
                float v = v_rms(iz);
                float t = sqrt(t0_temp * t0_temp + 4.0 * x_temp * x_temp * 1.0 / v * 1.0 / v);
                int it = int(t / dt) + nw / 2;
                if (it < nt)
                {
                    ref_cof_xt(it, ix) = 1;
                }
            }
        }
        for (int ix = 0; ix < nx; ix++)
        {
            arma::Col<float> vec_temp = conv(ref_cof_xt.col(ix), wavelet);
            seis_xt.col(ix) = vec_temp(span(nw - 1, nt + nw - 2));
        }


    return 0;
}


int generate_thetatgathers_modelbased_func(arma::Mat<float> &seis_thetat,
        arma::Col<float> &model,
        arma::Col<float> &wavelet, int nw,
        float dz, int nvz, float d_theta,
        int n_theta, int nt,float dt)
{
        //Theta-T合成（基于褶积运算）
        //参数说明：
        //      ref_cof_thetat :反射系数矩阵

        arma::Mat<float> ref_cof_thetat(nt, n_theta, fill::zeros);
        arma::Col<float> ref_location(nvz,fill::zeros);
        arma::Col<float> ref_location_t0(nvz,fill::zeros);
        int n_ref=0;

        for(int ivz=0; ivz<nvz-1; ivz++)
        {
            if(model(ivz+1) -model(ivz) > 100)
            {
                ref_location(n_ref)=(ivz+1)*dz;
                n_ref++;
            }
        }
        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ivz = 0; ivz < ref_location(ir) / dz; ivz++)
            {
                ref_location_t0(ir) = ref_location_t0(ir) + dz * 1.0 / model(ivz);
            }
        }

        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int i_theta = 0; i_theta < n_theta; i_theta++)
            {
                float theta_temp=(i_theta-n_theta/2)*d_theta;
                float t0_temp = ref_location_t0(ir);
                float t = t0_temp*1.0/cos(theta_temp);
                int it = int(t / dt) + nw / 2;
                if (it < nt)
                {
                    ref_cof_thetat(it, i_theta) = 1;
                }
            }
        }
        for (int i_theta = 0; i_theta < n_theta; i_theta++)
        {
            arma::Col<float> vec_temp = conv(ref_cof_thetat.col(i_theta), wavelet);
            seis_thetat.col(i_theta) = vec_temp(span(nw - 1, nt + nw - 2));
        }


    return 0;
}


int generate_cmpgathers_zerooffset_modelbased_func(arma::Mat<float> &seis_xt_zerooffset,
        arma::Col<float> &model,
        arma::Col<float> &wavelet, int nw,
        float dz, int nvz,
        int nx, int nt, float dt)
{
        //零偏移距时间-空间域剖面（CMP）合成（基于褶积运算）
        //参数说明：
        //      ref_cof_xt :反射系数矩阵
        //      seis_xt    :合成地震记录（CMP）

        arma::Mat<float> ref_cof_xt(nt, nx, fill::zeros);
        arma::Col<float> ref_location(nvz,fill::zeros);
        arma::Col<float> ref_location_t0(nvz,fill::zeros);
        int n_ref=0;

        for(int ivz=0; ivz<nvz-1; ivz++)
        {
            if(model(ivz+1) -model(ivz) > 100)
            {
                ref_location(n_ref)=(ivz+1)*dz;
                n_ref++;
            }
        }
        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ivz = 0; ivz < ref_location(ir) / dz; ivz++)
            {
                ref_location_t0(ir) = ref_location_t0(ir) + dz * 1.0 / model(ivz);
            }
        }

        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                float t0_temp = ref_location_t0(ir);
                float t = t0_temp;
                int it = int(t / dt) + nw / 2;
                if (it < nt)
                {
                    ref_cof_xt(it, ix) = 1;
                }
            }
        }
        for (int ix = 0; ix < nx; ix++)
        {
            arma::Col<float> vec_temp = conv(ref_cof_xt.col(ix), wavelet);
            seis_xt_zerooffset.col(ix) = vec_temp(span(nw - 1, nt + nw - 2));
        }


    return 0;
}



int generate_thetatgathers_zerotheta_modelbased_func(arma::Mat<float> &seis_thetat_zerotheta,
        arma::Col<float> &model,
        arma::Col<float> &wavelet, int nw,
        float dz, int nvz,
        int n_theta, int nt,float dt)
{
        //零角度Theta-T合成（基于褶积运算）
        //参数说明：
        //      ref_cof_thetat :反射系数矩阵

        arma::Mat<float> ref_cof_thetat(nt, n_theta, fill::zeros);
        arma::Col<float> ref_location(nvz,fill::zeros);
        arma::Col<float> ref_location_t0(nvz,fill::zeros);
        int n_ref=0;

        for(int ivz=0; ivz<nvz-1; ivz++)
        {
            if(model(ivz+1) -model(ivz) > 100)
            {
                ref_location(n_ref)=(ivz+1)*dz;
                n_ref++;
            }
        }
        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int ivz = 0; ivz < ref_location(ir) / dz; ivz++)
            {
                ref_location_t0(ir) = ref_location_t0(ir) + dz * 1.0 / model(ivz);
            }
        }

        for (int ir = 0; ir < n_ref; ir++)
        {
            for (int i_theta = 0; i_theta < n_theta; i_theta++)
            {
                float t0_temp = ref_location_t0(ir);
                float t = t0_temp;
                int it = int(t / dt) + nw / 2;
                if (it < nt)
                {
                    ref_cof_thetat(it, i_theta) = 1;
                }
            }
        }
        for (int i_theta = 0; i_theta < n_theta; i_theta++)
        {
            arma::Col<float> vec_temp = conv(ref_cof_thetat.col(i_theta), wavelet);
            seis_thetat_zerotheta.col(i_theta) = vec_temp(span(nw - 1, nt + nw - 2));
        }


    return 0;
}

