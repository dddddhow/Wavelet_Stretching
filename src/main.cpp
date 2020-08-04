//----------------------------------------------------------------------------//
//程序目的： 1、分析 时间-偏移距 空间NMO带来的子波拉伸问题
//           2、分析 时间-角度 空间NMO子波拉伸问题
//程序原理：
//          拿到CMP道集后，直接转变为Theta-T域（基于平层假设），然后分别在CMP域和
//          Theta-T域进行速度谱求取，分别自动拾取速度（均方根速度），根据均方根速度
//          进行叠加成像。
//
//程序参数说明：
//          nx：横向样点数（个）
//          dx：横向取样间隔（米）
//          nt：时间样点数（个）
//          dt：时间采样（秒）
//          nw：子波长度 （个）
//          v ：介质速度（米/秒）
//          nh：反射面深度（个）
//          h ：2倍反射面深度（米）
//          t0：分界面零偏移距时间（秒）
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 07 26
//----------------------------------------------------------------------------//

#include <armadillo>
#include <cmath>
#include <ctime>
#include <iostream>
#include "../lib/TJWPI_ShengShen.h"

using namespace std;
using namespace arma;

int main()
{
    cout << "--------------------------------------------------" << endl;
    cout << "START" << endl;
    clock_t start, finish;
    start = clock();
    cout << "--------------------------------------------------" << endl;

    //=========================================================================//
    //                               定义速度模型                              //
    //=========================================================================//
    cout << "--------------------------------------------------" << endl;
    //  model[][]:速度模型（米/秒）
    //  v_rms[][]:均方根速度模型（米/秒）
    //  nvx:模型横向距离（点）
    //  nvz:模型纵向距离（点）
    //  nx:单炮接收道数（双边接收）
    //  nt:单炮接收时间（点）
    //  nw:子波长度（点）
    //  dx:横向采样间距（米）
    //  dz:纵向采样间距（米）
    //  dt:时间采样间距（秒）

    //文件保存路径 fn
    string fn = "../file/";

    int nvx=200;
    int nvz=400;

    int nt = 2000;
    float dt = 0.001;

    int nw = 81;
    float v = 2000;

    int h = 200;
    float t0 = 2.0 * h * 1.0 / v;
    int nh = t0 * 1.0 / dt;

    int nx = 201;
    float dx = 5.0;
    float dz = 5.0;
    int nz = int(nt * dt * v / dz / 2);


    //速度模型加载
    //      速度模型model
    //      均方根速度v_rms
    //      均匀速度v_mean
    arma::Mat<float> model(nvz, nvx, fill::zeros);
    arma::Mat<float> v_rms(nt, nvx, fill::zeros);
    arma::Mat<float> v_mean(nt, nvx, fill::zeros);
    model.load("../data/model_nz400_nx200.dat", raw_binary);
    v_rms.load("../data/vrms_nt2000_nx200.dat", raw_binary);
    v_mean.load("../data/vmean_nt2000_nx200.dat", raw_binary);
    v_rms.reshape(nt, nvx);
    v_mean.reshape(nt, nvx);
    model.reshape(nvz, nvx);

    //由速度模型搜索反射点位置和个数
    //ref_num:      每一道有多少反射界面
    //ref_location: 每一道的每一个反射界面的深度（米）
    //ref_location_t0: 每一道的每一个反射界面的到达时（s）
    arma::Col<float> ref_num(nvx, fill::zeros);
    arma::Mat<float> ref_location(nvz, nvx, fill::zeros);
    arma::Mat<float> ref_location_t0(nvz, nvx, fill::zeros);

    for (int ivx = 0; ivx < nvx; ivx++)
    {
        int n_ref = 0;
        for (int ivz = 0; ivz < nvz - 1; ivz++)
        {
            if (model(ivz + 1, ivx) - model(ivz, ivx) > 100)
            {
                ref_location(n_ref, ivx) = (ivz + 1) * dz;
                n_ref++;
            }
        }
        ref_num(ivx) = n_ref;
    }

    for (int ivx = 0; ivx < nvx; ivx++)
    {
        for (int ir = 0; ir < ref_num(ivx); ir++)
        {
            for (int ivz = 0; ivz < ref_location(ir, ivx) / dz; ivz++)
            {
                ref_location_t0(ir, ivx) = ref_location_t0(ir, ivx) + dz * 1.0 / model(ivz, ivx);
            }
        }
    }

    //arma::Mat<float> v_rms(nt,nvx,fill::zeros);
    //for(int ivx=0; ivx<nvx; ivx++)
    //{
    //for(int it=0 ;it<nt; it++)
    //{
    //v_rms(it,ivx)=v;
    //}
    //}
    arma::Col<float> v_rms_col=v_rms.col(0);


    if (nh < nw / 2)
    {
        cout << "Error ! The Depth Minimum is " << (nw / 2) * dz << " m" << endl;
        return 0;
    }
    cout << "--------------------------------------------------" << endl;

    //=========================================================================//
    //                         生成子波（Ricker）                              //
    //=========================================================================//
    cout << "--------------------------------------------------" << endl;
    cout << "Calaulate the Wavelet" << endl;
    arma::Col<float> wavelet(nw);
    rickerwavelet_func(nw,dt, wavelet);
    cout << "--------------------------------------------------" << endl;

    //=========================================================================//
    //                         总循环（每个反射点进行一次）                    //
    //=========================================================================//
    //seis_xt_stack_all:X-T域NMO&Stack结果剖面
    //seis_thetat_stack_all:Theta-T域NMO&Stack结果剖面
    arma::Mat<float> seis_ori_stack_all(nt, nvx, fill::zeros);
    arma::Mat<float> seis_xt_stack_all(nt, nvx, fill::zeros);
    arma::Mat<float> seis_thetat_stack_all(nt, nvx, fill::zeros);
    int flag_save_nx = 0;
    cout << "     X     : " << nvx * dx << " m" << endl;
    cout << "     Z     : " << nvz * dz << " m" << endl;
    cout << "     Dt    : " << dt << " s" << endl;
    cout << "     Dx    : " << dx << " m" << endl;
    cout << "     Dz    : " << dz << " m" << endl;
    cout << "======================================================================" << endl;
    for (int ivx = 0; ivx < 1; ivx++)
    {

        //=============================道集合成===============================//
        cout << "--------------------------------------------------" << endl;
        cout << "Synthesize the Seismic Sections : " << endl;
        cout << "     Synthesize the CMP Section" << endl;

        //生成CMP(X-T)道集
        arma::Mat<float> seis_xt(nt, nx, fill::zeros);

        //generate_cmpgathers_func(seis_xt, v_rms_col,
                //wavelet, nw,
                //nx, nt, dx, dt , t0);

        arma::Col<float> model_col=model.col(ivx);
        generate_cmpgathers_modelbased_func(seis_xt, v_rms_col,
                model_col,
                wavelet, nw,
                dz, nvz,
                nx, nt, dx, dt);

        //X-T 转为 Theta-T道集
        int n_theta=nx;
        float theta_max = atan(tan(80*3.1415926/180));
        float d_theta = theta_max * 1.0 / ((n_theta - 1) / 2);
        arma::Mat<float> seis_thetat(nt,n_theta,fill::zeros);
        xt_to_thetat_func(seis_xt, seis_thetat,
                v_rms_col,
                nx,  nt, dx, dz, dt);

        //生成Theta-T道集（用来对比）
        arma::Mat<float> seis_thetat_modelbased(nt, nx, fill::zeros);

        generate_thetatgathers_modelbased_func(seis_thetat_modelbased,
        model_col,
        wavelet,  nw,
        dz, nvz, d_theta,
        n_theta, nt, dt);

        //generate_thetatgathers_func(seis_thetat_modelbased, v_rms_col,
                //wavelet, nw,
                //nx, nt,  dx, dt , t0, h);
        cout << "--------------------------------------------------" << endl;

        //=======================零偏移集合成=================================//
        cout << "--------------------------------------------------" << endl;
        cout << "     Synthesize the Zero_Offset Section" << endl;
        //零偏移距剖面zero_off_CMP
        arma::Mat<float> zero_off(nt, nx, fill::zeros);

        generate_cmpgathers_zerooffset_modelbased_func(zero_off,
                model_col,
                wavelet, nw,
                dz, nvz,
                nx, nt, dt);


        //零偏移距剖面zero_theta
        arma::Mat<float> zero_theta(nt, n_theta, fill::zeros);

        generate_thetatgathers_zerotheta_modelbased_func(zero_theta,
                model_col,
                wavelet,  nw,
                dz,  nvz,
                n_theta, nt, dt);

        cout << "--------------------------------------------------" << endl;

        //=========================================================================//
        //                              动校正（NMO）                              //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "NMO : " << endl;
        cout << "     NMO Correction of X-T Domain" << endl;

        //===================时间-空间域（CMP）动校正=========================//
        arma::Mat<float> seis_xt_after_nmo(nt, nx, fill::zeros);
        xt_nmo_func(seis_xt,seis_xt_after_nmo,
                v_rms_col,  nx,  nt,  dx,  dt);

        //op2
        arma::Mat<float> seis_xt_op2 = seis_xt_after_nmo - zero_off;


        cout << "       NMO Correction of T-Theta Domain" << endl;
        //====================时间-角度域动校正================================//
        arma::Mat<float> seis_thetat_after_nmo(nt, n_theta, fill::zeros);
        thetat_nmo_func(seis_thetat,seis_thetat_after_nmo,
                v_rms_col,
                nx, n_theta,  nt,
                dx, d_theta, dt);

        //OP2
        arma::Mat<float> seis_thetat_op2 = seis_thetat_after_nmo - zero_theta;
        cout << "--------------------------------------------------" << endl;



        cout << "--------------------------------------------------" << endl;
        //Theta-T域数据转换为X-T域
        // 参数说明：
        // seis_x_theta_xt_after_nmo: Theta-T域反变换到X-T域的剖面
        arma::Mat<float> seis_x_theta_xt_after_nmo(nt, nx, fill::zeros);
        arma::Mat<float> seis_xt_mappedfromthetat(nt, nx, fill::zeros);

        thetat_to_xt_func(seis_thetat_modelbased, seis_xt_mappedfromthetat,
                v_rms_col,
                nx, n_theta, nt, dx, d_theta, dz, dt);

        thetat_to_xt_func(seis_thetat_after_nmo, seis_x_theta_xt_after_nmo,
                v_rms_col,
                nx, n_theta, nt, dx, d_theta, dz, dt);


        //=========================================================================//
        //                              叠加成像                                   //
        //=========================================================================//
        //seis_thetat_stack :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        //seis_xt_stack     :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        arma::Mat<float> seis_thetat_stack(nt, 3, fill::zeros);

        seis_thetat_stack.col(0) = seis_thetat.col(n_theta / 2);
        for (int it = 1; it < nt; it++)
        {
            float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / (it * dt * v / 2));
            for (int ix = 0; ix < n_theta; ix++)
            {
                float theta = (ix - n_theta / 2) * d_theta;
                if (abs(theta) < theta_max_here)
                {
                    seis_thetat_stack(it, 1) = seis_thetat_stack(it, 1) + seis_thetat_after_nmo(it, ix);
                }
            }
            int n_theta_here = theta_max_here * 1.0 / d_theta * 2.0;
            seis_thetat_stack(it, 1) = seis_thetat_stack(it, 1) / n_theta_here;
        }
        seis_thetat_stack.col(2) = seis_thetat_stack.col(0) - seis_thetat_stack.col(1);

        arma::Mat<float> seis_xt_stack(nt, 3, fill::zeros);
        seis_xt_stack.col(0) = seis_xt.col(nx / 2);
        for (int ix = 0; ix < nx; ix++)
        {
            seis_xt_stack.col(1) = seis_xt_stack.col(1) + seis_xt_after_nmo.col(ix);
        }
        seis_xt_stack.col(1) = seis_xt_stack.col(1) / nx;
        seis_xt_stack.col(2) = seis_xt_stack.col(0) - seis_xt_stack.col(1);

        seis_ori_stack_all.col(ivx) = seis_xt_stack.col(0);
        seis_xt_stack_all.col(ivx) = seis_xt_stack.col(1);
        seis_thetat_stack_all.col(ivx) = seis_thetat_stack.col(1);

        cout << "--------------------------------------------------" << endl;

        //=========================================================================//
        //                              频谱分析                                   //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "Here We Begin to Do the Fourier-Transform" << endl;

        float fre_max = 1.0/dt;
        float df      = fre_max*1.0/nt;
        cout<<"     The Maxmial Frequence is : "<<fre_max<<" Hz"<<endl;
        cout<<"     The D-frequence is       : "<<df<<" Hz"<<endl;

        arma::Col<cx_float> seis_ori_fourier=fft(seis_xt.col(0));
        arma::Col<cx_float> seis_xt_stack_fourier=fft(seis_xt_stack.col(1));
        arma::Col<cx_float> seis_thetat_stack_fourier=fft(seis_thetat_stack.col(1));

        arma::Col<float> seis_ori_fourier_abs=abs(seis_ori_fourier);
        arma::Col<float> seis_xt_stack_fourier_abs=abs(seis_xt_stack_fourier);
        arma::Col<float> seis_thetat_stack_fourier_abs=abs(seis_thetat_stack_fourier);

        arma::Mat<float> seis_xt_thetat_fourier_abs(nt,3,fill::zeros);
        seis_xt_thetat_fourier_abs.col(0)=seis_ori_fourier_abs;
        seis_xt_thetat_fourier_abs.col(1)=seis_xt_stack_fourier_abs;
        seis_xt_thetat_fourier_abs.col(2)=seis_thetat_stack_fourier_abs;


        cout << "--------------------------------------------------" << endl;


        //=========================================================================//
        //                              文件保存                                   //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "Save Files" << endl;


        if (ivx == flag_save_nx)
        {
            cout << "     No." << flag_save_nx << " Trace" << endl;


            string fn_seis_xt = fn + "seis_x_t.dat";
            seis_xt.save(fn_seis_xt, raw_binary);

            string fn_seis_theta = fn + "seis_theta_t.dat";
            seis_thetat.save(fn_seis_theta, raw_binary);

            string fn_zero_off = fn + "seis_x_t_zeros_off.dat";
            zero_off.save(fn_zero_off, raw_binary);

            string fn_zero_theta = fn + "seis_theta_t_zeros_theta.dat";
            zero_theta.save(fn_zero_theta, raw_binary);

            string fn_seis_xt_after_nmo = fn + "seis_x_t_after_nmo.dat";
            seis_xt_after_nmo.save(fn_seis_xt_after_nmo, raw_binary);

            string fn_seis_xt_op2 = fn + "seis_x_t_op2.dat";
            seis_xt_op2.save(fn_seis_xt_op2, raw_binary);

            string fn_theta_after_nmo = fn + "seis_theta_t_after_nmo.dat";
            seis_thetat_after_nmo.save(fn_theta_after_nmo, raw_binary);

            string fn_thetat_op2 = fn + "seis_theta_t_op2.dat";
            seis_thetat_op2.save(fn_thetat_op2, raw_binary);

            string fn_seis_thetat_stack = fn + "seis_theta_t_stack.dat";
            seis_thetat_stack.save(fn_seis_thetat_stack, raw_binary);

            string fn_seis_xt_stack = fn + "seis_x_t_stack.dat";
            seis_xt_stack.save(fn_seis_xt_stack, raw_binary);


            string fn_seis_ori_fourier_abs=fn+"seis_ori_fourier_abs.dat";
            seis_ori_fourier_abs.save(fn_seis_ori_fourier_abs,raw_binary);

            string fn_seis_xt_stack_fourier_abs=fn+"seis_x_t_stack_fourier_abs.dat";
            seis_xt_stack_fourier_abs.save(fn_seis_xt_stack_fourier_abs,raw_binary);

            string fn_seis_thetat_stack_fourier_abs=fn+"seis_theta_t_stack_fourier_abs.dat";
            seis_thetat_stack_fourier_abs.save(fn_seis_thetat_stack_fourier_abs,raw_binary);

            string fn_seis_xt_thetat_fourier_abs=fn+"seis_ori_x_t_theta_t_fourier_abs.dat";
            seis_xt_thetat_fourier_abs.save(fn_seis_xt_thetat_fourier_abs,raw_binary);

            string fn_seis_x_theta_xt_after_nmo = fn + "seis_x_theta_xt_after_nmo.dat";
            seis_x_theta_xt_after_nmo.save(fn_seis_x_theta_xt_after_nmo, raw_binary);

            string fn_thetat_modelbased = fn + "seis_theta_t_modelbased.dat";
            seis_thetat_modelbased.save(fn_thetat_modelbased, raw_binary);

            string fn_seis_xt_mappedfromthetat = fn + "seis_x_t_mappedfromthetat.dat";
            seis_xt_mappedfromthetat.save(fn_seis_xt_mappedfromthetat, raw_binary);

        }

        cout << "--------------------------------------------------" << endl;
    }

    string fn_seis_xt_stack_all = fn + "/seis_x_t_stack_all.dat";
    seis_xt_stack_all.save(fn_seis_xt_stack_all, raw_binary);

    string fn_seis_ori_stack_all = fn + "seis_ori_stack_all.dat";
    seis_ori_stack_all.save(fn_seis_ori_stack_all, raw_binary);

    string fn_seis_thetat_stack_all = fn + "seis_theta_t_stack_all.dat";
    seis_thetat_stack_all.save(fn_seis_thetat_stack_all, raw_binary);

    cout << "======================================================================" << endl;

    cout << "--------------------------------------------------" << endl;
    cout << "FINISH" << endl;
    finish = clock();
    cout << "Time = " << double(finish - start) / CLOCKS_PER_SEC << "s" << endl;
    cout << "--------------------------------------------------" << endl;

    return 0;
}
