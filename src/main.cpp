//----------------------------------------------------------------------------//
//程序目的： 1、分析 时间-偏移距 空间NMO带来的子波拉伸问题
//           2、分析 时间-角度 空间NMO子波拉伸问题
//程序原理：
//          拿到CMP道集后，1、进行NMO 2、转变为Theta-T域（基于平层假设）进行NMO，
//          并进行叠加成像。
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

int main(int argc, char **argv)
{
    if (argc<=1)
    {
        cout<<"======================================================================"<<endl;
        cout<<" Wavelet Strecting Analysis In NMO & PSTM "<<endl;
        cout<<" 1. Analysis the shape change owing to NMO in XT domain(CMP) "<<endl;
        cout<<" 2. Analysis the shape change owing to NMO in Theta domain "<<endl<<endl;
        cout<<" Required environment : lib: armadillo"<<endl;
        cout<<"                        compiler:c++11 or c++14 c++17"<<endl;
        cout<<" Copyright：2020-   "<<endl;
        cout<<"            WPI TONGJI University"<<endl;
        cout<<" Author   ：ShengShen"<<endl;
        cout<<" Time      :2020 07 06"<<endl;
        cout<<"======================================================================"<<endl<<endl;;
        cout<<"ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR!"<<endl<<endl;;
        cout<<"            Please input the parameter_card path !"<<endl<<endl;;
        cout<<"ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR! ERROR!"<<endl<<endl;;
        cout<<"======================================================================"<<endl;
        return 0;
    }

    const string fn_par(argv[1]);



    cout<<"======================================================================"<<endl;
    cout<<" Wavelet Strecting Analysis In NMO & PSTM "<<endl;
    cout<<" 1. Analysis the shape change owing to NMO in XT domain(CMP) "<<endl;
    cout<<" 2. Analysis the shape change owing to NMO in Theta domain "<<endl<<endl;
    cout<<" Required environment : lib: armadillo"<<endl;
    cout<<"                        compiler:c++11 or c++14 c++17"<<endl;
    cout<<" Copyright：2020-   "<<endl;
    cout<<"            WPI TONGJI University"<<endl;
    cout<<" Author   ：ShengShen"<<endl;
    cout<<" Time      :2020 07 06"<<endl;
    cout<<"======================================================================"<<endl;
    cout << "START" << endl;
    clock_t start, finish;
    start = clock();
    //=========================================================================//
    //                         参数卡读取 & 定义速度模                         //
    //=========================================================================//
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

    //参数卡读取
    Parameter par;
    readpar_func(fn_par, par);
    printpar_func(par);
    string fn = par.fn_path_of_output;
    int nvx  = par.nvx;
    int nvz  = par.nvz;
    int nt   = par.nt;
    float dt = par.dt;
    int nw   = par.nw;
    int nx   = par.nx;
    float dx = par.dx;
    float dz = par.dz;

    //速度模型加载
    //      速度模型model
    //      均方根速度v_rms
    //      均匀速度v_mean
    arma::Mat<float> model(nvz, nvx, fill::zeros);
    arma::Mat<float> v_rms(nt, nvx, fill::zeros);
    model.load(par.fn_path_of_model_in, raw_binary);
    v_rms.load(par.fn_path_of_v_rms_in, raw_binary);
    v_rms.reshape(nt, nvx);
    model.reshape(nvz, nvx);

    //由速度模型搜索反射点位置和个数
    //ref_num:      每一道有多少反射界面
    //ref_location: 每一道的每一个反射界面的深度（米）
    //ref_location_t0: 每一道的每一个反射界面的到达时（s）
    arma::Col<float> ref_num(nvx, fill::zeros);
    arma::Mat<float> ref_location(nvz, nvx, fill::zeros);
    arma::Mat<float> ref_location_t0(nvz, nvx, fill::zeros);

    //=========================================================================//
    //                         生成子波（Ricker）                              //
    //=========================================================================//
    cout<<"======================================================================"<<endl;
    cout << "Calaulate the Wavelet" << endl;
    arma::Col<float> wavelet(nw);
    rickerwavelet_func(nw,dt, wavelet);
    cout<<"======================================================================"<<endl;

    //=========================================================================//
    //                         总循环（每个反射点进行一次）                    //
    //=========================================================================//
    cout<<"======================================================================"<<endl;
    cout<<"Generate Gathers & NMO & Stack"<<endl;
    //seis_xt_stack_all:X-T域NMO&Stack结果剖面
    //seis_thetat_stack_all:Theta-T域NMO&Stack结果剖面
    // flag_save_nx :保存剖面的地震位置ivx

    int flag_save_nx = 0;
    int nvx_end = 10;
    cout << "Only No." << flag_save_nx << " trace will be saved" << endl;

    arma::Mat<float> seis_ori_stack_all(nt, nvx, fill::zeros);
    arma::Mat<float> seis_xt_stack_all(nt, nvx, fill::zeros);
    arma::Mat<float> seis_thetat_stack_all(nt, nvx, fill::zeros);

    for (int ivx = 0; ivx < nvx_end; ivx++)
    {
        //**进度条显示**//
        //cout.width(3);
        cout<<"["<<(ivx+1)*100/nvx_end<<"]"<< "%";
        fflush(stdout);
        std::cout << "\b\b\b\b\b";


        arma::Col<float> v_rms_col=v_rms.col(0);
        arma::Col<float> model_col=model.col(ivx);
        //====================道集合成========================//
        DEBUG_PRINTF_("Synthesize the Seismic Sections : \n");
        DEBUG_PRINTF_("     Synthesize the CMP Section\n");

        //生成CMP(X-T)道集
        arma::Mat<float> seis_xt(nt, nx, fill::zeros);
        generate_cmpgathers_modelbased_func(seis_xt, v_rms_col,
                model_col,
                wavelet, nw,
                dz, nvz,
                nx, nt, dx, dt);

        //X-T 转为 Theta-T道集
        int n_theta=nx;
        float theta_max = atan(tan(80*3.1415926/180));
        float d_theta = theta_max * 1.0 / (n_theta/ 2);
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

        //======================零偏移合成===================//
        DEBUG_PRINTF_("     Synthesize the Zero_Offset Section\n");
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

        //====================================================================//
        //                              动校正（NMO）                         //
        //====================================================================//
        DEBUG_PRINTF_("NMO:\n");
        DEBUG_PRINTF_("     NMO Correction of X-T Domain\n");
        //=====================X-T域（CMP）动校正=========================//
        //      xt_nmo_type = 0:常规XT域NMO校正（校正到零偏移距处）
        //      xt_nmo_type = 1:非常规XT域NMO校正（校正到指定偏移距处）
        //      xt_nmo_type = 2:非常规XT域NMO校正（考虑ThetaT域校正）
        //      xt_nmo_type = 3:非常规XT域NMO校正（考虑ThetaT域校正）
        int xt_nmo_type= 1;

        arma::Mat<float> seis_xt_after_nmo(nt, nx, fill::zeros);

        if(xt_nmo_type == 0)
        {
            xt_nmo_func(seis_xt,seis_xt_after_nmo,
                    v_rms_col,  nx,  nt,  dx,  dt);
        }
        if(xt_nmo_type == 1)
        {
            int pilottrace = nx;
            //NMO到指定参考道集位置(pilottrace)处 seis_xt_after_nmo
            xt_nmo_pilot_trace_based_func(seis_xt, seis_xt_after_nmo,
                    v_rms_col, nx, nt, dx, dt, pilottrace);

            //生成指定参考道集处(pilottrace)的零相位子波剖面zero_off
            generate_cmpgathers_zerooffset_pilot_trace_based_func(
                    zero_off,
                    model_col, v_rms_col,
                    wavelet,  nw,
                    dz, nvz,  pilottrace,
                    nx, nt,   dx,  dt);
        }
        if(xt_nmo_type == 2)
        {
            thetat_nmo_xtbased_func(seis_xt,seis_xt_after_nmo,
                    v_rms_col, nx, n_theta, nt, dx, dt, d_theta);
        }
        if(xt_nmo_type == 3)
        {
            xt_nmo_thetatbased_func( seis_xt,seis_xt_after_nmo,
                    v_rms_col, nx, nt, dx, dt);
        }

        //op2
        arma::Mat<float> seis_xt_op2 = seis_xt_after_nmo - zero_off;


        DEBUG_PRINTF_("       NMO Correction of T-Theta Domain\n");
        //====================Theta-T域动校正=============================//
        arma::Mat<float> seis_thetat_after_nmo(nt, n_theta, fill::zeros);
        thetat_nmo_func(seis_thetat,seis_thetat_after_nmo,
                v_rms_col,
                nx, n_theta,  nt,
                dx, d_theta, dt);

        //OP2
        arma::Mat<float> seis_thetat_op2 = seis_thetat_after_nmo - zero_theta;


        //Theta-T域数据转换为X-T域
        // 参数说明：
        // seis_x_theta_xt_after_nmo: Theta-T域反变换到X-T域的剖面
        arma::Mat<float> seis_x_theta_xt_after_nmo(nt, nx, fill::zeros);
        arma::Mat<float> seis_xt_mappedfromthetat(nt, nx, fill::zeros);

        thetat_to_xt_func(seis_thetat, seis_xt_mappedfromthetat,
                v_rms_col,
                nx, n_theta, nt, dx, d_theta, dz, dt);

        thetat_to_xt_func(seis_thetat_after_nmo, seis_x_theta_xt_after_nmo,
                v_rms_col,
                nx, n_theta, nt, dx, d_theta, dz, dt);


        //====================================================================//
        //                              叠加成像                              //
        //====================================================================//
        //seis_thetat_stack :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        //seis_xt_stack     :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        arma::Mat<float> seis_thetat_stack(nt, 3, fill::zeros);

        seis_thetat_stack.col(0) = zero_theta.col(n_theta / 2);
        for (int it = 1; it < nt; it++)
        {
            float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 /
                    (it * dt * v_rms_col(it) / 2));
            for (int ix = 0; ix < n_theta; ix++)
            {
                float theta = (ix - n_theta / 2) * d_theta;
                if (abs(theta) < theta_max_here)
                {
                    seis_thetat_stack(it, 1) = seis_thetat_stack(it, 1)
                        + seis_thetat_after_nmo(it, ix);
                }
            }
            int n_theta_here = theta_max_here * 1.0 / d_theta * 2.0;
            seis_thetat_stack(it, 1) = seis_thetat_stack(it, 1) / n_theta_here;
        }
        seis_thetat_stack.col(2) = seis_thetat_stack.col(0)- seis_thetat_stack.col(1);

        arma::Mat<float> seis_xt_stack(nt, 3, fill::zeros);
        seis_xt_stack.col(0) = zero_off.col(nx / 2);
        for (int ix = 0; ix < nx; ix++)
        {
            seis_xt_stack.col(1) = seis_xt_stack.col(1)
                + seis_xt_after_nmo.col(ix);
        }
        seis_xt_stack.col(1) = seis_xt_stack.col(1) / nx;
        seis_xt_stack.col(2) = seis_xt_stack.col(0) - seis_xt_stack.col(1);

        seis_ori_stack_all.col(ivx)    = seis_xt_stack.col(0);
        seis_xt_stack_all.col(ivx)     = seis_xt_stack.col(1);
        seis_thetat_stack_all.col(ivx) = seis_thetat_stack.col(1);


        //====================================================================//
        //                              频谱分析                              //
        //====================================================================//
        //cout << "Here We Begin to Do the Fourier-Transform" << endl;

        float fre_max = 1.0/dt;
        float df      = fre_max*1.0/nt;
        DEBUG_PRINTF_("     The Maxmial Frequence is : %f Hz \n", &fre_max);
        DEBUG_PRINTF_("     The D-frequence is       : %f Hz \n", &df);

        arma::Col<cx_float> seis_ori_fourier=fft(zero_off.col(nx/2));
        arma::Col<cx_float> seis_xt_stack_fourier=fft(seis_xt_stack.col(1));
        arma::Col<cx_float> seis_thetat_stack_fourier=fft(seis_thetat_stack.col(1));

        arma::Col<float> seis_ori_fourier_abs=abs(seis_ori_fourier);
        arma::Col<float> seis_xt_stack_fourier_abs=abs(seis_xt_stack_fourier);
        arma::Col<float> seis_thetat_stack_fourier_abs=abs(seis_thetat_stack_fourier);

        arma::Mat<float> seis_xt_thetat_fourier_abs(nt,3,fill::zeros);
        seis_xt_thetat_fourier_abs.col(0)=seis_ori_fourier_abs;
        seis_xt_thetat_fourier_abs.col(1)=seis_xt_stack_fourier_abs;
        seis_xt_thetat_fourier_abs.col(2)=seis_thetat_stack_fourier_abs;

        //====================================================================//
        //                              文件保存                              //
        //====================================================================//
        DEBUG_PRINTF_("Save Files\n");

        if (ivx == flag_save_nx)
        {
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

    }

    string fn_seis_xt_stack_all = fn + "/seis_x_t_stack_all.dat";
    seis_xt_stack_all.save(fn_seis_xt_stack_all, raw_binary);

    string fn_seis_ori_stack_all = fn + "seis_ori_stack_all.dat";
    seis_ori_stack_all.save(fn_seis_ori_stack_all, raw_binary);

    string fn_seis_thetat_stack_all = fn + "seis_theta_t_stack_all.dat";
    seis_thetat_stack_all.save(fn_seis_thetat_stack_all, raw_binary);

    printf("\n");
    cout<<"======================================================================"<<endl;
    cout << "FINISH" << endl;
    finish = clock();
    cout << "Time = " << double(finish - start) / CLOCKS_PER_SEC << "s" << endl;
    cout<<"======================================================================"<<endl;

    return 0;
}
