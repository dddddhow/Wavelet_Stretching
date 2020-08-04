//----------------------------------------------------------------------------//
//程序目的： 1、分析 时间-偏移距 空间NMO带来的子波拉伸问题
//           2、分析 时间-角度 空间NMO子波拉伸问题
//程序原理：
//          从模型出发，给定一个模型，解析计算其CMP道集和Theta-T道集，
//          然后分别进行NMO校正，并对比结果
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
//程序依赖关系说明：
//          1、需要armadillo库
//          2、-std=c++11 或者 更高
//
//Copyright：2020-
//          WPI TONGJI University
//Author  ：ShengShen
//Time    ：2020 07 22
//----------------------------------------------------------------------------//
#include <armadillo>
#include <ctime>
#include <iostream>

using namespace std;
using namespace arma;

//-------------------------------子函数---------------------------------------//
double rickerwavelet_func(int nw, arma::Col<float> &signal);
//程序说明：
//      生成地震子波（Ricker Wavelet）
//参数说明：
//      nw：子波长度
//      signal：子波（向量）
//----------------------------------------------------------------------------//
int main(int argc, char **argv)
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

    int nvx = 200;
    int nvz = 400;

    int nx = 201;
    int nt = 2000;

    int nw = 81;

    float dx = 5.0;
    float dz = 5.0;
    float dt = 0.001;

    //文件保存路径 fn
    string fn = "../file_xtTothetatToxt/";

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
    ref_location_t0 = ref_location_t0 * 2.0;

    /*
    for (int ir = 0; ir < ref_num(0); ir++)
    {
        float t = 0;
        for (int iz = 0; iz < ref_location(ir, 0) / dz; iz++)
        {
            t = t + dz * 1.0 / model(iz, 0);
        }
        cout << " No." << ir << "  is : " << t*2 - ref_location_t0(ir, 0) << endl;
        cout << " t is :" << t << endl;
        cout << " No." << ir << " by  v_mean  is : " << ref_location(ir, 0) - v_rms(int(t*1.0 / dt), 0) * t  << endl;
    }

    return 0;
*/
    cout << "--------------------------------------------------" << endl;

    //=========================================================================//
    //                         生成子波（Ricker）                              //
    //=========================================================================//
    cout << "--------------------------------------------------" << endl;
    cout << "Here We Begin to Calaulate the Wavelet" << endl;
    //生成子波 ricker wavelet
    arma::Col<float> wavelet(nw);
    rickerwavelet_func(nw, wavelet);
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
    //cout<<"     V     : "<< v <<" m/s"<<endl;
    cout << "======================================================================" << endl;
    for (int ivx = 0; ivx < 1; ivx++)
    {
        cout << "--------------------------------------------------" << endl;
        //=========================================================================//
        //                              正演模拟（基于褶积）                       //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "Here We Begin to Synthesize the Seismic Sections" << endl;
        cout << "     Here We Begin to Synthesize the CMP Section" << endl;
        //时间-空间域剖面（CMP）合成（基于褶积运算）
        //参数说明：
        //      ref_cof_xt :反射系数矩阵
        //      seis_xt    :合成地震记录（CMP）

        arma::Mat<float> ref_cof_xt(nt, nx, fill::zeros);
        arma::Mat<float> seis_xt(nt, nx, fill::zeros);

        for (int ir = 0; ir < ref_num(ivx); ir++)
        {
            for (int ix = 0; ix < nx; ix++)
            {
                float x_temp = (ix - nx / 2) * dx;
                float t0_temp = ref_location_t0(ir, ivx);
                int iz = ref_location(ir, ivx) / dz;
                float v = v_rms(iz, ivx);
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

        cout << "--------------------------------------------------" << endl;

        // CMP Domain-> Theta-T Domain
        //时间-角度域剖面（T-Theta）合成（基于褶积运算）
        //参数说明：
        //      ref_cof_thetat :反射稀疏矩阵
        //      seis_thetat    :合成地震记录（CMP）
        int n_theta = nx;
        float theta_max = atan(tan(80 * 3.1415926 * 1.0 / 180));
        float d_theta = theta_max * 1.0 / (n_theta / 2);
        cout << "         Theta max is      : " << theta_max << " radians";
        cout << " : " << theta_max * 180 * 1.0 / 3.1415926 << " degrees" << endl;
        cout << "         D_Theta  is       : " << d_theta << " radinans";
        cout << " : " << d_theta * 180 * 1.0 / 3.1415926 << " degrees" << endl;

        arma::Mat<float> ref_cof_thetat(nt, n_theta, fill::zeros);
        arma::Mat<float> seis_thetat(nt, n_theta, fill::zeros);

        //CMP道集加密
        cout << "     Here We Encrypted CMP gathers" << endl;
        arma::Mat<float> seis_xt_interp;
        {
            float x_interp_times = 50;
            float t_interp_times = 50;
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

        cout << "         The NX of Interpolated CMP is :" << nx_interp << endl;
        cout << "         The dx of Interpolated CMP is :" << dx_interp << " m" << endl;
        cout << "         The NT of Interpolated CMP is :" << nt_interp << endl;
        cout << "         The dt of Interpolated CMP is :" << dt_interp << " s" << endl;

        cout << "     Here We Transform CMP Gathers to Theta-T Gathers" << endl;
        for (int it = 0; it < nt; it++)
        {
            //给出t 计算h
            float h_temp = 0;
            float t_htemp = 0;
            float v_htemp = 0;
            for (int ir = 0; ir < ref_num(ivx); ir++)
            {
                if (it * dt >= ref_location_t0(ir, ivx))
                {
                    h_temp = ref_location(ir, ivx);
                    t_htemp = it * dt - ref_location_t0(ir, ivx);
                    v_htemp = model(ref_location(ir, ivx) / dz + 1, ivx);
                }
            }
            h_temp = h_temp + t_htemp * v_htemp;

            for (int i_theta = 0; i_theta < n_theta; i_theta++)
            {
                float theta = (i_theta - n_theta / 2) * d_theta;
                float theta_max_here = atan((nx / 2) * dx * 1.0 / h_temp);
                if (abs(theta) < theta_max_here)
                {
                    //给定 theta t0，确定x
                    float t0_temp = it * dt;
                    float x = v_mean(it, ivx) * 1.0 / 2 * t0_temp * tan(theta);
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
        //=========================================================================//
        //                              动校正（NMO）                              //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "Here We Begin to calculate the NMO Correction" << endl;

        cout << "     Here We Synthesize the Zero_Theta Section" << endl;
        //零偏移距剖面zero_theta
        /*
        arma::Mat<float> zero_theta(nt, n_theta, fill::zeros);
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iz = 0; iz < nw; iz++)
            {
                float theta = (ix - n_theta / 2) * d_theta;
                float theta_max_here = atan((n_theta / 2) * dx * 1.0 / h);
                if (abs(theta) < theta_max_here)
                {
                    zero_theta(iz + nh - nw / 2 + 1, ix) = wavelet(iz);
                }
            }
        }
    */
        cout << "     Here We Synthesize the Zero_Offset Section" << endl;
        //零偏移距剖面zero_off_CMP
        arma::Mat<float> zero_off(nt, nx, fill::zeros);
        for (int ir = 0; ir < ref_num(ivx); ir++)
        {
            int nh = ref_location_t0(ir, ivx) * 1.0 / dt;
            for (int ix = 0; ix < nx; ix++)
            {
                for (int iz = 0; iz < nw; iz++)
                {
                    zero_off(iz + nh - nw / 2 + 1, ix) = wavelet(iz);
                }
            }
        }

        cout << "     Here We Begin to Calculate the NMO Correction of CMP Section" << endl;
        //时间-空间域（CMP）动校正
        //参数说明：
        //      t_xt_norm ：按照时距关系公式计算的每一点的时间
        //      t_xt_nmo  ：动校正量
        //      seis_xt_after_nmo：动校正后的剖面
        //      seis_xt_op2：动校正后域零偏移距剖面做差结果
        arma::Mat<float> t_xt_norm(nt, nx, fill::zeros);
        arma::Mat<float> t_xt_nmo(nt, nx, fill::zeros);
        arma::Mat<float> seis_xt_after_nmo(nt, nx, fill::zeros);

        for (int ix = 0; ix < nx; ix++)
        {
            for (int it = 0; it < nt; it++)
            {
                float t0_temp = it * dt;
                float x_temp = (ix - nx / 2) * dx;
                float t = sqrt(t0_temp * t0_temp + 4.0 * x_temp * x_temp * 1.0 / v_rms(it, ivx) * 1.0 / v_rms(it, ivx));
                t_xt_norm(it, ix) = t;
                t_xt_nmo(it, ix) = t - t0_temp;
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

                int it_ori = int(abs(max_sup));
                //cout<<it_ori<<endl;
                if (it_ori + 1 < nt)
                {
                    seis_xt_after_nmo(it, ix) = b * seis_xt(it_ori, ix) + a * seis_xt(it_ori + 1, ix);
                }
            }
        }
        //OP2
        arma::Mat<float> seis_xt_op2 = seis_xt_after_nmo - zero_off;

        cout << "--------------------------------------------------" << endl;
        cout << "     Here We Begin to Calculate the NMO Correction of T-Theta Section" << endl;
        //时间-角度域动校正
        //参数说明：
        //      t_thetat_norm ：按照时距关系公式计算的每一点的时间
        //      t_thetat_nmo  ：动校正量
        //      seis_thetat_after_nmo：动校正后的剖面
        //      seis_thetat_op2：动校正后域零偏移距剖面做差结果
        arma::Mat<float> t_thetat_norm(nt, n_theta, fill::zeros);
        arma::Mat<float> t_thetat_nmo(nt, n_theta, fill::zeros);
        arma::Mat<float> seis_thetat_after_nmo(nt, n_theta, fill::zeros);

        for (int ix = 0; ix < n_theta; ix++)
        {
            for (int it = 1; it < nt; it++)
            {
                //给出t 计算h
                float h_temp = 0;
                float t_htemp = 0;
                float v_htemp = 0;
                for (int ir = 0; ir < ref_num(ivx); ir++)
                {
                    if (it * dt >= ref_location_t0(ir, ivx))
                    {
                        h_temp = ref_location(ir, ivx);
                        t_htemp = it * dt - ref_location_t0(ir, ivx);
                        v_htemp = model(ref_location(ir, ivx) / dz + 1, ivx);
                    }
                }
                h_temp = h_temp + t_htemp * v_htemp;

                float theta = (ix - (nx - 1) / 2) * d_theta;
                float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / h_temp);
                if (abs(theta) < theta_max_here)
                {
                    // float h_temp = v_mean(it, ivx) * it * dt / 2;
                    float t0_temp = it * dt;
                    float t = sqrt(t0_temp * t0_temp + 4.0 * h_temp * h_temp * sin(theta) * sin(theta) * 1.0 / v_rms(it, ivx) * 1.0 / v_rms(it, ivx) * 1.0 / cos(theta) * 1.0 / cos(theta));

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
                float theta = (ix - (nx - 1) / 2) * d_theta;
                float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / (it * dt * v_mean(it, ivx) / 2));
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

        //OP2
        //arma::Mat<float> seis_thetat_op2 = seis_thetat_after_nmo - zero_theta;

        cout << "--------------------------------------------------" << endl;

        //=========================================================================//
        //                              叠加成像                                   //
        //=========================================================================//
        //seis_thetat_stack :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        //seis_xt_stack     :叠加数据   第一道表示零偏移距数据（即真实地震数据）
        //                              第二道表示叠加后数据
        //                              第三道为前两道相减
        cout << "--------------------------------------------------" << endl;
        cout << "Here We Begin to Stack the Gathers" << endl;
        arma::Mat<float> seis_thetat_stack(nt, 3, fill::zeros);

        seis_thetat_stack.col(0) = seis_thetat.col(n_theta / 2);
        for (int it = 1; it < nt; it++)
        {
            float theta_max_here = atan(((nx - 1) / 2) * dx * 1.0 / (it * dt * v_mean(it, ivx) / 2));
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
        cout << "--------------------------------------------------" << endl;

        //=========================================================================//
        //                              频谱分析                                   //
        //=========================================================================//
        cout << "--------------------------------------------------" << endl;
        cout << "Here We Begin to Do the Fourier-Transform" << endl;
        float fre_max = 1.0 / dt;
        float df = fre_max * 1.0 / nt;
        cout << "     The Maxmial Frequence is : " << fre_max << " Hz" << endl;
        cout << "     The D-frequence is       : " << df << " Hz" << endl;

        arma::Col<cx_float> seis_ori_fourier = fft(seis_xt.col(0));
        arma::Col<cx_float> seis_xt_stack_fourier = fft(seis_xt_stack.col(1));

        arma::Col<float> seis_ori_fourier_abs = abs(seis_ori_fourier);
        arma::Col<float> seis_xt_stack_fourier_abs = abs(seis_xt_stack_fourier);

        arma::Mat<float> seis_xt_thetat_fourier_abs(nt, 2, fill::zeros);
        seis_xt_thetat_fourier_abs.col(0) = seis_ori_fourier_abs;
        seis_xt_thetat_fourier_abs.col(1) = seis_xt_stack_fourier_abs;

        cout << "--------------------------------------------------" << endl;

        seis_ori_stack_all.col(ivx) = seis_xt_stack.col(0);
        seis_xt_stack_all.col(ivx) = seis_xt_stack.col(1);

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
            //zero_theta.save(fn_zero_theta, raw_binary);

            string fn_seis_xt_nmo = fn + "seis_x_t_nmo.dat";
            t_xt_nmo.save(fn_seis_xt_nmo, raw_binary);

            string fn_seis_xt_after_nmo = fn + "seis_x_t_after_nmo.dat";
            seis_xt_after_nmo.save(fn_seis_xt_after_nmo, raw_binary);

            string fn_seis_xt_op2 = fn + "seis_x_t_op2.dat";
            seis_xt_op2.save(fn_seis_xt_op2, raw_binary);

            string fn_thetat_nmo = fn + "seis_theta_t_nmo.dat";
            t_thetat_nmo.save(fn_thetat_nmo, raw_binary);

            string fn_theta_after_nmo = fn + "seis_theta_t_after_nmo.dat";
            seis_thetat_after_nmo.save(fn_theta_after_nmo, raw_binary);

            string fn_thetat_op2 = fn + "seis_theta_t_op2.dat";
            //seis_thetat_op2.save(fn_thetat_op2, raw_binary);

            string fn_seis_thetat_stack = fn + "seis_theta_t_stack.dat";
            seis_thetat_stack.save(fn_seis_thetat_stack, raw_binary);

            string fn_seis_xt_stack = fn + "seis_x_t_stack.dat";
            seis_xt_stack.save(fn_seis_xt_stack, raw_binary);
        }
        cout << "--------------------------------------------------" << endl;
    }

    string fn_seis_xt_stack_all = fn + "/seis_x_t_stack_all.dat";
    seis_xt_stack_all.save(fn_seis_xt_stack_all, raw_binary);

    string fn_seis_ori_stack_all = fn + "seis_ori_stack_all.dat";
    seis_ori_stack_all.save(fn_seis_ori_stack_all, raw_binary);

    cout << "======================================================================" << endl;
    cout << "--------------------------------------------------" << endl;
    cout << "FINISH" << endl;
    finish = clock();
    cout << "Time = " << double(finish - start) / CLOCKS_PER_SEC << "s" << endl;
    cout << "--------------------------------------------------" << endl;

    return 0;
}

double rickerwavelet_func(int nw, arma::Col<float> &signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
    float t, t1;
    float temp[nw];
    float freq = 28;
    float PI = 3.1415926;
    float dt = 0.001;

    for (i = 0; i < nw; i++)
    {
        t = dt * i;
        t1 = 1.0 / freq; //双边雷克子波
        //t1=0;//单边雷克子波
        temp[i] = 10 * (1 - 2 * PI * PI * freq * freq * (t - t1) * (t - t1)) * exp(-PI * PI * freq * freq * (t - t1) * (t - t1));
    }

    //导数形式
    /*    signal[0]=0;*/
    //for (i=1;i<param->Lw;i++)
    //{
    //signal[i]=(temp[i]-temp[i-1])*1.0/param->dt;
    /*    }*/
    //雷克子波形式
    for (i = 0; i < nw; i++)
    {
        signal[i] = temp[i];
    }
    return 0.0;
}
