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

#include <iostream>
#include <armadillo>
#include <cmath>
#include <ctime>
#include "ShengShen.h"
using namespace std;
using namespace arma;

//-------------------------------子函数---------------------------------------//
double rickerwavelet_func(int nw, arma::Col<float> & signal);
//程序说明：
//      生成地震子波（Ricker Wavelet）
//参数说明：
//      nw：子波长度
//      signal：子波（向量）

int ReadPar_func(const string & fn_par, Parameter & par);

int PrintPar_func(Parameter & Par);

//----------------------------------------------------------------------------//

int main(int argc, char **argv )
{
    cout<<"--------------------------------------------------"<<endl;
    cout<<"START"<<endl;
    clock_t start,finish;
    start = clock();
    cout<<"--------------------------------------------------"<<endl;

    //参数卡读取
    cout<<"--------------------------------------------------"<<endl;
    string fn_in="./NMO_Ana.par";
    Parameter par;
    cout<<"Read Parameters From File"<<endl;
    ReadPar_func( fn_in , par);
    cout<<"Read Parameters From File"<<endl;
    PrintPar_func(par);

    int nt   = 2000;
    float dt = 0.0005;

    int nw   = 81;
    float v  = 2000;

    int h    = 200;
    float t0 = 2.0*h*1.0/v;
    int nh   = t0*1.0/dt;

    int nx   = 101;
    float dx = 5.0;
    float dz = 5.0;
    int nz   = int(nt*dt*v/dz/2);

    cout<<"     X     : "<< nx*dx <<" m"<<endl;
    cout<<"     Z     : "<< nz*dz <<" m"<<endl;
    cout<<"     Dt    : "<< dt <<" s"<<endl;
    cout<<"     Dx    : "<< dx <<" m"<<endl;
    cout<<"     Dz    : "<< dz <<" m"<<endl;
    cout<<"     V     : "<< v <<" m/s"<<endl;
    cout<<"     Depth : "<< h <<" m"<<endl;

    if(nh<nw/2)
    {
        cout<<"Error ! The Depth Minimum is "<<(nw/2)*dz<<" m"<<endl;
        return 0;
    }
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Here We Begin to Calaulate the Wavelet"<<endl;
    //生成子波 ricker wavelet
    arma::Col<float> wavelet(nw);
    rickerwavelet_func(nw,wavelet);


    cout<<"--------------------------------------------------"<<endl;
    cout<<"Here We Begin to Synthesize the Seismic Sections"<<endl;
    cout<<"     Here We Begin to Synthesize the CMP Section"<<endl;
    //时间-空间域剖面（CMP）合成（基于褶积运算）
    //参数说明：
    //      ref_cof_xt :反射系数矩阵
    //      seis_xt    :合成地震记录（CMP）
    arma::Mat<float> ref_cof_xt(nt,nx,fill::zeros);
    arma::Mat<float> seis_xt(nt,nx,fill::zeros);

    for(int ix=0; ix<nx; ix++)
    {
        float x_temp=(ix-nx/2)*dx;
        float t0_temp=t0;
        float t=sqrt(t0_temp*t0_temp+4.0*x_temp*x_temp*1.0/v*1.0/v);
        int it=int(t/dt)+nw/2;
        if(it<nt)
        {
            ref_cof_xt(it,ix)=1;
        }
    }

    for(int ix=0; ix<nx; ix++)
    {
        arma::Col<float> vec_temp=conv(ref_cof_xt.col(ix),wavelet);
        seis_xt.col(ix)=vec_temp(span(nw-1,nt+nw-2));
        //seis_xt.col(ix)=vec_temp(span(0,nt-1));
    }

    cout<<"     Here We Begin to Sythnesize the T-Theta Section"<<endl;


    // CMP Domain-> Theta-T Domain
    //时间-角度域剖面（T-Theta）合成（基于褶积运算）
    //参数说明：
    //      ref_cof_thetat :反射稀疏矩阵
    //      seis_thetat    :合成地震记录（CMP）
    int n_theta=nx;
    float theta_max=atan((nx-1)*dx*1.0/2*1.0/(h/2));
    float theta_max_here=atan((nx-1)*dx*1.0/2*1.0/h);
    float d_theta=theta_max*1.0/((n_theta-1)/2);
    cout<<"         Theta max is      : "<<theta_max<<" radians";
    cout<<" : "<<theta_max*180*1.0/3.1415926<<" degrees"<<endl;
    cout<<"         Theta max here is : "<<theta_max_here<<" radians";
    cout<<" : "<<theta_max_here*180*1.0/3.1415926<<" degrees"<<endl;
    cout<<"         D_Theta  is       : "<<d_theta<<" radinans";
    cout<<" : "<<d_theta*180*1.0/3.1415926<<" degrees"<<endl;

    arma::Mat<float> ref_cof_thetat(nt,n_theta,fill::zeros);
    arma::Mat<float> seis_thetat(nt,n_theta,fill::zeros);

    //CMP道集加密
    cout<<"     Here We Encrypted CMP gathers"<<endl;
    arma::Mat<float> seis_xt_interp;
    {
        float x_interp_times=50;
        float t_interp_times=50;
        arma::fvec x=regspace<fvec>(0,nx-1);
        arma::fvec y=regspace<fvec>(0,nt-1);

        arma::fvec xi=regspace<fvec>(0, 1.0/x_interp_times, nx-1);
        arma::fvec yi=regspace<fvec>(0, 1.0/t_interp_times, nt-1);

        interp2(x, y, seis_xt, xi, yi, seis_xt_interp);
    }
    //seis_xt_interp.save("../file_DataBased/seis_x_t_interp.dat",raw_binary);

    int nx_interp=seis_xt_interp.n_cols;
    int nt_interp=seis_xt_interp.n_rows;
    float dx_interp=nx*1.0/nx_interp*dx;
    float dt_interp=nt*1.0/nt_interp*dt;

    cout<<"         The NX of Interpolated CMP is :"<<nx_interp<<endl;
    cout<<"         The dx of Interpolated CMP is :"<<dx_interp<<" m"<<endl;
    cout<<"         The NT of Interpolated CMP is :"<<nt_interp<<endl;
    cout<<"         The dt of Interpolated CMP is :"<<dt_interp<<" s"<<endl;

    cout<<"     Here We Transform CMP Gathers to Theta-T Gathers"<<endl;
    for(int it=0; it<nt; it++)
    {
        for(int i_theta=0; i_theta<n_theta; i_theta++)
        {
            float theta=(i_theta-n_theta/2)*d_theta;
            float theta_max_here=atan((nx/2)*dx*1.0/(it*dt*v/2));
            if(abs(theta)<theta_max_here)
            {
                //给定 theta t0，确定x
                float t0_temp=it*dt;
                float x=v*1.0/2*t0_temp*tan(theta);
                float x_coordinate_float=x*1.0/dx_interp;
                int x_coordinate=int(x_coordinate_float)+nx_interp/2;

                //给定 Theta t0，确定t
                float t=t0_temp*1.0/cos(theta);
                float t_coordinate_float_interp=t*1.0/dt_interp;
                int t_coordinate_interp=int(t_coordinate_float_interp);

                float t_coordinate_float=t*1.0/dt;
                int t_coordinate=int(t_coordinate_float);

                if(x_coordinate+1< nx_interp && x_coordinate>=0 &&
                        t_coordinate_interp+1< nt_interp && t_coordinate_interp>=0 &&
                        t_coordinate+1< nt && t_coordinate>=0)
                {
                    float a_x=(int(x_coordinate_float*1000000)%1000000)*1.0/1000000;
                    float b_x=1-a_x;

                    float a_t=(int(t_coordinate_float_interp*1000000)%1000000)*1.0/1000000;
                    float b_t=1-a_t;

                    seis_thetat(t_coordinate,i_theta)=
                        a_t*
                        (
                         a_x*seis_xt_interp(t_coordinate_interp,x_coordinate)+
                         b_x*seis_xt_interp(t_coordinate_interp,x_coordinate+1)
                        )
                        +b_t*
                        (
                         a_x*seis_xt_interp(t_coordinate_interp+1,x_coordinate)+
                         b_x*seis_xt_interp(t_coordinate_interp+1,x_coordinate+1)
                        );
                }
            }
        }
    }

    cout<<endl;
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Here We Begin to calculate the NMO Correction"<<endl;
    //动校正（NMO）
    cout<<"     Here We Synthesize the Zero_Offset Section"<<endl;
    //零偏移距剖面zero_off_CMP
    arma::Mat<float> zero_off(nt,nx,fill::zeros);
    for(int ix=0; ix<nx; ix++)
    {
        for(int iz=0; iz<nw; iz++)
        {
            zero_off(iz+nh-nw/2+1,ix)=wavelet(iz);
        }
    }

    //零偏移距剖面zero_theta
    arma::Mat<float> zero_theta(nt,n_theta,fill::zeros);
    for(int ix=0; ix<nx; ix++)
    {
        for(int iz=0; iz<nw; iz++)
        {
            float theta=(ix-n_theta/2)*d_theta;
            float theta_max_here=atan((n_theta/2)*dx*1.0/h);
            if(abs(theta)<theta_max_here)
            {
                zero_theta(iz+nh-nw/2+1,ix)=wavelet(iz);
            }
        }
    }


    cout<<"     Here We Begin to Calculate the NMO Correction of CMP Section"<<endl;
    //时间-空间域（CMP）动校正
    //参数说明：
    //      t_xt_norm ：按照时距关系公式计算的每一点的时间
    //      t_xt_nmo  ：动校正量
    //      seis_xt_after_nmo：动校正后的剖面
    //      seis_xt_op2：动校正后域零偏移距剖面做差结果
    arma::Mat<float> t_xt_norm(nt,nx,fill::zeros);
    arma::Mat<float> t_xt_nmo(nt,nx,fill::zeros);
    arma::Mat<float> seis_xt_after_nmo(nt,nx,fill::zeros);

    for(int ix=0; ix<nx; ix++)
    {
        for(int it=0; it<nt; it++)
        {
            float t0_temp=it*dt;
            float x_temp=abs((ix-(nx-1)/2)*dx);
            float t=sqrt(t0_temp*t0_temp+4.0*x_temp*x_temp*1.0/v*1.0/v);
            t_xt_norm(it,ix)=t;
            t_xt_nmo(it,ix)=t_xt_norm(it,ix)-t0_temp;
        }
    }

    for(int ix=0; ix<nx; ix++)
    {
        for(int it=0; it<nt; it++)
        {
            float t0_temp=it*dt;
            float t=t0_temp+t_xt_nmo(it,ix);
            float max_sup=t*1.0/dt;
            float a=(int(max_sup*1000000)%1000000)*1.0/1000000;
            float b=1-a;

            int it_ori=int(max_sup);
            if(it_ori+1<nt)
            {
                seis_xt_after_nmo(it,ix)=b*seis_xt(it_ori,ix)+a*seis_xt(it_ori+1,ix);
            }
        }
    }
    //op2
    arma::Mat<float> seis_xt_op2=seis_xt_after_nmo-zero_off;

    cout<<"     Here We Begin to Calculate the NMO Correction of T-Theta Section"<<endl;
    //时间-角度域动校正
    //参数说明：
    //      t_thetat_norm ：按照时距关系公式计算的每一点的时间
    //      t_thetat_nmo  ：动校正量
    //      seis_thetat_after_nmo：动校正后的剖面
    //      seis_thetat_op2：动校正后域零偏移距剖面做差结果
    arma::Mat<float> t_thetat_norm(nt,n_theta,fill::zeros);
    arma::Mat<float> t_thetat_nmo(nt,n_theta,fill::zeros);
    arma::Mat<float> seis_thetat_after_nmo(nt,n_theta,fill::zeros);

    for(int ix=0; ix<n_theta; ix++)
    {
        for(int it=1; it<nt; it++)
        {
            float theta=(ix-(nx-1)/2)*d_theta;
            float theta_max_here=atan(((nx-1)/2)*dx*1.0/(it*dt*v/2));
            if(abs(theta)<theta_max_here)
            {
                float h_temp=v*it*dt/2;
                float t0_temp=it*dt;
                float t=sqrt(t0_temp*t0_temp+4.0*h_temp*h_temp
                        *sin(theta)*sin(theta)*1.0/v*1.0/v
                        *1.0/cos(theta)*1.0/cos(theta));

                if(int(t*1.0/dt)<nt)
                {
                    t_thetat_norm(it,ix)=t;
                    t_thetat_nmo(it,ix)=t-t0_temp;
                }
            }
        }
    }

    for(int ix=0; ix<n_theta; ix++)
    {
        for(int it=0; it<nt; it++)
        {
            float theta=(ix-(nx-1)/2)*d_theta;
            float theta_max_here=atan(((nx-1)/2)*dx*1.0/(it*dt*v/2));
            if(abs(theta)<theta_max_here)
            {
                float t0_temp=it*dt;
                float t=t0_temp+t_thetat_nmo(it,ix);
                float max_sup=t*1.0/dt;
                float a=(int(max_sup*1000000)%1000000)*1.0/1000000;
                float b=1-a;

                int it_ori=int(max_sup);
                //cout<<it_ori<<endl;
                if(it_ori+1<nt)
                {
                    seis_thetat_after_nmo(it,ix)=b*seis_thetat(it_ori,ix)
                        +a*seis_thetat(it_ori+1,ix);
                }
            }
        }
    }

    //OP2
    arma::Mat<float> seis_thetat_op2=seis_thetat_after_nmo-zero_theta;

    //=========================================================================//
    //                              叠加成像                                   //
    //=========================================================================//
    //seis_thetat_stack :叠加数据   第一道表示零偏移距数据（即真实地震数据）
    //                              第二道表示叠加后数据
    //                              第三道为前两道相减
    //seis_xt_stack     :叠加数据   第一道表示零偏移距数据（即真实地震数据）
    //                              第二道表示叠加后数据
    //                              第三道为前两道相减
    arma::Mat<float> seis_thetat_stack(nt,3,fill::zeros);

    seis_thetat_stack.col(0)=seis_thetat.col(n_theta/2);
    for(int it=1; it<nt; it++)
    {
        for(int ix=0;ix<n_theta;ix++)
        {
            float theta=(ix-n_theta/2)*d_theta;
            float theta_max_here=atan(((nx-1)/2)*dx*1.0/(it*dt*v/2));
            if(abs(theta)<theta_max_here)
            {
                seis_thetat_stack(it,1)=seis_thetat_stack(it,1)+seis_thetat_after_nmo(it,ix);
            }
        }
        int n_theta_here=theta_max_here*1.0/d_theta*2.0;
        seis_thetat_stack(it,1)=seis_thetat_stack(it,1)/n_theta_here;
    }
    seis_thetat_stack.col(2)=seis_thetat_stack.col(0)-seis_thetat_stack.col(1);

    arma::Mat<float> seis_xt_stack(nt,3,fill::zeros);
    seis_xt_stack.col(0)=seis_xt.col(nx/2);
    for(int ix=0;ix<nx;ix++)
    {
        seis_xt_stack.col(1)=seis_xt_stack.col(1)+seis_xt_after_nmo.col(ix);
    }
    seis_xt_stack.col(1)=seis_xt_stack.col(1)/nx;
    seis_xt_stack.col(2)=seis_xt_stack.col(0)-seis_xt_stack.col(1);


    cout<<"--------------------------------------------------"<<endl;

    //=========================================================================//
    //                              文件保存                                   //
    //=========================================================================//
    cout<<"--------------------------------------------------"<<endl;
    cout<<"Save Files"<<endl;

    //文件保存路径 fn
    string fn="../file_DataBased/";

    string fn_seis_xt=fn+"seis_x_t.dat";
    seis_xt.save(fn_seis_xt,raw_binary);

    string fn_seis_theta=fn+"seis_theta_t.dat";
    seis_thetat.save(fn_seis_theta,raw_binary);

    string fn_zero_off=fn+"seis_x_t_zeros_off.dat";
    zero_off.save(fn_zero_off,raw_binary);

    string fn_zero_theta=fn+"seis_theta_t_zeros_theta.dat";
    zero_theta.save(fn_zero_theta,raw_binary);

    string fn_seis_xt_nmo=fn+"seis_x_t_nmo.dat";
    t_xt_nmo.save( fn_seis_xt_nmo,raw_binary);

    string fn_seis_xt_after_nmo=fn+"seis_x_t_after_nmo.dat";
    seis_xt_after_nmo.save(fn_seis_xt_after_nmo,raw_binary);

    string fn_seis_xt_op2=fn+"seis_x_t_op2.dat";
    seis_xt_op2.save(fn_seis_xt_op2,raw_binary);

    string fn_thetat_nmo=fn+"seis_theta_t_nmo.dat";
    t_thetat_nmo.save(fn_thetat_nmo,raw_binary);

    string fn_theta_after_nmo=fn+"seis_theta_t_after_nmo.dat";
    seis_thetat_after_nmo.save(fn_theta_after_nmo,raw_binary);

    string fn_thetat_op2=fn+"seis_theta_t_op2.dat";
    seis_thetat_op2.save(fn_thetat_op2,raw_binary);

    string fn_seis_thetat_stack=fn+"seis_theta_t_stack.dat";
    seis_thetat_stack.save(fn_seis_thetat_stack,raw_binary);

    string fn_seis_xt_stack=fn+"seis_x_t_stack.dat";
    seis_xt_stack.save(fn_seis_xt_stack,raw_binary);

    cout<<"--------------------------------------------------"<<endl;

    cout<<"--------------------------------------------------"<<endl;
    cout<<"FINISH"<<endl;
    finish = clock();
    cout<<"Time = "<<double(finish-start)/CLOCKS_PER_SEC<<"s"<<endl;
    cout<<"--------------------------------------------------"<<endl;


    return 0;
}


double rickerwavelet_func(int nw, arma::Col<float> & signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
    float t, t1;
    float temp[nw];
    float freq=28;
    float PI=3.1415926;
    float dt=0.001;

    for(i=0;i<nw;i++)
    {
        t=dt*i;
        t1=1.0/freq;//双边雷克子波
        //t1=0;//单边雷克子波
        temp[i]=10*(1-2*PI*PI*freq*freq*(t-t1)*(t-t1))
            *exp(-PI*PI*freq*freq*(t-t1)*(t-t1));
    }

    //导数形式
    /*    signal[0]=0;*/
    //for (i=1;i<param->Lw;i++)
    //{
    //signal[i]=(temp[i]-temp[i-1])*1.0/param->dt;
    /*    }*/

    //雷克子波形式
    for (i=0; i<nw; i++)
    {
        signal[i]=temp[i];
    }

    return 0.0;
}

//string = std::string;
int ReadPar_func(const string & fn_par, Parameter & par)
{
    auto Par = & par;
    //TJ_NMO_Ana_Par *Par = (TJ_NMO_Ana_Par *)calloc(1,sizeof(TJ_NMO_Ana_Par));
    int LineCount = 0;
    FILE *fp  = fopen(fn_par.c_str(),"r");
    if (fp==NULL)
    {
        printf("failed to open %s !!!\n",fn_par.c_str());
        return 0;
    }

    // read a complete line into line[]
    char line[TJU_LENGTH_FILENAME];
    while (1)
    {
        if (NULL==fgets (line,TJU_LENGTH_FILENAME,fp))
        {break;}
        if (line[0]=='*')
        {continue;}
        LineCount++;
        if (LineCount > TJU_NLINE_PAR)
        {break;}
        switch(LineCount){
            case 1 : sscanf(line,"%s",Par->fn_cmp_in );break;
            case 2 : sscanf(line,"%s",Par->fn_cmp_out);break;
            case 3 : sscanf(line,"%s",Par->fn_RadonSpectrum);break;
            case 4 : sscanf(line,"%s",Par->fn_tvpairs);break;
            case 5 : sscanf(line,"%d %f",&Par->nt,&Par->dt);break;
            case 6 : sscanf(line,"%d %f %f",&Par->nx,&Par->offset_min,&Par->doffset); break;
            case 7 : sscanf(line,"%d %f %f",&Par->nv,&Par->v_min,&Par->dv);break;
            case 8 : sscanf(line,"%f",&Par->freq_max);break;
            case 9 : sscanf(line,"%d",&Par->taper_length);break;
            case 10: sscanf(line,"%f %f",&Par->perc_under,&Par->perc_over);break;
            case 11: sscanf(line,"%d"   ,&Par->Nthd);break;
        }
    }

}

/*
 * Print Par
 */
int PrintPar_func(Parameter & par)
{
    auto Par = & par;
    printf("======================================================================\n");
    printf("***** PrintPar :   *****\n");
    // primary par
    printf("* 1  fn_cmp_in         = %s \n",Par->fn_cmp_in );
    printf("* 2  fn_cmp_out        = %s \n",Par->fn_cmp_out);
    printf("* 3  fn_RadonSpectrum  = %s \n",Par->fn_RadonSpectrum );
    printf("* 4  fn_tvpairs        = %s \n",Par->fn_tvpairs);
    printf("* 5  nt= %d \t,dt    = %f \n",Par->nt,Par->dt);
    printf("* 6  nx= %d \t,offmin= %f ,doff= %f\n",Par->nx,Par->offset_min,Par->doffset);
    printf("* 7  nv= %d \t,v_min = %f ,dv  = %f\n",Par->nv,Par->v_min,Par->dv);
    printf("* 8  freq_max       = %f \n",Par->freq_max);
    printf("* 9  taper_length   = %d \n",Par->taper_length);
    printf("* 10 perc_under=%f ,perc_over=%f\n",Par->perc_under,Par->perc_over);
    // secondary par
    printf("* Num of thread     = %d \n",Par->Nthd);

    printf("***** RadonDeMulti_PrintPar END *****\n");
    printf("======================================================================\n");
    return 0;
}

