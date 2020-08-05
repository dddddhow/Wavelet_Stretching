//----------------------------------------------------------------------------//
//程序目的： 1、生成子波(Rickr子波)
//
//程序原理：
//
//程序参数说明：
//      nw：子波长度
//      dt: 时间采样间隔（秒）
//      signal：子波（向量）
//
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

#include <armadillo>
#include <iostream>
#include "../lib/TJWPI_ShengShen.h"

using namespace arma;

double rickerwavelet_func(int nw, float dt, arma::Col<float> &signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
    float t, t1;
    float temp[nw];
    float freq = 28;
    float PI = 3.1415926;

    for (i = 0; i < nw; i++)
    {
        t = dt * i;
        t1 = 1.0 / freq; //双边雷克子波
        //t1=0;//单边雷克子波
        temp[i] = 10 * (1 - 2 * PI * PI * freq * freq * (t - t1) * (t - t1)) *
            exp(-PI * PI * freq * freq * (t - t1) * (t - t1));
    }

    //雷克子波形式
    for (i = 0; i < nw; i++)
    {
        signal[i] = temp[i];
    }

    return 0.0;
}

