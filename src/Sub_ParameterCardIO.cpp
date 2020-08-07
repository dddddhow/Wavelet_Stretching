//----------------------------------------------------------------------------//
//程序目的： 1、参数卡读写函数
//
//程序原理：
//
//
//程序参数说明：
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

int readpar_func(const string &fn_par, Parameter &par)
{
    auto Par = &par;
    int LineCount = 0;
    FILE *fp = fopen(fn_par.c_str(), "r");
    if (fp == NULL)
    {
        printf("failed to open %s !!!\n", fn_par.c_str());
        return 0;
    }

    // read a complete line into line[]
    char line[TJU_LENGTH_FILENAME];
    while (1)
    {
        if (NULL == fgets(line, TJU_LENGTH_FILENAME, fp))
        {
            break;
        }
        if (line[0] == '*')
        {
            continue;
        }
        LineCount++;
        if (LineCount > TJU_NLINE_PAR)
        {
            break;
        }
        switch (LineCount)
        {
            case 1:
                sscanf(line, "%s", Par->fn_path_of_model_in);
                break;
            case 2:
                sscanf(line, "%s", Par->fn_path_of_v_rms_in);
                break;
            case 3:
                sscanf(line, "%s", Par->fn_path_of_output);
                break;
            case 4:
                sscanf(line, "%d  %d %f %f", &Par->nvx, &Par->nvz, &Par->dx, &Par->dz);
                break;
            case 5:
                sscanf(line, "%d %d %f", &Par->nx, &Par->nt, &Par->dt);
                break;
            case 6:
                sscanf(line, "%d", &Par->nw);
                break;
        }
    }

    fclose(fp);
    return 0;
}




/*
 * Print Par
 */
int printpar_func(Parameter &par)
{
    auto Par = &par;
    printf("======================================================================\n");
    printf("PrintPar : \n");
    // primary par
    printf("* 1  Path of model     : %s \n", Par->fn_path_of_model_in);
    printf("* 2  Path of V_rms     : %s \n", Par->fn_path_of_v_rms_in);
    printf("* 3  Path of Output    : %s \n", Par->fn_path_of_output);
    printf("* 5  Model Size        : (nvz,nvx) -> (%d,%d)\n", Par->nvz, Par->nvx);
    printf("* 6  V_rms Size        : (nt,nx)   -> (%d,%d)\n", Par->nt, Par->nx);
    printf("* 7  Sample Interval   : (dz,dx,dt)-> (%f m,%f m,%f s)\n",Par->dz, Par->dx, Par->dt);
    printf("* 8  Number of Wavelet : %d \n", Par->nw);
    printf("======================================================================\n");
    return 0;
}

