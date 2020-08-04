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

int ReadPar_func(const string &fn_par, Parameter &par)
{
    auto Par = &par;
    //TJ_NMO_Ana_Par *Par = (TJ_NMO_Ana_Par *)calloc(1,sizeof(TJ_NMO_Ana_Par));
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
                sscanf(line, "%s", Par->fn_cmp_in);
                break;
            case 2:
                sscanf(line, "%s", Par->fn_cmp_out);
                break;
            case 3:
                sscanf(line, "%s", Par->fn_RadonSpectrum);
                break;
            case 4:
                sscanf(line, "%s", Par->fn_tvpairs);
                break;
            case 5:
                sscanf(line, "%d %f", &Par->nt, &Par->dt);
                break;
            case 6:
                sscanf(line, "%d %f %f", &Par->nx, &Par->offset_min, &Par->doffset);
                break;
            case 7:
                sscanf(line, "%d %f %f", &Par->nv, &Par->v_min, &Par->dv);
                break;
            case 8:
                sscanf(line, "%f", &Par->freq_max);
                break;
            case 9:
                sscanf(line, "%d", &Par->taper_length);
                break;
            case 10:
                sscanf(line, "%f %f", &Par->perc_under, &Par->perc_over);
                break;
            case 11:
                sscanf(line, "%d", &Par->Nthd);
                break;
        }
    }
}

/*
 * Print Par
 */
int PrintPar_func(Parameter &par)
{
    auto Par = &par;
    printf("======================================================================\n");
    printf("***** PrintPar :   *****\n");
    // primary par
    printf("* 1  fn_cmp_in         = %s \n", Par->fn_cmp_in);
    printf("* 2  fn_cmp_out        = %s \n", Par->fn_cmp_out);
    printf("* 3  fn_RadonSpectrum  = %s \n", Par->fn_RadonSpectrum);
    printf("* 4  fn_tvpairs        = %s \n", Par->fn_tvpairs);
    printf("* 5  nt= %d \t,dt    = %f \n", Par->nt, Par->dt);
    printf("* 6  nx= %d \t,offmin= %f ,doff= %f\n", Par->nx, Par->offset_min, Par->doffset);
    printf("* 7  nv= %d \t,v_min = %f ,dv  = %f\n", Par->nv, Par->v_min, Par->dv);
    printf("* 8  freq_max       = %f \n", Par->freq_max);
    printf("* 9  taper_length   = %d \n", Par->taper_length);
    printf("* 10 perc_under=%f ,perc_over=%f\n", Par->perc_under, Par->perc_over);
    // secondary par
    printf("* Num of thread     = %d \n", Par->Nthd);

    printf("***** RadonDeMulti_PrintPar END *****\n");
    printf("======================================================================\n");
    return 0;
}

