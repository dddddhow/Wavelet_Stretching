#define TJU_NLINE_PAR           11      // num of effective lines in par file
#define TJU_LENGTH_FILENAME     300

typedef struct Parameter
{
    //====================== Start Of Par File =================================//

    char    fn_cmp_in [TJU_LENGTH_FILENAME];            // line 01
    char    fn_cmp_out[TJU_LENGTH_FILENAME];            // line 02
    char    fn_RadonSpectrum[TJU_LENGTH_FILENAME];      // line 03
    char    fn_tvpairs[TJU_LENGTH_FILENAME];            // line 04
    int     nt;  float   dt;                            // line 05
    int     nx;  float   offset_min;  float   doffset;  // line 06
    float   *XArray;                       // extra
    int     nv;  float   v_min     ;  float   dv     ;  // line 07
    float   freq_max;                                   // line 08
    int     taper_length;                               // line 09
    float   perc_under,perc_over;                       // line 10
    int     Nthd;                                       // line 11
    //====================== End   Of Par File =================================//

} Parameter;
