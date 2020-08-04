set -x
n1=2000
n2=3
x1end=150
f1=0
d1=0.5
width=1000
height=200
linewidth=2

xgraph <seis_ori_x_t_theta_t_fourier_abs.dat n1=$n1 d1=$d1 n2=$n2 width=$width height=$height x1end=$x1end label1=Frequence label2=Amplitude windowtitle=Ori_XT_ThetaT_FTAna style=normal  linewidth=$linewidth  x2end=200&

xgraph <seis_x_t_stack_fourier_abs.dat n1=$n1 d1=$d1 width=$width height=$height x1end=$x1end label1=Frequence label2=Amplitude windowtitle=X_T style=normal  linewidth=$linewidth linecolor=3  x2end=200 &


xgraph <seis_theta_t_stack_fourier_abs.dat n1=$n1 d1=$d1 width=$width height=$height x1end=$x1end label1=Frequence label2=Amplitude windowtitle=Theta_T style=normal  linewidth=$linewidth  linecolor=4   x2end=200 &

xgraph <seis_ori_fourier_abs.dat n1=$n1 d1=$d1 width=$width height=$height x1end=$x1end label1=Frequence label2=Amplitude windowtitle=Ori style=normal  linewidth=$linewidth linecolor=2   x2end=200&

