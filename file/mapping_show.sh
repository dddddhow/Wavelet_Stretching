n1=2000
n2=201
d1=0.001
d2theta=0.8
d2x=5
f1=0
f2theta=-80
f2x=-500

xwigb <seis_x_t.dat                  n1=$n1 d1=$d1 d2=$d2x     f1=$f1 f2=$f2x     title=XT_ModelBased         xbox=100 ybox=100&
xwigb <seis_x_t_mappedfromthetat.dat n1=$n1 d1=$d1 d2=$d2x     f1=$f1 f2=$f2x     title=XT_MappedTFromThetaT  xbox=700 ybox=100&

xwigb <seis_theta_t_modelbased.dat   n1=$n1 d1=$d1 d2=$d2theta f1=$f1 f2=$f2theta title=ThetaT_ModelBased     xbox=100  ybox=1000&
xwigb <seis_theta_t.dat              n1=$n1 d1=$d1 d2=$d2theta f1=$f1 f2=$f2theta title=ThetaT_MappedFromXT   xbox=700  ybox=1000&
