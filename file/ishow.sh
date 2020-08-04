n1=2000
n2=201
d1=0.001
d2=5
f1=0
f2=-500


xwigb <seis_x_t.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-OFFSET  xbox=100&
xwigb <seis_x_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-OFFSET-Nmo xbox=700&
xwigb <seis_x_theta_xt_after_nmo.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=X-THETA-XT-NMO xbox=1300 &
xwigb <seis_x_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=2100&

ximage <seis_x_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-OFFSET_NMO xbox=100 ybox=1000 legend=1&
ximage <seis_x_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=700 ybox=1000 legend=1&

xwigb <seis_x_t_stack.dat n1=$n1 d1=$d1 d2=$d2 title=X-T-STACK xbox=1300 ybox=1000 &
