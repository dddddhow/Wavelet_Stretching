
n1=2000
n2=201
d1=0.001
d2=0.8
f1=0
f2=-80

xwigb <seis_theta_t.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA  xbox=100 &
xwigb <seis_theta_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA-NMO xbox=700 &
#xwigb <seis_x_theta_xt_after_nmo.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=X-THETA-XT-NMO xbox=1300 &
xwigb <seis_theta_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=1300 &

ximage <seis_theta_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA-NMO xbox=100 ybox=700 legend=1&
ximage <seis_theta_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=700 ybox=700 legend=1&

xwigb <seis_theta_t_stack.dat n1=$n1 d1=$d1 d2=$d2 title=Theta-T-STACK xbox=1300 ybox=700 &
