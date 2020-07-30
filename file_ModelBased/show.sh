
n1=2000
n2=101
d1=0.0005
d2=1.36397
f1=0
f2=-68.1986

xwigb <seis_theta_t.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA  xbox=10 &
xwigb <seis_theta_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA-NMO xbox=500 &
xwigb <seis_theta_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=1000 &

ximage <seis_theta_t_after_nmo.dat  n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=T-THETA-NMO xbox=500 ybox=100 legend=1&
ximage <seis_theta_t_op2.dat n1=$n1 d1=$d1 d2=$d2 f1=$f1 f2=$f2 title=OP2 xbox=1000 ybox=100 legend=1&
xwigb <seis_theta_t_stack.dat n1=$n1 d1=$d1 d2=$d2 title=Theta-T-STACK xbox=10 ybox=100  x1beg=0.1 x1end=0.3 &
