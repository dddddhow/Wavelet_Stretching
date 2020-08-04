set -x

n1=2000
n2=201
d1=0.001
f1=0
f2x=-500
f2theta=-80
d2x=5
d2theta=0.8


xwigb <seis_ori_stack_all.dat n1=$n1 d1=$d1 d2=$d2x f1=$f1 f2=$f2x title=Ori  xbox=100&
xwigb <seis_x_t_stack_all.dat  n1=$n1 d1=$d1 d2=$d2x f1=$f1 f2=$f2x title=X_T xbox=700&
xwigb <seis_theta_t_stack_all.dat n1=$n1 d1=$d1 d2=$d2theta f1=$f1 f2=$f2theta title=Theta_T xbox=1300&

