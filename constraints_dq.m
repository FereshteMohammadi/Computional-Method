function Cq = constraints_dq(sys, q)

Cq = [
    constraints_dq_revolute(sys, q)
    constraints_dq_simple(sys, q)
    constraints_dq_simple_driving(sys, q)
    constraints_dq_translational(sys, q)
    ];