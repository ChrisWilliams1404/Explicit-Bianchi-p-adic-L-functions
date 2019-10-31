from darmonpoints.sarithgroup import *
from darmonpoints.cohomology_arithmetic import *
from darmonpoints.homology import *
from sage.modular.pollack_stevens.manin_map import unimod_matrices_to_infty
from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial
from darmonpoints.integrals import integrate_H1, get_basic_integral
from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial

########################
p = 11
working_prec = 10
prec = working_prec
use_ps_dists = False
use_shapiro = False

E = EllipticCurve('11a1')

G = BigArithGroup(p, (1,1), E.conductor() / p, base=QQ,magma=magma,use_shapiro=False,matrix_group=True)
G1 = BigArithGroup(p, (1,1), E.conductor() / p, base=QQ,magma=magma,use_shapiro=False,matrix_group=False)

HH = ArithCoh(G, use_ps_dists=use_ps_dists)
HH1 = ArithCoh(G1, use_ps_dists=use_ps_dists)
phi = HH.get_cocycle_from_elliptic_curve(E,0)
phi1 = HH.get_cocycle_from_elliptic_curve(E,0)
Phi = get_overconvergent_class_quaternionic(p,phi,G,prec,0,E.ap(p),use_ps_dists=use_ps_dists)
Phi1 = get_overconvergent_class_quaternionic(p,phi1,G1,prec,0,E.ap(p),use_ps_dists=use_ps_dists)

## Debug parabolic property
Phi2 = Phi.parent()(Phi1.values())
eps = Phi - Phi2
HOC = Phi.parent()
Phi_cusp = eps
Gp = HOC.group()
g = Gp(Matrix(ZZ,2,2,[1,0,-11,1]))
A = HOC.coefficient_module().acting_matrix(g,prec+1).change_ring(Qp(p,prec))-1
b = Phi_cusp.evaluate(g)._moments.change_ring(Qp(p,prec))
print A.solve_right(b,check=False)

## Check difference is a couboundary
Phi2 = Phi.parent()(Phi1.values())
eps = Phi - Phi2
HOC = Phi.parent()
Phi_cusp = eps
Gp = HOC.group()
A = Matrix(Qp(p,prec),0,prec+1,0)
b = Matrix(Qp(p,prec),0,1,0)
for g in Gp.gens():
    A = A.stack(HOC.coefficient_module().acting_matrix(g,prec+1).change_ring(Qp(p,prec))-1)
    b = b.stack(Phi_cusp.evaluate(g)._moments.change_ring(Qp(p,prec)))
v = HOC.coefficient_module()(A.solve_right(b,check=False))
for g in Gp.gens():
    print (eps.evaluate(g) - (g * v - v)).valuation_list()


# # Check that the modular symbol is a lift of the corresponding cohomology class
# v = Phi1.cusp_boundary_element((1,0))
# for g0 in Gp.gens():
#     g = g0.quaternion_rep
#     cusp_list = [(1,(1,0)), (-1, (1,0).apply(g.list()))]
#     A = Phi1.evaluate_cuspidal_modsym_at_cusp_list(cusp_list)
#     B = Phi1.evaluate(g) + v - g * v
#     print (A-B).valuation_list()


### Pollack's approach
L = E.padic_lseries(p,implementation="pollackstevens",precision=prec)

# Compare Lseries
print '####################'
for n in range(1,5):
    A = L[n]
    B = 1/2*Phi1.get_Lseries_term(n)
    print A == B
    print '..'
print '####################'


##########################################################################
## Second Test
#######################################################################
from darmonpoints.sarithgroup import *
from darmonpoints.cohomology_arithmetic import *
from darmonpoints.homology import *
from sage.modular.pollack_stevens.manin_map import unimod_matrices_to_infty
from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial
from darmonpoints.integrals import integrate_H1, get_basic_integral

p = 5
working_prec = 10
prec = working_prec
use_ps_dists = False
use_shapiro = False
set_verbose(0)

E = EllipticCurve('15a1')

G = BigArithGroup(p, (1,1), E.conductor() / p, base=QQ,magma=magma,use_shapiro=False, matrix_group=False)
HH = ArithCoh(G, use_ps_dists=use_ps_dists)
phi = HH.get_cocycle_from_elliptic_curve(E,0)
Phi = get_overconvergent_class_quaternionic(p,phi,G,prec,0,E.ap(p),use_ps_dists=use_ps_dists)
Phi.elliptic_curve = E

### Pollack's approach
L = E.padic_lseries(p, implementation='pollackstevens',precision=prec)

# Compare Lseries
print '####################'
for n in range(1,10):
    %time A = 2 * L[n]
    %time B = Phi.get_Lseries_term(n)
    print A == B
    print '..'
print '####################'

