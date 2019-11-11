#### COMPUTING p-ADIC L-FUNCTIONS ##########
from darmonpoints.sarithgroup import *
from darmonpoints.cohomology_arithmetic import *
from darmonpoints.homology import *
from sage.modular.pollack_stevens.manin_map import unimod_matrices_to_infty
from darmonpoints.integrals import integrate_H1, get_basic_integral

########################
p = 11 # prime to work with
working_prec = 10 # precision of the p-adics
prec = working_prec # the precision of the distributions
use_ps_dists = False # whether to use Pollack's implementation of distributions.
use_shapiro = False # wheter to use Shapiro to work with a group without p in the level

E = EllipticCurve('11a1')

# The following initializes the S-arithmetic group data.
# INPUT:
# - p : prime
# - (a,b): invariants of quaternion algebra. For Matrix algebra, input (1,1)
# - level: the extra level. This must be coprime to p.
# - base: the base field
G = BigArithGroup(p, (1,1), E.conductor() / p, base=QQ, magma=magma, use_shapiro=use_shapiro)

# Create the cohomology group (with trivial coefficients). Only need to pass the group G created before.
HH = ArithCoh(G)

# Get a cohomology class attached to the elliptic curve the '0' parameter means: take phi^+ + phi^-.
phi = HH.get_cocycle_from_elliptic_curve(E,0)

# Lift the cocycle to an overconvergent class. Note that one needs to give the eigenvalue of Up (obtained via E.ap(p)
Phi = get_overconvergent_class_quaternionic(p,phi,G,prec,0,E.ap(p),use_ps_dists=use_ps_dists)

### Pollack's approach - from vanilla Sage
L = E.padic_lseries(p,implementation="pollackstevens",precision=prec)

# We compare the L-series obtained from Sage's approach with ours. Note the factor of 1/2...
print '####################'
for n in range(1,5):
    A = L[n]
    B = 1/2 * Phi.get_Lseries_term(n)
    print A == B
print '####################'


##########################################################################
## Second Test
#######################################################################
from darmonpoints.sarithgroup import *
from darmonpoints.cohomology_arithmetic import *
from darmonpoints.homology import *
from sage.modular.pollack_stevens.manin_map import unimod_matrices_to_infty
from darmonpoints.integrals import integrate_H1, get_basic_integral

p = 5
working_prec = 10
prec = working_prec
use_ps_dists = False
use_shapiro = False
set_verbose(0)

E = EllipticCurve('15a1')

G = BigArithGroup(p, (1,1), E.conductor() / p, base=QQ, magma=magma, use_shapiro=False)
HH = ArithCoh(G, use_ps_dists=use_ps_dists)
phi = HH.get_cocycle_from_elliptic_curve(E,0)
Phi = get_overconvergent_class_quaternionic(p,phi,G,prec,0,E.ap(p),use_ps_dists=use_ps_dists)
Phi.elliptic_curve = E

### Pollack's approach
L = E.padic_lseries(p, implementation='pollackstevens',precision=prec)

# Compare Lseries
print '####################'
for n in range(1,10):
    A = 2 * L[n]
    B = Phi.get_Lseries_term(n)
    print A == B
print '####################'

##########################################################################
## Bianchi Test
#######################################################################
from darmonpoints.sarithgroup import BigArithGroup
from darmonpoints.cohomology_arithmetic import *

page_path = '~/Dropbox/darmonpoints/darmonpoints/KleinianGroups-1.0/klngpspec'
magma.eval('SetSeed(%s)'%1231241)
magma.attach_spec(page_path)
magma.eval('Page_initialized := true')
set_verbose(1)

K.<a> = QuadraticField(-1)
E = EllipticCurve('11a1')
p = 5
assert len(K.ideal(p).factor()) == 2
P = K.ideal(p).factor()[0][0]
N = K.ideal(p * E.conductor())
assert K.ideal(p).divides(N)
M = N/P

implementation = 'coset_enum'
G = BigArithGroup(P, (1,1), M, base= K, magma = magma, use_shapiro=True,grouptype="PGL2", implementation=implementation) # needs magma

HH = ArithCoh(G) # needs darmonpoints & magma
phi = HH.get_cocycle_from_elliptic_curve(EllipticCurve('55a').change_ring(K)) # needs darmonpoints & magma

prec = 10
ap = ##
apbar = ##
Phi = get_overconvergent_class_bianchi(P,phi,G,prec,ap, apbar)
Phi.elliptic_curve = E

r, s = 1,1
print Phi.get_Lseries_term((r,s))
