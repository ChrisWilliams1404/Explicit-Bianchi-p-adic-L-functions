## Specific functionality to work with overconvergent Bianchi modular symbols
# from sarithgroup import BigArithGroup
# page_path = '/home/float/darmonpoints/darmonpoints/KleinianGroups-1.0/klngpspec'
# magma.attach_spec(page_path)
## Define number field and prime

K.<a> = QuadraticField(-1)
E = EllipticCurve('11a1')
assert len(K.ideal(p).factor()) == 2
P = K.ideal(p).factor()[0][0]
N = K.ideal(p * E.conductor())
assert K.ideal(p).divides(N)
M = N/P
G = BigArithGroup(P, (1,1), M, base= K, magma = magma, use_shapiro=True,grouptype="PGL2") # needs magma

K.<a>  = QuadraticField(-11)
P = K.ideal(5).factor()[0][0]
G = BigArithGroup(P, (1,1), K.ideal(5).factor()[1][0] * K.ideal(11).factor()[0][0], base= K, magma = magma, use_shapiro=False) # needs magma

from cohomology_arithmetic import * # needs darmonpoints
HH = ArithCoh(G) # needs darmonpoints & magma
phi = HH.get_cocycle_from_elliptic_curve(EllipticCurve('55a').change_ring(K)) # needs darmonpoints & magma

