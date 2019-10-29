
def get_P1List(self, N=False):
        from sage.modular.modsym.p1list import P1List
    	if N is False:
		N = self.level * self.ideal_p
        return P1List(N)

def cusp_reduction_table(self):
    r'''
    Returns a dictionary and a set (of cusps).
    The dictionary keys are the normalized elements of P1(Z/N), say (c:d),
    and the value is a triple (c',d',T), where (c',d') is one of the cusp
    representatives (that is, is in the returned set), and T is an element
    in the stabilizer of Infinity that satisfies
    (c': d') * T = (c: d) (as elements of P1(Z/N)).
    '''
    from sage.modular.modsym.p1list import lift_to_sl2z
    P = self.get_P1List()
    remaining_points = set(list(P))
    reduction_table = dict([])
    cusp_set = set([])
    while len(remaining_points) > 0:
        c = remaining_points.pop()
        new_cusp = Matrix(ZZ,2,2,lift_to_sl2z(c[0], c[1], P.N()))
        new_cusp.set_immutable()
        cusp_set.add(new_cusp)
        reduction_table[c]=(new_cusp,matrix(ZZ,2,2,1))
        for hh in Zmod(P.N()):
            h = hh.lift()
            for u in [-1, 1]:
                new_c = P.normalize(u * c[0], u**-1 * c[1] + h * c[0])
                if new_c not in reduction_table:
                    remaining_points.remove(new_c)
                    T = matrix(ZZ,2,2,[u,h,0,u**-1])
                    reduction_table[new_c]=(new_cusp, T)
                    assert P.normalize(*(vector(c) * T)) == new_c
    return reduction_table, cusp_set

def find_matrix_from_cusp(self, cusp):
    r'''
    Returns a matrix gamma and a cusp representative modulo Gamma0(N) (c2:d2),
    such that gamma * cusp = (c2:d2).
    '''
    from sage.modular.modsym.p1list import lift_to_sl2z
    a, c = cusp.numerator(), cusp.denominator()
    reduction_table, _ = self.cusp_reduction_table()
    P = self.get_P1List()

    # Find a matrix g = [a,b,c,d] in SL2(Z) such that g * a/c = oo
    # Define (c1:d1) to be the rep in P1(Z/N) such that (c1:d1) == (c:d).
    if c == 0:
        assert a == 1
        g = Matrix(ZZ,2,2,[1,0,0,1])
        c1, d1 = P.normalize(0, 1)
    else:
        g0, d, b = ZZ(a).xgcd(-c)
        assert g0 == 1
        g = Matrix(QQ,2,2,[d,-b,-c,a]) # the inverse
        c1, d1 = P.normalize(c, d)
    assert g.determinant() == 1

    A, T = reduction_table[(c1,d1)]
    gamma = A.parent()(A * T * g)

    # This block test correctness
    tst = Cusp(Gamma0(P.N())(gamma).acton(Cusp(a,c)))
    tst = (tst.numerator(), tst.denominator())
    assert tst == tuple(A.column(0))
    return gamma, A
