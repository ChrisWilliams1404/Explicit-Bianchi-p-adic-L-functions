
class arith_gp(object):
        
    def __init__(self,N):
        
        self.N = N
        self._get_P1List = get_P1List(None,N)
        self._cusp_reduction_table,self._cusp_set = cusp_reduction_table(self._get_P1List)

    def get_P1List(self):
        return self._get_P1List

    def cusp_reduction_table(self):
        return self._cusp_reduction_table

    def cusp_set(self):
        return self._cusp_set
       
    def find_matrix_from_cusp(self, cusp = None, a=1,c=0):
        r'''
        Returns a matrix gamma and a cusp representative modulo Gamma0(N) (c2:d2),
        represented as a matrix (a,b;c,d), such that gamma * cusp = (c2:d2).

        HACK: inputting a,c manually as we don't have cusp package
        '''

        ## Initialise a,c if we have cusp to input
        if cusp is not None:            
            a, c = cusp.numerator(), cusp.denominator()

        reduction_table = self.cusp_reduction_table()
        P = self.get_P1List()
        if hasattr(P.N(),'number_field'):
            K = P.N().number_field()
        else:
            K = QQ

        # Find a matrix g = [a,b,c,d] in SL2(O_K) such that g * a/c = oo
        # Define (c1:d1) to be the rep in P1(O_K/N) such that (c1:d1) == (c:d).
        if c == 0: ## case cusp infinity: (a,c) should equal (1,0)
            assert a == 1
            g = Matrix(2,2,[1,0,0,1])
            c1, d1 = P.normalize(0, 1)
        else:
            if K == QQ:
                g0, d, b = ZZ(a).xgcd(-c)
                assert g0 == 1
            else:
                d,b = xgcd_F(a,-c)

            g = Matrix(2,2,[d,-b,-c,a]) # the inverse
            c1, d1 = P.normalize(c, d)
        assert g.determinant() == 1

        A, T = reduction_table[(c1,d1)]
        gamma = A.parent()(A * T * g)

        # This block test correctness
        #tst = Cusp(Gamma0(P.N())(gamma).acton(Cusp(a,c)))
        #tst = (tst.numerator(), tst.denominator())
        #assert tst == tuple(A.column(0))
        return gamma, A



def get_P1List(self, N=None):
    """
    Generates the projective line of O_F/N, where N is an ideal specified   
    in the input, or computed from a parent object (e.g. arithmetic group).
    """

    ## If N is not specified, then compute it from self
    if N is None:
        N = self.level * self.ideal_p

    ## Return object representing Projective line over O_F/N
    if hasattr(N,'number_field'): ## Base field not Q   
        from sage.modular.modsym.p1list_nf import P1NFList
        return P1NFList(N)
    else:	## Base field Q
        from sage.modular.modsym.p1list import P1List
        return P1List(N)
		


def cusp_reduction_table(P):
    r'''
    Returns a dictionary and the set of cusps.

    Assumes we have a finite set surjecting to the cusps (namely, P^1(O_F/N)). Runs through
    and computes a subset which represents the cusps, and shows how to go from any element 
    of P^1(O_F/N) to the chosen equivalent cusp.

    Takes as input the object representing P^1(O_F/N), where F is a number field
    (that is possibly Q), and N is some ideal in the field.  Runs the following algorithm:
	    - take a remaining element C = (c:d) of P^1(O_F/N);
	    - add this to the set of cusps, declaring it to be our chosen rep;
	    - run through every translate C' = (c':d') of C under the stabiliser of infinity, and
		    remove this translate from the set of remaining elements;
	    - store the matrix T in the stabiliser such that C' * T = C (as elements in P^1)
		    in the dictionary, with key C'.
    '''
    if hasattr(P.N(),'number_field'):
        K = P.N().number_field()
    else:
        K = QQ
   
    from sage.modular.modsym.p1list_nf import lift_to_sl2_Ok
    from sage.modular.modsym.p1list import lift_to_sl2z
    ## Define new function on the fly to pick which of Q/more general field we work in
    ## lift_to_matrix takes parameters c,d, then lifts (c:d) to a 2X2 matrix over the NF representing it
    lift_to_matrix = lambda c, d: lift_to_sl2z(c,d,P.N()) if K.degree() == 1 else lift_to_sl2_Ok(P.N(), c, d)

    ## Put all the points of P^1(O_F/N) into a list; these will corr. to our dictionary keys
    remaining_points = set(list(P)) if K == QQ else set([c.tuple() for c in P])
    reduction_table = dict([])
    cusp_set = list([])

    ## Loop over all points of P^1(O_F/N)
    while len(remaining_points) > 0:
        ## Pick a new cusp representative
        c = remaining_points.pop()
        ## c is an MSymbol so not hashable. Create tuple that is
        ## Represent the cusp as a matrix, add to list of cusps, and add to dictionary
        new_cusp = Matrix(2,2,lift_to_matrix(c[0], c[1])) 
        new_cusp.set_immutable()
        cusp_set.append(new_cusp)
        reduction_table[c]=(new_cusp,matrix(2,2,1)) ## Set the value to I_2
        ## Now run over the whole orbit of this point under the stabiliser at infinity.
        ## For each elt of the orbit, explain how to reduce to the chosen cusp.

        ## Run over lifts of elements of O_F/N:
        if K == QQ:
            residues = Zmod(P.N())
            units = [1, -1]
        else:
            residues = P.N().residues()
            units = K.roots_of_unity()

        for hh in residues:
            h = K(hh) ## put into the number field
            ## Run over all finite order units in the number field
            for u in units:
                ## Now have the matrix (u,h; 0,u^-1).
                ## Compute the action of this matrix on c
                new_c = P.normalize(u * c[0], u**-1 * c[1] + h * c[0])
                if K != QQ: 
                    new_c = new_c.tuple()
                if new_c not in reduction_table:
                    ## We've not seen this point before! But it's equivalent to c, so kill it!
                    ## (and also store the matrix we used to get to it)
                    remaining_points.remove(new_c)
                    T = matrix(2,2,[u,h,0,u**-1]) ## we used this matrix to get from c to new_c
                    reduction_table[new_c]=(new_cusp, T) ## update dictionary with the new_c + the matrix
                    if K != QQ:
                        assert P.normalize(*(vector(c) * T)).tuple() == new_c ## sanity check
                    else:
                        assert P.normalize(*(vector(c) * T)) == new_c ## sanity check


    return reduction_table, cusp_set



def xgcd_F(a,c):
    """
    Compute gcd if a,c are coprime in F, and x,y such that
        ax+cy = 1.
    """
    if a.parent() != c.parent():
        raise ValueError('a,c not in the same field.')
    else:
        OK = a.parent()
    if gcd(a,c) != 1:    
        raise ValueError('a,c not coprime.')

    for xbar in OK.ideal(c).residues():
        if a*xbar - 1 in OK.ideal(c):
            x = xbar
            break

    y = (1 - a*x)/c
    return x,y




X = var('X')
K.<a> = NumberField(X^2+1)
N = K.primes_above(13)[0]
G = arith_gp(N)
H = arith_gp(5)