##============================================================
##
##              PRESENTING ARITHMETIC GROUPS
##              
##============================================================


class arith_gp(object):
        
    def __init__(self,N):
        
        self.N = N
        self._get_P1List = get_P1List(None,N)
        self._cusp_reduction_table,self._cusp_set = cusp_reduction_table(self._get_P1List)
        self.P = self._get_P1List
        if hasattr(self.P.N(),'number_field'):
            self.K = self.P.N().number_field()
        else:
            self.K = QQ

        ## HARD CODING SL2(Z) FOR NOW SINCE THIS IS ALREADY IN SAGE
        self.level_1_gp = Gamma0(1) ## Insert your favourite group here
        self.level_1_word_problem = self.level_1_gp.farey_symbol()
        self.level_1_gens = self.level_1_word_problem.generators() ## WARNING!!! FAREY_SYMBOL USES DIFFERENT GENS TO GP.gens()!!!!
       
        
    ##=================================================
    ## BASIC FUNCTIONALITY
    def __repr__(self):
        return 'Class for working with arithmetic group with index {}'.format(self.get_P1List())

    def get_P1List(self):
        return self._get_P1List

    def cusp_reduction_table(self):
        return self._cusp_reduction_table

    def cusp_set(self):
        return self._cusp_set


    ##=================================================
    ## COMPUTE GENS OF SUBGROUP

    @cached_method
    def coset_reps(self):
        """
        Compute generators of the subgroup Gamma_0(N), where N is the specified level.
    
        Write down representatives of the cosets
        for Gamma_0(N) in Gamma(1), which we identify with P^1(O_F/N). We already have
        code to compute with this: namely, cusp_reduction_table does precisely this.
        """
        ## Retrieve the cusp reduction table. Recall that this is a dictionary with keys
        ## given by tuples (a,c) representing the element (a:c) in P^1(O_F/N). The entries
        ## are C, A, where c is the corresponding cusp (from cusp_set) and A is a matrix 
        ## taking C to a matrix with bottom row (a:c). In particular: 
        crt = self.cusp_reduction_table()

        ## Generate the coset representatives: this is given by taking A*C as we range 
        ## over all (A,C) in the values of cusp_reduction_table
        coset_reps = {}
        for key in crt.keys():
            coset_reps[key] = crt[key][0] * crt[key][1]

        ## coset_reps is now a dictionary: keys are elements of P^1(O_F/N), and values are
        ## matrices which are coset reps for Gamma_0(N) in Gamma(1) cor. to these elements

        return coset_reps

    
    def represent_in_coset(self,g):
        """
        g is an element of Gamma(1). Represent it as h.p, where h in Gamma_0 and p is a rep.
        """
        ## We can read off the class from the bottom row, computing in P(O_F/N)
        c,d = g[1]
        if self.K == QQ:
            coset_class = self.P.normalize(c,d)
        else:
            coset_class = self.P.normalize(c,d).tuple()
        representative = self.coset_reps()[coset_class]
        h = g*representative^(-1)
        
        ## Now check that h really is in Gamma_0(N)
        if self.K == QQ:
            assert h[1][0] in ZZ.ideal(self.P.N())
        else:
            assert h[1][0] in self.P.N()
        
        return (h,representative)

    @cached_method
    def compute_generators(self):
        """
        Compute generators for Gamma_0(N) via its right coset representatives in Gamma(1).

        Returns:
            - small_gens_matrices, a dictionary: the keys are matrices A, 
                which are generators of the small group, and the values are integers;

            - small_gens_words, a list: the D[A]-th entry is the matrix A written 
                as a word in the generators of Gamma(1). 

        The words are written in Tietze form, i.e. [1,2,1,-2,-1,2] corr. to 
        g * h * g * h^(-1) * g^(-1) * h, where g,h = (ordered) gens of Gamma(1).

        """
        big_gp_gens = tuple(self.level_1_gens) + tuple([~g for g in self.level_1_gens]) 
        small_gens_matrices_dict = {} ## This will contain the output: matrix form, dictionary with keys the matrices
        small_gens_matrices = [] ## list of the matrices in order
        small_gens_words = [] ## This will contain the output: word form. 
        current_index = 0

        ## Loop over all gens of big group
        for g in big_gp_gens:
            ## Loop over all coset reps
            for key in self.coset_reps().keys():
                p = self.coset_reps()[key]
                ## compute p*g, represent as h * p_prime for h in subgroup
                product = p*g
                (h,p_prime) = self.represent_in_coset(product)
                h = self.level_1_gp(h)
                if not h.is_one():
                    ## check h is not 1
                    if not h in small_gens_matrices_dict: ## HACK: check we're not repeating gens (maybe better to kill this)
                        if not h^(-1) in small_gens_matrices_dict:
                        
                            ## This is new. Add h to the dictionary and add one to the index for next time
                            small_gens_matrices_dict[h] = current_index
                            current_index += 1
                    
                            ## also add h to the list of matrices
                            small_gens_matrices.append(h)

                            ##============================ [ REPLACE ORACLE ]================================
                            ## HACK! WE'RE ONLY DOING GAMMA0(1) in SL2(Z)
                            ## Now solve the word problem
                            word = self.level_1_word_problem.word_problem(h)
                            small_gens_words.append(word)
                    
        return small_gens_matrices_dict, small_gens_matrices, small_gens_words
        

    def level_N_word_problem(self,h):
        """
        Solve the word problem in the small group Gamma_0(N) in the list of generators output 
        by compute_generators.

        Firstly, we write this as h = 1.h. Then we write h = gh', where g in Gens(G) (so we must be
        able to solve the word problem for G). Then write 1.g = zp', so that
            h = z.p'h'. Now iterate. We will end up with z_1 z_2 ... z_t p_0, where p_0 = id rep.

        Outputs a list of ints in {-t,-t+1,...,t-1,t}, where the output of compute_generators is
        [a_1,...,a_t].

        EXAMPLE:
           h = abc in H, a,b,c in Gens(G)
           h = 1.abc
             = 1ap^(-1) . pbc
             = 1ap^(-1) . pbq^(-1) . qc
             = 1ap^(-1) . pbq^(-1) . qc1^(-1)
             = z1 . z2 . z3, where each zi is in the generating set.
        """
        ## Maybe write check to make sure h is in the right group?
        h = self.level_1_gp(h)
 
        ## compute the generators of H
        gens_dict, gens_matrices, gens_words = self.compute_generators()
        gens_G = self.level_1_gens ## gens of G

        ## Initialise final output
        h_level_N_wp = []

        ## Write h in the generators of g
        h_level_1_wp = self.level_1_word_problem.word_problem(h)
    
        ## We start with p_0 = id representative
        current_p = matrix([[1,0],[0,1]])

        ## loop through every generator that appears in the word of h (in G)
        for gen_ind in h_level_1_wp:

            ## Compute the generator we're processing
            current_gen = gens_G[ZZ(gen_ind).abs()-1]**ZZ(gen_ind).sign()
           
            ## Compute the generator and update p_i to p_{i+1}
            (h_current, current_p) = self.represent_in_coset(current_p * current_gen)
            h_current = self.level_1_gp(h_current)

            ## h_current should be one of our generators! As it is of form p'gp^(-1)
            if not h_current.is_one(): ## check not identity
                if not h_current in gens_dict:
                    ## h_current^(-1) should be in the dictionary
                    assert h_current^(-1) in gens_dict ## sanity check
                    ## Great, we've found a generator. Let's move on
                    ## Compute the index corresponding to this generator
                    gen_number = gens_dict[~h_current]
    
                    ## The generator appearing is an inverse. So append negative the index
                    h_level_N_wp.append(-gen_number-1)
              
                else:
                    assert h_current in gens_dict ## sanity check
                    ## Great, we've found a generator. Let's move on
                    ## Compute the index corresponding to this generator
                    gen_number = gens_dict[h_current]

                    ## The generator appearing is not inverse. So append the index
                    h_level_N_wp.append(gen_number+1)
        
        ## Check that we have actually solved the word problem correctly
        check_h = matrix([[1,0],[0,1]])
        for i in h_level_N_wp:
            check_h *= gens_matrices[ZZ(i).abs() - 1]**(ZZ(i).sign())
        assert check_h == h
            
        return h_level_N_wp


       
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


    def compute_cusp_stabiliser(self,cusp_matrix):
        """
        Compute (a finite index subgroup of) the stabiliser of a cusp 
        in Q or a quadratic imaginary field.
        
        We know the stabiliser of infinity is given by matrices of form 
        (u, a; 0, u^-1), so a finite index subgroup is generated by (1, alpha; 0, 1)
        and (1, 1; 0, 1) for K = Q(alpha). Given the cusp, we use a matrix
        sending infinty to that cusp, and the conjugate by it, before taking powers
        to ensure the result is integral and lies in Gamma_0(N).

        Input: 
            - a cusp (in matrix form: as output by cusp_set)
            - N (the level: an ideal in K).
        
        Outputs a list of the generators (as matrices).
        """

        P = self.get_P1List()
        if hasattr(P.N(),'number_field'):
            K = P.N().number_field()
            ## Write down generators of a finite index subgroup in Stab_Gamma(infinity)
            infinity_gens = [matrix(K,[[1,1],[0,1]]), matrix(K,[[1,K.gen()],[0,1]])]
            N_ideal = P.N()
        else:
            K = QQ
            infinity_gens = [matrix([[1,1],[0,1]])]
            N_ideal = ZZ.ideal(P.N())

        ## Initilise (empty) list of generators of Stab_Gamma(cusp)
        cusp_gens = []

        ## Loop over all the generators of stab at infinity, conjugate into stab at cusp
        for T in infinity_gens:
            T_conj = cusp_matrix * T * cusp_matrix^(-1)
            gen = T_conj

            ## Now take successive powers until the result is in Gamma_0(N)
            while not gen[1][0] in N_ideal:
                 gen *= T_conj

            ## We've found an element in Stab_Gamma(cusp): add to our list of generators
            cusp_gens.append(gen)
   
        return cusp_gens


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

    (The normalisation here might be wrong)
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
                ## Now have the matrix T = (u,h; 0,u^-1).
                ## Compute the action of this matrix on c.
                ## lift c = (c,d) to A' = (a,b;c,d), compute action of A'*T, normalise bottom row in P(O_F/N)
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


    

#for k in crt.keys():
#     print 'k = \n{}'.format(k)
#     print 'cusp = \n{}'.format(crt[k][0])
#     print 'matrix = \n{}'.format(crt[k][1])
#     print 'product c*m =\n{}'.format(crt[k][0]*crt[k][1])
#     print '======================'

    

X = var('X')
K.<a> = NumberField(X^2+1)
N = K.primes_above(13)[0]
#G = arith_gp(N)
G = arith_gp(5)
H = arith_gp(N)
h = G.level_1_gp(matrix([[6,1],[5,1]]))
