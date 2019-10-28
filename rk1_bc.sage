#F: Imaginary quadratic field
#D: F.ideal(1)
#M: Ideal Conductor of E (over F) in F 


def getKey(item):
    return item[6]


def fwrite(string, outfile,erase_file=False,newline = True):
    if erase_file:
        code = "w"
    else:
        code ="a"
        
    if outfile is None:
        fout = sys.stdout
        if newline:
            fout.write(string + '\n')
        else:
            fout.write(string)
    else:
        with open(outfile,code) as fout:
            if newline:
                fout.write(string + '\n')
            else:
                fout.write(string)
        return


def covolume(F,M,D=1,prec = None,zeta = None):
    
    from sage.symbolic.constants import pi
    n = F.degree()
    if prec is None:
        prec = 53
    disc = ZZ(F.discriminant())
    if n > 1:
        if zeta is None:
            zetaf = RealField(prec)(magma_free('Evaluate(LSeries(QuadraticField(%s)),2)'%disc))
        else:
            zetaf = zeta
    else:
        from sage.functions.transcendental import Function_zeta
        if zeta is None:
            zetaf = Function_zeta()(RealField(prec)(2))
        else:
            zetaf = zeta

    if F != QQ:
        D = F.ideal(D)
        M = F.ideal(M)
        Phi = QQ(1)
        for P,_ in D.factor():
            Phi *= QQ(P.norm().abs() - 1)
        Psi = QQ(M.norm()).abs()
        for P,e in M.factor():
            np = QQ(P.norm())
            Psi *= np**(ZZ(e)-1) * (np + 1)
    else:
        M = ZZ(M)
        Phi = ZZ(D)
        for np,_ in D.factor():
            Phi *= QQ(1)-QQ(1)/np
        Psi = ZZ(M).abs()
        for np,e in M.factor():
            Psi *= np**(ZZ(e)-1) * (np + 1)
    RR = RealField(prec)
    pi = RR(pi)
    covol =  (RR(disc).abs()**(3.0/2.0) * zetaf * Phi)/((4 * pi**2)**(F.degree()-1))
    index = RR(Psi)
    indexunits = 1 # There is a factor missing here, due to units.
    return covol * index / indexunits
    
    
    
    
class Rank_List(object):
    
    from array import array
    
    def __init__(self, rank, upper_bound=100):
        """
        Stores
    
        """
        self.rank = rank;                               #You want to produce elliptic curves, that are base-change over Q of this rank.
        self.cremona = CremonaDatabase();               # Loads the Cremona database
        self.upper_bound = upper_bound                  #You want to search for elliptic curves over Q up to this conductor as the upper bound
        self.Disc_List = [-3,-1,-7,-2,-11]
        self.list = []


    def __repr__(self):
            
            """
            Represent when called.
            """
            
            return 'List of Elliptic curves whose base changes are the given rank'
            
    def generate_ell_curve_1 (self):
            """
            Main function.
            """
        
                      
            
            ind = 0;
            
            for Disc in self.Disc_List:
            
                
                 
            
                K.<gg> = QuadraticField(Disc);                 # Initialize the Number Field
                
                while True:
                    ans = magma_free('Evaluate(LSeries(QuadraticField(%s)),2)'%Disc)
                    if ans != '':
                        K_zetaf = RealField(20)(ans)
                        break
                    sleep(60)                    
                
                
                for cond in range(11,self.upper_bound):
                    for curve_label in self.cremona.curves(cond):
                        
                        E = self.cremona.elliptic_curve(str(cond)+str(curve_label));
                        rank_E = E.rank()
                    
                        if (E.rank()==1):
                        
                            E_twist = E.quadratic_twist(Disc);
                            rank_E_twist = E_twist.rank()
                        
                            if (rank_E_twist == 0):
                            
                                
                                temp_list= [];
                                
                                
                                temp_list.append(E);                               # Elliptic curve over Q
                                
                                temp_list.append(K);                               # The discriminant
                                temp_list.append(K.defining_polynomial())         # Defining polynomial for the imaginary quadratic field
                                temp_list.append([rank_E,rank_E_twist])
                                E_cond = E.conductor()
                                temp_list.append(E_cond)
                                
                                for p in primes (3,100):
                                
                                    if (len(K.ideal(p).factor())==1):
                                        continue;
                                    
                                    if (Disc % p ==0 ):
                                        continue;
                                
                                    if (E_cond  % p == 0):
                                        continue;
                                
                                    Ered = E.reduction(p)
                                
                                    if (Ered.is_ordinary()):
                                        continue;
                                
                                    if (E.ap(p) !=0 ):
                                        continue;
                                
                                    EK = E.base_extend(K);
                                    base_change_conductor = EK.conductor()
                                    
                                    p_temp_list = []
                                    for ii in temp_list: p_temp_list.append(ii);
                                    
                                    p_temp_list.append(p)
                                    p_temp_list.append(covolume(K,p*base_change_conductor,zeta=K_zetaf));
                                    self.list.append(p_temp_list);

                        if (E.rank()==0):
                        
                            E_twist = E.quadratic_twist(Disc);
                            rank_E_twist = E_twist.rank()
                        
                            if (rank_E_twist == 1):
                            
                                 
                                temp_list= [];
                                
                                
                                temp_list.append(E);                               # Elliptic curve over Q
                               
                                temp_list.append(K);                               # The discriminant
                                temp_list.append(K.defining_polynomial())         # Defining polynomial for the imaginary quadratic field
                                temp_list.append([rank_E,rank_E_twist])
                                E_cond = E.conductor()
                                temp_list.append(E_cond)
                                
                                for p in primes (3,100):
                                
                                    if (len(K.ideal(p).factor())==1):
                                        continue;
                                    
                                    if (Disc % p ==0 ):
                                        continue;
                                
                                    if (E_cond  % p == 0):
                                        continue;
                                
                                    Ered = E.reduction(p)
                                
                                    if (Ered.is_ordinary()):
                                        continue;
                                
                                    if (E.ap(p) !=0 ):
                                        continue;
                                
                                    EK = E.base_extend(K);
                                    base_change_conductor = EK.conductor()
                                
                                    p_temp_list = []
                                    for ii in temp_list: p_temp_list.append(ii); 
                                    
                                    p_temp_list.append(p)
                                    p_temp_list.append(covolume(K,p*base_change_conductor,zeta=K_zetaf));
                                    self.list.append(p_temp_list);
                                                                   
            self.list = sorted(self.list, key =getKey)
            return;
                        
## Example code


Y = Rank_List(1)
Y.generate_ell_curve_1()  
print "\n\n The first curve in the list with its attributes:", Y.list[0]
print "\n\n The second curve in the list with its attributes: ", Y.list[1]
#print Y.list[2]
#print Y.list[3]                        
                            
                        
                    
        
            
