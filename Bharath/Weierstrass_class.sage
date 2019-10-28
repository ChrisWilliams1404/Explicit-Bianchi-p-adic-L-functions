r'''
TESTS

sage: p=5
sage: Zpp = Zp(5, prec = 36, type = 'fixed-mod', print_mode = 'series', show_prec=False)
sage: R.<x>  = PowerSeriesRing(Zpp, default_prec=41);
sage: S.<y> = PowerSeriesRing(R, default_prec=41); 
sage: U = 1/(1+x^2 + y^3 + p^2)
sage: poly  = x^2 + p + p^3*y + y^2
sage: F = U*poly
sage: print(Weierstrass(F))
5 + x^2 + 5^3*y + y^2
'''

# The following is the required precision for x, y and p 

#x_prec = 10
#y_prec = 10
#p_prec = 5
#max_lamb_prec = 6


# define the rings we will have to work with: they will need to have higher precision


class WeierstrassPrep(object):
    
    def __init__(self, f, xprec=10, yprec=10, pprec=5, maxlambprec=6):
        """
        Initialisation. Input two variable power series f and some bunch of precisions.
        
        # x_prec: precision for the polynomial in the variable x
        # y_prec: precision for the polynomial in the variable y
        # p_prec: precision for the power of p
        # max_lamb_prec : Estimate the maximal lambda invariant. 


        # To set the auxillary precision, the auxillary increment must be the sum of precisions required for p,x and y. 
        # This allows us to ignore terms of the form p^r*x^m*y^n everywhere, where r + m + n >= x_prec + y_prec + z_prec
    
        """
        
        self.f = f
        self.S = f.parent()
        self.R = self.S.base_ring()
        self.Zpp = self.R.base_ring()
        
        self.p = self.Zpp.prime()
    
        # x and y are the variables
        self.x = self.R.gen()
        self.y = self.S.gen()
        
        # precision
        self.aux_prec_inc = maxlambprec + xprec + yprec + pprec
    
    
        ## Low precision versions
        self.Zpplow = self.Zpp.change(prec=pprec);      
        self.Rlow  = PowerSeriesRing(self.Zpplow, default_prec=xprec, names = 'x');
        self.xlow = self.Rlow.gen()
        self.Slow = PowerSeriesRing(self.Rlow, default_prec=yprec, names = 'y');                              # S: Power Series ring over R = Z_p[[x]]
        self.ylow = self.Slow.gen()
        
        self.xprec = xprec
        self.yprec = yprec
        self.pprec = pprec
        self.maxlambprec = maxlambprec


    def __repr__(self):
            """
            Represent when called.
            """
            return 'Weierstrass preparation class: runs Weierstrass preparation on a chosen power series.'
            


    def Weierstrass(self):
            """
            Main function.
            """
    
            return self.weierstrass_prep()
    



    def shift_by_n(self,g,n):
            # Input polynomial: g
            # Input shift: n
        
            deg_g = g.degree(); 
        
        
            g_coeff  = g.list();
        
            h = self.S([0]); 
        
            for i in range(0,deg_g+1):
                if (i+n > deg_g):
                    break; 
                h = h + g_coeff[i+n]*self.y^(i);    
                
            return h;
    
    
    
    
    
        #Computing Lambda invariant
    
    def lamb_invariant(self):
            lamb = -infinity;
    
            f_coeff = self.f.list();
        
            for i in range(self.f.degree()+1):
                if (f_coeff[i].is_unit()):
                    lamb = i;
                    break;
            
            return lamb;
        
        
        
    #Weierstrass preparation over Z_p[[x]][y]]. This can be modified to work over any ring 
    
    def weierstrass_prep(self):
          
            weier_prec =  self.xprec + self.yprec + self.pprec + self.maxlambprec+ 4*(self.xprec+self.yprec+self.pprec+self.maxlambprec) ;                  
        
            f_i = self.f;
        
            lambda_f = self.lamb_invariant();
            if (lambda_f == -infinity):
                print "This is not a Weierstrass power series! Boo"
                return;
            
            for i in range(0,weier_prec+1):
                f_i = f_i/(self.shift_by_n(f_i,lambda_f));
        
            # If you want to work with rings other than Z_p[[x]][[y]], then at this stage, you can output f_i.truncate(lambda_f+1)
    
        
            # First, let's truncate the terms required in the variable y
            ff_trun_y =  f_i.truncate(lambda_f+1);
        
            # Now, let's truncate each of the coefficients of ff_trun_y in the variable x
            ff_trun_y_x = self.Slow([0]); 
            for i in range(0,lambda_f+1):
                ff_trun_y_x =  ff_trun_y_x + (self.Rlow(ff_trun_y.list()[i])).truncate(self.xprec)*self.ylow^i;
        
        
            #R.set_default_prec(xprec);
            #S.set_default_prec(yprec); 
            # Change the precision again
        
        
            # ff_trun_y_x is the required polynomial with the required precision.
            return self.S(ff_trun_y_x)
            
            
## EXAMPLE CODE


#p=5
#Zpp = Zp(5, prec = 36, type = 'fixed-mod', print_mode = 'series', show_prec=False)
#R.<T>  = PowerSeriesRing(Zpp, default_prec=41);
#S.<Z> = PowerSeriesRing(R, default_prec=41); 
#U = 1/(1+T^2 + Z^3 + p^2)
#poly  = T^2 + p^2*T + p^3*Z + Z^6
#F = U*poly
#W = WeierstrassPrep(F)
#print W.Weierstrass()
