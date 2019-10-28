from sage.modular.pollack_stevens.padic_lseries import log_gamma_binomial

def basic_integral(self, a, j):
    ap = self._sign_ap
    p = self.parent().S_arithgroup().prime()
    g, a0, b = xgcd(a, -p)
    symb = self.evaluate(matrix(ZZ,2,2,[a,b,p,a0]))
    K = self.parent().coefficient_module().base_ring()

    return sum(ZZ(j).binomial(r)
                 * ((a - ZZ(K.teichmuller(a))) ** (j - r))
                 * (p ** r)
                 * (-1)**(r+1) * symb.moment(r) for r in range(j + 1)) / ap

def get_Lseries_term(self, n):
    r"""
    Return the `n`-th coefficient of the `p`-adic `L`-series

     """

    if n in self._Lseries_coefficients:
        return self._Lseries_coefficients[n]
    else:
        p = self.parent().S_arithgroup().prime()
        ap = self._sign_ap
        if p == 2:
            gamma = 1 + 4
        else:
            gamma = 1 + p
        K = self.parent().coefficient_module().base_ring()
        precision = K.precision_cap()

        S = QQ[['z']]
        z = S.gen()
        M = precision
        dn = 0
        if n == 0:
            precision = M
            lb = [1] + [0 for a in range(M - 1)]
        else:
            lb = log_gamma_binomial(p, gamma, z, n, 2 * M)
            if precision is None:
                precision = min([j + lb[j].valuation(p)
                                 for j in range(M, len(lb))])
            lb = [lb[a] for a in range(M)]

        for j in range(len(lb)):
            cjn = lb[j]
            temp = sum((ZZ(K.teichmuller(a)) ** (-j))
                       * self.basic_integral(a, j) for a in range(1, p))
            dn = dn + cjn * temp
        self._Lseries_coefficients[n] = dn.add_bigoh(precision)
        # self._Lseries_coefficients[n] /= self._cinf # DEBUG
        return self._Lseries_coefficients[n]
