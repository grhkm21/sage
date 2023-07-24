def torsion_sample(Et, l):
    # assert Et.order() % l**2 == 0
    S = Et.order()
    s = S // l**valuation(S, l)
    while (P := Et.random_point() * s) == Et(0):
        pass
    while P * l != Et(0):
        P *= l

    return P


# Computes (Psi_l(x, j) mod p, roots j
def algo1(j, l, Kq, Kp2x):
    # Construct objects and parameters
    Et = EllipticCurve(Kq, j=j)
    assert Et.is_supersingular()

    # TODO: Replace with P, Q = Et.torsion_basis(l)
    # When that method is fixed, as it is currently too slow
    # P, Q = Et.torsion_basis(l)

    P = Q = torsion_sample(Et, l)
    while P.weil_pairing(Q, l) == 1:
        Q = torsion_sample(Et, l)

    gens = [Q] + [P + (i - 1) * Q for i in range(1, l + 1)]
    jinvs = [Et.isogeny_codomain(G).j_invariant() for G in gens]

    t = Kq['t'].gen()
    mpoly = prod(t - j for j in jinvs)
    return mpoly, jinvs


def number_of_supersingular(p):
    eps = [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2]
    return p // 12 + eps[p % 12]


# Computes Psi_l(x, y) mod p using supersingular curves
def modular_polynomial_mod_p(l, p, names="x, y"):
    if number_of_supersingular(p) <= l:
        raise NotImplementedError(f"{number_of_supersingular(p) = } <= {l = }")

    from sage.structure.category_object import normalize_names
    vx, vy = normalize_names(2, names)

    Kq = FiniteField(p**(6 * (l - 1)), names="t")
    Kp2 = FiniteField(p**2, names="u")
    Kp2x = Kp2[vx]

    # Embedding from Kq -> Kp2
    # It is defined as the "inverse" (categorically, section) of Kp2 -> Kq
    # There are two such maps but it doesn't matter for us
    # This is not canonical since we specified the `names` argument above and
    # hence the modulus of Kq and Kp2 is randomised
    phi = list(Hom(Kp2, Kq))[0]
    iphi = phi.section()

    if p % 4 == 3:
        j0 = Kp2(1728)
    else:
        d = 1
        while legendre_symbol(d, p) != -1:
            d += 1
        # TODO: Replace hilbert_class_polynomial with mod p version
        j0 = hilbert_class_polynomial(-4 * d).any_root(Kp2)

    # TODO: Benchmark with Lagrange interpolation
    seen = set()
    roots = set([j0])
    coefs = Matrix(Kp2x, l + 1, l + 1)
    rhs = vector(Kp2x, l + 1)
    for i in range(l + 1):
        _roots = list(roots)
        while (j := _roots.pop()) in seen:
            pass
        roots = set(_roots)

        mpoly, jinvs = algo1(phi(j), l, Kq, Kp2x)
        for t in range(l + 1):
            coefs[i, t] = j**t
        # Compute polynomial in Kq[x] then cast back to Kp2[x]
        mpoly = Kp2x(list(map(iphi, mpoly)))
        rhs[i] = mpoly - j**(l + 1)
        seen.add(j)

        for u in map(iphi, set(jinvs)):
            if u not in seen:
                roots.add(u)

    assert coefs.right_nullity() == 0
    coef = coefs.solve_right(rhs)

    Kp2xy = Kp2x[vy]
    return Kp2xy.flattening_morphism()(Kp2xy(coef.list() + [1]))


"""
Benchmark:

23 Jul 22:30
================
sage: time modular_polynomial_mod_p(5, 307)
CPU times: user 530 ms, sys: 6.32 ms, total: 537 ms

sage: time modular_polynomial_mod_p(7, 307)
CPU times: user 1.09 s, sys: 12.9 ms, total: 1.11 s

sage: time modular_polynomial_mod_p(11, 307)
CPU times: user 4.69 s, sys: 62.1 ms, total: 4.75 s

sage: time modular_polynomial_mod_p(13, 307)
CPU times: user 7.76 s, sys: 61.5 ms, total: 7.82 s

sage: time modular_polynomial_mod_p(17, 307)
CPU times: user 21.7 s, sys: 298 ms, total: 22 s

sage: time modular_polynomial_mod_p(23, 307)
CPU times: user 7min 21s, sys: 6.91 s, total: 7min 28s

24 Jul 02:39

Changelog:
- The finite fields GF(p^l) are replaced with GF(p^l, names='t') that
    uses a random irreducible modulus instead of Conway polynomials

sage: time modular_polynomial_mod_p(11, 149)
CPU times: user 3.5 s, sys: 35.1 ms, total: 3.53 s

sage: time modular_polynomial_mod_p(11, 307)
CPU times: user 3.79 s, sys: 43.3 ms, total: 3.84 s

sage: time modular_polynomial_mod_p(13, 307)
CPU times: user 6.53 s, sys: 53.9 ms, total: 6.59 s

sage: time modular_polynomial_mod_p(17, 307)
CPU times: user 16.4 s, sys: 157 ms, total: 16.5 s

sage: time modular_polynomial_mod_p(23, 307)
CPU times: user 42.1 s, sys: 523 ms, total: 42.6 s
"""

