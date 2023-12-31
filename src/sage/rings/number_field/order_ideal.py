r"""
Ideals of (Not Necessarily Maximal) Orders in Number Fields

This module implements (integral) ideals of orders in number fields.

.. NOTE::

    Currently, Sage only offers very limited functionality for
    ideals of non-maximal orders (compared to the maximal case).
    This should hopefully change in the future.

TESTS:

This module is currently experimental::

    sage: import sage.rings.number_field.order_ideal
    doctest:warning ...

EXAMPLES::

    sage: O = QuadraticField(-1).order(5*i)
    sage: I = O.ideal([13, 5*i-1]); I
    Ideal (60*a + 1, 65*a) of Order in Number Field in a with defining polynomial x^2 + 1 with a = 1*I

An ideal of an order in a relative number field::

    sage: K.<a,b> = NumberField([x^2 + 1, x^2 - 3])
    sage: O = K.order([3*a,2*b])
    sage: I = O.ideal((-6*b + 6)*a + 6*b + 18); I
    Ideal ((-60*b + 180)*a + 72, (-54*b + 174)*a - 6*b + 54, (-72*b + 288)*a + 72, 1872*a) of Relative Order in Number Field in a with defining polynomial x^2 + 1 over its base field

Perhaps the most useful functionality at this time is mapping ideals
of *quadratic* orders to corresponding binary quadratic forms::

    sage: K.<t> = QuadraticField(-21463)
    sage: O = K.order(t)
    sage: I = O.ideal([123457, t + 45259]); I
    Ideal (23058*t + 1, 123457*t) of Order in Number Field in t with defining polynomial x^2 + 21463 with t = 146.5025597046004?*I
    sage: I.quadratic_form()
    123457*x^2 - 90518*x*y + 16592*y^2

.. TODO::

    Generalize more functionality (such as primality testing and
    factoring) from
    :class:`~sage.rings.number_field.number_field_ideal.NumberFieldFractionalIdeal`
    to ideals of not necessarily maximal orders.

AUTHORS:

- Lorenz Panny (2022)
"""

# ****************************************************************************
#       Copyright (C) 2022 Lorenz Panny
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.free_module import Module_free_ambient
from sage.modules.free_module_element import vector
from sage.rings.ideal import Ideal_generic
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygens
from sage.rings.rational_field import QQ
from sage.structure.factorization import Factorization
from sage.structure.richcmp import richcmp
from sage.structure.sequence import Sequence

from sage.misc.superseded import experimental_warning
experimental_warning(34198, 'Ideals of non-maximal orders are an experimental feature. Be wary of bugs.')

import sage.rings.number_field.order

#TODO I*u works when u lies in I.ring().number_field(), but u*I doesn't

def NumberFieldOrderIdeal(O, *args, **kwds):
    r"""
    Construct either a :class:`NumberFieldOrderIdeal_generic`
    or a :class:`NumberFieldOrderIdeal_quadratic` from the
    given arguments.

    EXAMPLES::

        sage: from sage.rings.number_field.order_ideal import NumberFieldOrderIdeal
        sage: R.<x> = QQ[]
        sage: K.<t> = NumberField(x^3 - 40)
        sage: O = K.order(t)
        sage: I = NumberFieldOrderIdeal(O, [13, t-1]); I
        Ideal (12*t^2 + 1, 12*t^2 + t, 13*t^2) of Order in Number Field in t with defining polynomial x^3 - 40
        sage: K.absolute_degree()
        3
        sage: type(I)
        <class 'sage.rings.number_field.order_ideal.NumberFieldOrderIdeal_generic'>

    ::

        sage: L.<u> = QuadraticField(-3)
        sage: J = NumberFieldOrderIdeal(L.maximal_order(), [(u+5)/2])
        sage: L.absolute_degree()
        2
        sage: type(J)
        <class 'sage.rings.number_field.order_ideal.NumberFieldOrderIdeal_quadratic'>
    """
    if O.absolute_degree() == 2:
        cls = NumberFieldOrderIdeal_quadratic
    else:
        cls = NumberFieldOrderIdeal_generic
    return cls(O, *args, **kwds)

class NumberFieldOrderIdeal_generic(Ideal_generic):
    r"""
    An ideal of a not necessarily maximal order in a number field.
    """
    def __init__(self, O, gens, *, coerce=True):
        r"""
        Ideals of not necessarily maximal orders.

        The preferred way to construct objects of this class is via
        the :meth:`~sage.rings.number_field.order.Order.ideal` method.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<t> = NumberField(x^3 - 40)
            sage: O = K.order(t)
            sage: I = O.ideal([13, t-1]); I
            Ideal (12*t^2 + 1, 12*t^2 + t, 13*t^2) of Order in Number Field in t with defining polynomial x^3 - 40
            sage: I.norm()
            13

            sage: K.<a> = QuadraticField(-5)
            sage: O = K.order_of_conductor(37)
            sage: I = K.ideal([2, a + 1])
            sage: I.lift_order(O)
            Ideal (37*a + 1, 74*a) of Order in Number Field in a with defining polynomial x^2 + 5 with a = 2.236067977499790?*I
        """
        if not isinstance(O, sage.rings.number_field.order.Order):
            raise TypeError('not a number-field order')

        _, from_V, to_V = O.number_field().absolute_vector_space()
        if isinstance(gens, Module_free_ambient):
            # Intentional choice for code clarity
            pass

        else:
            if isinstance(gens, Ideal_generic):
                gens = gens.gens()

            if not isinstance(gens, (list, tuple)):
                gens = [gens]

            if coerce:
                gens = Sequence(gens, O)

            gens = [to_V(a*g) for a in O.basis() for g in gens]

        self._module = O.free_module().submodule(gens)
        basis = [O(from_V(v)) for v in self._module.basis()]

        super().__init__(O, basis, coerce=False)

    def __hash__(self):
        r"""
        Return a hash value for this ideal.

        EXAMPLES::

            sage: hash(QuadraticField(-1).order(5*i).ideal(7))  # random
            42
        """
        return hash((NumberFieldOrderIdeal_generic, self.ring(), self._module))

    def _richcmp_(self, other, op):
        r"""
        Implement comparison of ideals of number-field orders.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: g, = K.ring_of_integers().ring_generators()
            sage: O = K.order(567*g)
            sage: I = O.ideal([191, 567*t-27]); I
            Ideal (56133/2*t + 1/2, 108297*t) of Order in Number Field in t with defining polynomial x^2 + 123 with t = 11.09053650640942?*I
            sage: I == I
            True
            sage: I != I
            False
            sage: I == I^2
            False
        """
        if not isinstance(other, NumberFieldOrderIdeal_generic):
            return NotImplemented
        if self.ring() != other.ring():
            return NotImplemented
        return richcmp(self._module, other._module, op)

    def _contains_(self, x):
        r"""
        Test whether this ideal contains a given element.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: g, = K.ring_of_integers().ring_generators()
            sage: O = K.order(567*g)
            sage: I = O.ideal([191, 567*t-27]); I
            Ideal (56133/2*t + 1/2, 108297*t) of Order in Number Field in t with defining polynomial x^2 + 123 with t = 11.09053650640942?*I
            sage: 0 in I
            True
            sage: -191 in I
            True
            sage: 1 in I
            False
            sage: 567*t in I
            False
            sage: t-91 in I
            False
            sage: t-91 in K.maximal_order().ideal(I.gens(), future=1)
            True
        """
        x = self.ring()(x)
        return vector(QQ, x.list()) in self._module

    def __pari__(self):
        """
        Return PARI Hermite Normal Form representations of this ideal.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: K.<w> = NumberField(x^3 - 37)
            sage: I = K.maximal_order().ideal(K.class_group().0.ideal(), future=True); I
            Ideal (4/3*w^2 + 1/3*w + 1/3, w^2 + w, 2*w^2) of Maximal Order in Number Field in w with defining polynomial x^3 - 37
            sage: I.__pari__()
            [2, 1, 1; 0, 1, 0; 0, 0, 1]
        """
        return self.pari_hnf()

    def _pari_init_(self):
        """
        Return self in PARI Hermite Normal Form as a string.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: K.<w> = NumberField(x^3 - 37)
            sage: I = K.maximal_order().ideal(K.class_group().0.ideal(), future=True); I
            Ideal (4/3*w^2 + 1/3*w + 1/3, w^2 + w, 2*w^2) of Maximal Order in Number Field in w with defining polynomial x^3 - 37
            sage: I._pari_init_()
            '[2, 1, 1; 0, 1, 0; 0, 0, 1]'
        """
        return str(self.__pari__())

    def pari_hnf(self):
        """
        Return PARI's representation of this ideal in Hermite normal form.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: R.<x> = PolynomialRing(QQ)
            sage: K.<a> = NumberField(x^3 - 2)
            sage: I = K.ideal(2/(5+a))
            sage: I.pari_hnf()
            [2, 0, 50/127; 0, 2, 244/127; 0, 0, 2/127]
        """
        try:
            return self.__pari_hnf
        except AttributeError:
            nf = self.number_field().pari_nf()
            self.__pari_hnf = nf.idealhnf(0)
            hnflist = [nf.idealhnf(x) for x in self.gens()]
            for ideal in hnflist:
                self.__pari_hnf = nf.idealadd(self.__pari_hnf, ideal)
            return self.__pari_hnf

    def free_module(self):
        r"""
        Return the free `\ZZ`-module corresponding to this ideal
        as a submodule of the vector space associated to the ambient
        number field.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: g, = K.ring_of_integers().ring_generators()
            sage: O = K.order(567*g)
            sage: I = O.ideal([191, 567*t-27]); I
            Ideal (56133/2*t + 1/2, 108297*t) of Order in Number Field in t with defining polynomial x^2 + 123 with t = 11.09053650640942?*I
            sage: I.free_module()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [    1/2 56133/2]
            [      0  108297]
            sage: I.free_module().is_submodule(O.free_module())
            True

        .. SEEALSO::

            - :meth:`sage.rings.number_field.number_field.NumberField_absolute.absolute_vector_space`
            - :meth:`sage.rings.number_field.order.Order.free_module`
            - :meth:`sage.rings.number_field.number_field_ideal.NumberFieldIdeal.free_module`
        """
        return self._module

    def norm(self):
        r"""
        Return the norm of this ideal.

        The norm is defined as the index (as an abelian group)
        of the ideal in its order.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: g, = K.ring_of_integers().ring_generators()
            sage: O = K.order(567*g)
            sage: I = O.ideal([191, 567*t-27])
            sage: I.norm()
            191
            sage: (O.free_module() / I.free_module()).cardinality()
            191
        """
        return self.free_module().index_in(self.ring().free_module())

    def number_field(self):
        return self.ring().number_field()

    def factor(self):
        """
        Return factorization of this ideal in terms of prime ideals.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: K.<a> = NumberField(x^4 + 23); K
            Number Field in a with defining polynomial x^4 + 23
            sage: O_K = K.maximal_order(); O_K
            Maximal Order in Number Field in a with defining polynomial x^4 + 23
            sage: I = O_K.ideal(19, future=True); I
            Ideal (19/4*a^3 + 19/4*a^2 + 19/4*a + 19/4, 19/2*a^3 + 19/2*a, 19*a^2, 19*a^3) of Maximal Order in Number Field in a with defining polynomial x^4 + 23
            sage: F = I.factor(); F
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
        try:
            return self.__factorization
        except AttributeError:
            order = self.ring()
            K = self.number_field()
            F = K.pari_nf().idealfactor(self.pari_hnf())
            A = []
            for p, e in zip(F[0], F[1]):
                A.append((order.ideal(p, future=True), ZZ(e)))
            self.__factorization = Factorization(A)
            return self.__factorization

def _gens_from_bqf(O, Q):
    r"""
    Helper function to compute generators of the ideal associated to
    the given binary quadratic form.

    The resulting map induces a left-inverse of
    :meth:`NumberFieldOrderIdeal_quadratic.quadratic_form`
    on equivalence classes.

    REFERENCES:

    The correspondence itself is classical.
    Implemented after [Coh1993]_, §5.2.

    EXAMPLES::

        sage: K.<a> = QuadraticField(-18511)
        sage: O = K.order([83/2*a+1/2,83*a])
        sage: I = O.ideal([7868666196280491149/2*a + 107/2, 12934030292704639127*a])
        sage: F = I.quadratic_form(); F
        16673990859269835983*x^2 - 19847240181321592929*x*y + 5906098697962163510*y^2
        sage: J = O.ideal(F); J
        Ideal (7868666196280491149/2*a + 107/2, 12934030292704639127*a) of Order in Number Field in a with defining polynomial x^2 + 18511 with a = 136.0551358824797?*I
        sage: G = J.quadratic_form(); G
        16673990859269835983*x^2 - 19847240181321592929*x*y + 5906098697962163510*y^2
        sage: u = -479159/2*a + 665/2
        sage: J *= u; J
        Ideal (19274211614964208316630114412627/2*a + 43763/2, 33600036954602803764386186323012*a) of Order in Number Field in a with defining polynomial x^2 + 18511 with a = 136.0551358824797?*I
        sage: H = J.quadratic_form()
        sage: H.is_equivalent(F)
        True

    TESTS:

    Randomized testing::

        sage: from sage.rings.number_field.order_ideal import _random_for_testing
        sage: O, random_ideal = _random_for_testing()
        sage: I = random_ideal()
        sage: F = I.quadratic_form()
        sage: M = matrix.random_unimodular(MatrixSpace(ZZ,2))
        sage: assert M.det() == 1
        sage: G = F.matrix_action_left(M)
        sage: G.discriminant() == F.discriminant()
        True
        sage: J = O.ideal(G)
        sage: J.quadratic_form().is_equivalent(G)
        True
    """
    D = O.discriminant()
    if Q.discriminant() != O.discriminant():
        raise ValueError('order and form must have the same discriminant')
    a, b, c = Q
    sqrtD, = (s for s in O(D).sqrt(all=True) if s.real() > 0 or s.imag() > 0)
    g = (- b + sqrtD) / 2
    t = sqrtD if a < 0 else 1
    return a*t, g*t

class NumberFieldOrderIdeal_quadratic(NumberFieldOrderIdeal_generic):
    r"""
    An ideal of a not necessarily maximal order in a *quadratic* number field.
    """
    def __init__(self, O, gens, *, coerce=True):
        r"""
        Ideals of *quadratic* orders are implemented by a specialized
        class because they have some extra features not present in
        general number fields.

        The preferred way to construct objects of this class is via
        the :meth:`~sage.rings.number_field.order.Order.ideal` method.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-100)
            sage: O = K.order(t)
            sage: type(O.ideal(7))
            <class 'sage.rings.number_field.order_ideal.NumberFieldOrderIdeal_quadratic'>
        """
        from sage.quadratic_forms.binary_qf import BinaryQF
        if isinstance(gens, BinaryQF):
            gens = _gens_from_bqf(O, gens)
            coerce = False
        super().__init__(O, gens, coerce=coerce)

    def conjugate(self):
        r"""
        Return the conjugate of this ideal, defined by conjugating
        the generators.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-123)
            sage: g, = K.ring_of_integers().ring_generators()
            sage: O = K.order(567*g)
            sage: I = O.ideal([191, 567*t-27])
            sage: I.norm()
            191
            sage: I.norm() in I.conjugate() * I
            True
            sage: I.conjugate() * I == I.norm() * O
            True
        """
        conj_gens = [g.conjugate() for g in self.gens()]
        return NumberFieldOrderIdeal(self.ring(), conj_gens)

    def gens_two(self):
        r"""
        Express this ideal using exactly two generators, the first of
        which is a generator for the intersection of the ideal with `\ZZ`.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-100)
            sage: O = K.order(t)
            sage: I = O.ideal([123, 131-t, 21+23*t])
            sage: I.gens_two()
            (41, t - 8)
            sage: I == O.ideal(I.gens_two())
            True

        The second generator is zero if and only if the ideal is
        generated by an integer::

            sage: J = O.ideal([-33*t, 11*t-6589])
            sage: J.gens_two()
            (11, 0)
            sage: J == O.ideal(11)
            True

        .. WARNING::

            The returned generators do *not* necessarily
            form a `\ZZ`-basis of the ideal.

        TESTS:

        Random testing::

            sage: from sage.rings.number_field.order_ideal import _random_for_testing
            sage: O, random_ideal = _random_for_testing()
            sage: I = random_ideal()
            sage: gs2 = I.gens_two()
            sage: len(gs2)
            2
            sage: all(g in O for g in gs2)
            True
            sage: O.ideal(gs2) == I
            True
        """
        O = self.ring()
        if self.is_zero():
            return (O.zero(),)*2
        M = self._module & (ZZ**2).submodule([(1,0)])
        (N,_), = M.gens()
        NOgens = [N*g for g in O.free_module().basis()]
        Q = self._module / self._module.submodule(NOgens)
        if Q.invariants():
            assert len(Q.invariants()) == 1
            _, from_V, _ = O.number_field().absolute_vector_space()
            alpha, = (from_V(g.lift()) for g in Q.gens())
        else:
            alpha = 0
        return tuple(map(O, (N, alpha)))

    def is_principal(self):
        r"""
        Determine whether or not this ideal is principal.

        .. SEEALSO::

            To find a generator, use :meth:`gens_reduced`.

        EXAMPLES:

            sage: K.<a> = QuadraticField(-163)
            sage: O = K.order(7*a)
            sage: O.class_number()
            24
            sage: order = lambda v: next(e for e in range(1,99) if (v^e).is_principal())
            sage: I = O.ideal([47, 7*a-35])
            sage: order(I)
            24
            sage: J = O.ideal([71, 7*a-65])
            sage: order(J)
            12
            sage: next(e for e in range(99) if (I^e * J.conjugate()).is_principal())
            10
            sage: (I^10 * J.conjugate()).is_principal()
            True

        ::

            sage: K.<a> = QuadraticField(229)
            sage: O = K.order(7*a)
            sage: I = O.ideal([3, 7*a-2])
            sage: J = O.ideal([5, 7*a-4])
            sage: (I * J.conjugate()).is_principal()
            True
            sage: el = 104 + 7*a
            sage: el.norm()
            -405
            sage: (I * (J * el).conjugate()).is_principal()
            True
        """
        if self.is_zero():
            return True
        f = self.quadratic_form()
        sol = f.solve_integer(1)
        if sol is None:
            sol = f.solve_integer(-1)
        return sol is not None

    def gens_reduced(self):
        r"""
        Express this ideal in terms of at most two generators,
        and one if possible (i.e., if the ideal is principal).

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 + 11*x + 5)
            sage: O = K.order(7*a)
            sage: I = O.ideal([31915, -71145879*a - 32195694])
            sage: I.gens_reduced()
            (-63*a + 17,)

        ALGORITHM:

        Compute a reduction of the :meth:`quadratic_form` to see
        if it represents `1`, then use the transformation matrix
        to find an element in the ideal whose norm equals the
        norm of the ideal.

        TESTS::

            sage: from sage.rings.number_field.order_ideal import _random_for_testing
            sage: O, random_ideal = _random_for_testing()
            sage: I = random_ideal()
            sage: gs = I.gens_reduced()
            sage: len(gs) in (1,2)
            True
            sage: all(g in O for g in gs)
            True
            sage: (len(gs) == 1) == I.is_principal()
            True
            sage: O.ideal(gs) == I
            True
        """
        if self.is_zero():
            return (self.ring().zero(),)
        f, bas = self.quadratic_form(basis=True)
        sol = f.solve_integer(1)
        if sol is None:
            sol = f.solve_integer(-1)
        if sol is None:
            return self.gens_two()
        gen = sum(c*g for c,g in zip(sol, bas))
        assert NumberFieldOrderIdeal(self.ring(), gen) == self
        return (gen,)

    def is_equivalent(self, other, narrow=False):
        r"""
        Determine whether this ideal is equivalent to another ideal
        in the same order.

        If ``narrow`` is ``True``, test narrow equivalence instead.

        (Two ideals are equivalent if they differ by multiplication
        by a non-zero element. They are narrowly equivalent if they
        differ by multiplication by an element of positive norm.)

        EXAMPLES::

            sage: K.<a> = QuadraticField(-163)
            sage: O = K.order(7*a)
            sage: I = O.ideal([47, 7*a-35])
            sage: J = O.ideal([71, 7*a-65])
            sage: I.is_equivalent(J)
            False
            sage: (I^10).is_equivalent(J)
            True

        ::

            sage: K.<a> = QuadraticField(229)
            sage: O = K.order(7*a)
            sage: O.class_number()
            3
            sage: I = O.ideal([3, 7*a-2])
            sage: J = O.ideal([5, 7*a-4])
            sage: I.is_equivalent(J)
            True

        ::

            sage: K.<a> = QuadraticField(273)
            sage: O = K.order(11*a)
            sage: O.class_number()
            20
            sage: I = O.ideal([17, 11*a-11])
            sage: J = O.ideal([19, 11*a-12])
            sage: I.is_equivalent(J)
            False
            sage: (I^3).is_equivalent(J)
            False
            sage: (I^6).is_equivalent(J^2)
            True
            sage: el = 177 + 11*a
            sage: el.norm()
            -1704
            sage: (I^6).is_equivalent(J^2, narrow=True)
            True
            sage: (I^6).is_equivalent(J^2*el, narrow=True)
            False
        """
        if self.ring() != other.ring():
            raise ValueError('ideals must lie in the same order')
        if self.is_zero():
            return other.is_zero()
        if other.is_zero():
            return False
        gs = (self * other.conjugate()).gens_reduced()
        assert len(gs) in (1,2)
        if len(gs) > 1:
            return False
        elif narrow:
            return gs[0].norm() > 0
        return True

    def quadratic_form(self, *, basis=False):
        r"""
        Return the binary quadratic form associated to this ideal.

        This map induces an injective homomorphism from the narrow
        class group on ideals to the class group on quadratic forms.

        If ``basis`` is set to ``True`` (default: ``False``), the
        method additionally returns a `\ZZ`-basis `(a,b)` of this
        ideal `I` such that `f(x,y)` equals
        `\mathrm{norm}(xa+yb) / \mathrm{norm}(I)`,
        where `f` is the returned quadratic form.

        .. NOTE::

            The narrow class group is the group of invertible ideals
            modulo the principal ideals generated by an element of
            positive norm.

            - For *imaginary* quadratic orders, the narrow class group
              is identical to the class group.

            - For *real* quadratic orders, identifying the classes of
              `f(x,y)` and `-f(y,x)` recovers a correspondence with the
              standard class group.

        REFERENCES:

        The correspondence itself is classical.
        Implemented after [Coh1993]_, §5.2.

        EXAMPLES::

            sage: K.<t> = QuadraticField(-419)
            sage: O = K.order(t)
            sage: O.discriminant().factor()
            -1 * 2^2 * 419
            sage: I = O.ideal([t-1, 105]); I
            Ideal (104*t + 1, 105*t) of Order in Number Field in t with defining polynomial x^2 + 419 with t = 20.46948949045873?*I
            sage: f = I.quadratic_form(); f
            105*x^2 - 208*x*y + 107*y^2
            sage: f.discriminant().factor()
            -1 * 2^2 * 419
            sage: power(f,3).reduced_form()
            x^2 + 419*y^2

        ::

            sage: u = 23*t - 45
            sage: J = I*u
            sage: g = J.quadratic_form(); g
            23485980*x^2 - 22795498*x*y + 5531329*y^2
            sage: f.is_equivalent(g)
            True

        The inverse operation (modulo equivalence) can be computed by
        passing a :class:`~sage.quadratic_forms.binary_qf.BinaryQF`
        to ``O.ideal()``::

            sage: II = O.ideal(f); II
            Ideal (104*t + 1, 105*t) of Order in Number Field in t with defining polynomial x^2 + 419 with t = 20.46948949045873?*I
            sage: II.quadratic_form().is_equivalent(f)
            True

        TESTS:

        Randomized testing (principal ideals generated by an element
        of positive norm map to the principal form)::

            sage: from sage.rings.number_field.order_ideal import _random_for_testing
            sage: O, _ = _random_for_testing()
            sage: D = O.discriminant()
            sage: gen = (O.random_element() for _ in iter(int,1))
            sage: g = next(el for el in gen if el.norm() > 0)
            sage: I = O.ideal(g)
            sage: F = I.quadratic_form()
            sage: F.is_equivalent(BinaryQF.principal(D))
            True

        Randomized testing (mapping is a homomorphism)::

            sage: from sage.rings.number_field.order_ideal import _random_for_testing
            sage: O, random_ideal = _random_for_testing()
            sage: I = random_ideal()
            sage: J = random_ideal()
            sage: F = I.quadratic_form()
            sage: G = J.quadratic_form()
            sage: H = (I*J).quadratic_form()
            sage: F.discriminant() == G.discriminant() == H.discriminant() == O.discriminant()
            True
            sage: H.is_equivalent(F*G)
            True

        Constructing an ideal from a form is indeed a one-sided inverse::

            sage: II = O.ideal(F)
            sage: II.quadratic_form().is_equivalent(F)
            True
        """
        if self.is_zero():
            if basis:
                return BinaryQF(0), (self.ring().zero(),)*2
            return BinaryQF(0)

        # find a ZZ-basis of the ideal with a in QQ
        M = self._module.basis_matrix()[:,::-1]
        M = M.row_space(ZZ).basis_matrix()[:,::-1]
        b,a = map(self.ring(), M.rows())
        if b.real() < 0 or b.imag() < 0:
            b = -b

        # compute the (reduced) norm form of the ideal
        A,B = (g.matrix() for g in (a,b))
        x,y = polygens(QQ, 'x,y')
        Q = (x*A - y*B).determinant() / self.norm()
        Q = Q.change_ring(ZZ)

        from sage.quadratic_forms.binary_qf import BinaryQF
        if basis:
            return BinaryQF(Q), (a,-b)
        return BinaryQF(Q)


def _random_for_testing():
    r"""
    A small helper function to produce somewhat random examples
    of ideals of (not necessarily maximal) quadratic orders.

    Used in doctests.

    EXAMPLES::

        sage: from sage.rings.number_field.order_ideal import _random_for_testing
        sage: O, random_ideal = _random_for_testing()
        sage: random_ideal().ring() is O
        True
    """
    from sage.misc.prandom import choice
    from sage.misc.prandom import randrange
    from sage.rings.number_field.number_field import QuadraticField
    from sage.arith.misc import primes
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as Zmod

    while True:
        d = ZZ(choice((-1,+1)) * randrange(1,10**5))
        if not d.is_square():
            break

    K = QuadraticField(d)
    g, = K.ring_of_integers().ring_generators()

    f = randrange(1, 100)
    O = K.order(f*g)
    assert O.discriminant() == f**2 * K.discriminant()

    poly = (f*g).minpoly()
    base = []
    for l in primes(1000):
        vs = poly.roots(ring=Zmod(l), multiplicities=False)
        base += [NumberFieldOrderIdeal(O, [l, f*g-ZZ(v)]) for v in vs]
    def random_ideal():
        I = NumberFieldOrderIdeal(O, [1])
        for _ in range(randrange(20)):
            J = choice(base)
            I = NumberFieldOrderIdeal(O, [x*y for x in I.gens() for y in J.gens()])
        return I
    return O, random_ideal
