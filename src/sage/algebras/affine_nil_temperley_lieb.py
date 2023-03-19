"""
Affine nilTemperley Lieb Algebra of type A
"""
# ****************************************************************************
#  Copyright (C) 2010 Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.rings.ring import Ring
from sage.rings.integer_ring import ZZ
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method


class AffineNilTemperleyLiebTypeA(CombinatorialFreeModule):
    r"""
    Construct the affine nilTemperley Lieb algebra of type `A_{n-1}^{(1)}` as used in [Pos2005]_.

    INPUT:

     - ``n`` -- a positive integer

    The affine nilTemperley Lieb algebra is generated by `a_i` for `i=0,1,\ldots,n-1`
    subject to the relations `a_i a_i = a_i a_{i+1} a_i = a_{i+1} a_i a_{i+1} = 0` and
    `a_i a_j = a_j a_i` for `i-j \not \equiv \pm 1`, where the indices are taken modulo `n`.

    EXAMPLES::

        sage: A = AffineNilTemperleyLiebTypeA(4)
        sage: a = A.algebra_generators(); a
        Finite family {0: a0, 1: a1, 2: a2, 3: a3}
        sage: a[1]*a[2]*a[0] == a[1]*a[0]*a[2]
        True
        sage: a[0]*a[3]*a[0]
        0
        sage: A.an_element()
        2*a0 + 1 + 3*a1 + a0*a1*a2*a3
    """

    def __init__(self, n, R=ZZ, prefix='a'):
        """
        Initiate the affine nilTemperley Lieb algebra over the ring `R`.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3, prefix="a"); A
            The affine nilTemperley Lieb algebra A3 over the ring Integer Ring
            sage: TestSuite(A).run()
            sage: A = AffineNilTemperleyLiebTypeA(3, QQ); A
            The affine nilTemperley Lieb algebra A3 over the ring Rational Field
        """
        if not isinstance(R, Ring):
            raise TypeError("argument R must be a ring")
        self._cartan_type = CartanType(['A', n - 1, 1])
        self._n = n
        W = WeylGroup(self._cartan_type)
        self._prefix = prefix
        self._index_set = W.index_set()
        self._base_ring = R
        category = AlgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(self, R, W, category=category)

    def _element_constructor_(self, w):
        """
        Constructs a basis element from an element of the Weyl group.

        If `w = w_1 ... w_k` is a reduced word for `w`, then `A(w)` returns
        zero if `w` contains braid relations.
        TODO: Once the functorial construction is in sage, perhaps this should be
        handled constructing the affine nilTemperley Lieb algebra as a quotient algebra.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3, prefix="a")
            sage: W = A.weyl_group()
            sage: w = W.from_reduced_word([2,1,2])
            sage: A(w)
            0
            sage: w = W.from_reduced_word([2,1])
            sage: A(w)
            a2*a1
        """
        W = self.weyl_group()
        assert w in W
        word = w.reduced_word()
        if all(self.has_no_braid_relation(W.from_reduced_word(word[:i]), word[i]) for i in range(len(word))):
            return self.monomial(w)
        return self.zero()

    @cached_method
    def one_basis(self):
        """
        Return the unit of the underlying Weyl group, which index
        the one of this algebra, as per
        :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: A.one_basis()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A.one_basis() == A.weyl_group().one()
            True
            sage: A.one()
            1
        """
        return self.weyl_group().one()

    def _repr_(self):
        """
        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3); A
            The affine nilTemperley Lieb algebra A3 over the ring Integer Ring
        """
        return "The affine nilTemperley Lieb algebra A%s over the ring %s"%(self._n, self._base_ring)

    def weyl_group(self):
        """
        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: A.weyl_group()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root space)
        """
        return self.basis().keys()

    def index_set(self):
        """
        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: A.index_set()
            (0, 1, 2)
        """
        return self._index_set

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators `a_i` for `i=0,1,2,\ldots,n-1`.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: a = A.algebra_generators();a
            Finite family {0: a0, 1: a1, 2: a2}
            sage: a[1]
            a1
        """
        return self.weyl_group().simple_reflections().map(self.monomial)

    def algebra_generator(self, i):
        """
        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: A.algebra_generator(1)
            a1
            sage: A = AffineNilTemperleyLiebTypeA(3, prefix = 't')
            sage: A.algebra_generator(1)
            t1
        """
        return self.algebra_generators()[i]

    def product_on_basis(self, w, w1):
        """
        Return `a_w a_{w1}`, where `w` and `w1` are in the Weyl group
        assuming that `w` does not contain any braid relations.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(5)
            sage: W = A.weyl_group()
            sage: s = W.simple_reflections()
            sage: [A.product_on_basis(s[1],x) for x in s]
            [a1*a0, 0, a1*a2, a3*a1, a4*a1]

            sage: a = A.algebra_generators()
            sage: x = a[1] * a[2]
            sage: x
            a1*a2
            sage: x * a[1]
            0
            sage: x * a[2]
            0
            sage: x * a[0]
            a1*a2*a0

            sage: [x * a[1] for x in a]
            [a0*a1, 0, a2*a1, a3*a1, a4*a1]

            sage: w = s[1]*s[2]*s[1]
            sage: A.product_on_basis(w,s[1])
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert self(w) != self.zero()
        for i in w1.reduced_word():
            if self.has_no_braid_relation(w, i):
                w = w.apply_simple_reflection(i)
            else:
                return self.zero()
        return self.monomial(w)

    @cached_method
    def has_no_braid_relation(self, w, i):
        """
        Assuming that `w` contains no relations of the form `s_i^2` or `s_i s_{i+1} s_i` or
        `s_i s_{i-1} s_i`, tests whether `w s_i` contains terms of this form.

        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(5)
            sage: W = A.weyl_group()
            sage: s=W.simple_reflections()
            sage: A.has_no_braid_relation(s[2]*s[1]*s[0]*s[4]*s[3],0)
            False
            sage: A.has_no_braid_relation(s[2]*s[1]*s[0]*s[4]*s[3],2)
            True
            sage: A.has_no_braid_relation(s[4],2)
            True
        """
        if w == w.parent().one():
            return True
        if i in w.descents():
            return False
        s = w.parent().simple_reflections()
        wi = w*s[i]
        adjacent = [(i-1)%w.parent().n, (i+1)%w.parent().n]
        for j in adjacent:
            if j in w.descents():
                if j in wi.descents():
                    return False
                else:
                    return True
        return self.has_no_braid_relation(w*s[w.first_descent()],i)

    def _repr_term(self, t, short_display=True):
        """
        EXAMPLES::

            sage: A = AffineNilTemperleyLiebTypeA(3)
            sage: W = A.weyl_group()
            sage: A._repr_term(W.from_reduced_word([1,2,0]))
            'a1*a2*a0'
            sage: A._repr_term(W.from_reduced_word([1,2,0]), short_display = False)
            'a[1]*a[2]*a[0]'
        """
        redword = t.reduced_word()
        if len(redword) == 0:
            return "1"
        elif short_display:
            return "*".join("%s%d"%(self._prefix, i) for i in redword)
        else:
            return "*".join("%s[%d]"%(self._prefix, i) for i in redword)
