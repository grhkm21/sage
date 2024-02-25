r"""
Quotients of modules with basis
"""
# ****************************************************************************
#  Copyright (C) 2010-2015 Florent Hivert <Florent.Hivert@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.modules_with_basis import ModulesWithBasis


class QuotientModuleWithBasis(CombinatorialFreeModule):
    r"""
    A class for quotients of a module with basis by a submodule.

    INPUT:

    - ``submodule`` -- a submodule of ``self``
    - ``category`` -- a category (default: ``ModulesWithBasis(submodule.base_ring())``)

    ``submodule`` should be a free submodule admitting a basis in
    unitriangular echelon form. Typically ``submodule`` is a
    :class:`SubmoduleWithBasis` as returned by
    :meth:`Modules.WithBasis.ParentMethods.submodule`.

    The ``lift`` method should have a method
    ``.cokernel_basis_indices`` that computes the indexing set of a
    subset `B` of the basis of ``self`` that spans some supplementary
    of ``submodule`` in ``self`` (typically the non characteristic
    columns of the aforementioned echelon form). ``submodule`` should
    further implement a ``submodule.reduce(x)`` method that returns
    the unique element in the span of `B` which is equivalent to `x`
    modulo ``submodule``.

    This is meant to be constructed via
    :meth:`Modules.WithBasis.FiniteDimensional.ParentMethods.quotient_module`

    This differs from :class:`sage.rings.quotient_ring.QuotientRing`
    in the following ways:

    - ``submodule`` needs not be an ideal. If it is, the
      transportation of the ring structure is taken care of by the
      ``Subquotients`` categories.

    - Thanks to ``.cokernel_basis_indices``, we know the indices of a
      basis of the quotient, and elements are represented directly in
      the free module spanned by those indices rather than by wrapping
      elements of the ambient space.

    There is room for sharing more code between those two
    implementations and generalizing them. See :trac:`18204`.

    .. SEEALSO::

        - :meth:`Modules.WithBasis.ParentMethods.submodule`
        - :meth:`Modules.WithBasis.FiniteDimensional.ParentMethods.quotient_module`
        - :class:`SubmoduleWithBasis`
        - :class:`sage.rings.quotient_ring.QuotientRing`
    """
    @staticmethod
    def __classcall_private__(cls, submodule, category=None):
        r"""
        Normalize the input.

        TESTS::

            sage: from sage.modules.with_basis.subquotient import QuotientModuleWithBasis
            sage: X = CombinatorialFreeModule(QQ, range(3)); x = X.basis()
            sage: I = X.submodule( (x[0]-x[1], x[1]-x[2]) )
            sage: J1 = QuotientModuleWithBasis(I)
            sage: J2 = QuotientModuleWithBasis(I, category=Modules(QQ).WithBasis().Quotients())
            sage: J1 is J2
            True
        """
        default_category = ModulesWithBasis(submodule.category().base_ring()).Quotients()
        category = default_category.or_subcategory(category, join=True)
        return super().__classcall__(cls, submodule, category)

    def __init__(self, submodule, category):
        r"""
        Initialize this quotient of a module with basis by a submodule.

        TESTS::

            sage: from sage.modules.with_basis.subquotient import QuotientModuleWithBasis
            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: I = X.submodule( (x[0]-x[1], x[1]-x[2]) )
            sage: Y = QuotientModuleWithBasis(I)
            sage: Y.print_options(prefix='y')
            sage: Y
            Free module generated by {2} over Rational Field
            sage: Y.category()
            Join of Category of finite dimensional vector spaces with basis over Rational Field and Category of quotients of sets
            sage: Y.basis().list()
            [y[2]]
            sage: TestSuite(Y).run()
        """
        self._submodule = submodule
        self._ambient = submodule.ambient()
        embedding = submodule.lift
        indices = embedding.cokernel_basis_indices()
        CombinatorialFreeModule.__init__(self,
                                         submodule.base_ring(), indices,
                                         category=category)

    def ambient(self):
        r"""
        Return the ambient space of ``self``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: Y = X.quotient_module((x[0]-x[1], x[1]-x[2]))
            sage: Y.ambient() is X
            True
        """
        return self._ambient

    def lift(self, x):
        r"""
        Lift ``x`` to the ambient space of ``self``.

        INPUT:

        - ``x`` -- an element of ``self``

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: Y = X.quotient_module((x[0]-x[1], x[1]-x[2]));         y = Y.basis()
            sage: Y.lift(y[2])
            x[2]
        """
        assert x in self
        return self._ambient._from_dict(x._monomial_coefficients)

    def retract(self, x):
        r"""
        Retract an element of the ambient space by projecting it back to ``self``.

        INPUT:

        - ``x`` -- an element of the ambient space of ``self``

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: Y = X.quotient_module((x[0]-x[1], x[1]-x[2]));         y = Y.basis()
            sage: Y.print_options(prefix='y')
            sage: Y.retract(x[0])
            y[2]
            sage: Y.retract(x[1])
            y[2]
            sage: Y.retract(x[2])
            y[2]
        """
        return self._from_dict(self._submodule.reduce(x)._monomial_coefficients)


class SubmoduleWithBasis(CombinatorialFreeModule):
    r"""
    A base class for submodules of a ModuleWithBasis spanned by a
    (possibly infinite) basis in echelon form.

    INPUT:

    - ``basis`` -- a family of elements in echelon form in some
      :class:`module with basis <ModulesWithBasis>` `V`, or data that
      can be converted into such a family

    - ``support_order`` -- an ordering of the support of ``basis``
      expressed in ``ambient`` given as a list

    - ``unitriangular`` -- if the lift morphism is unitriangular

    - ``ambient`` -- the ambient space `V`

    - ``category`` -- a category

    Further arguments are passed down to
    :class:`CombinatorialFreeModule`.

    This is meant to be constructed via
    :meth:`Modules.WithBasis.ParentMethods.submodule`.

    .. SEEALSO::

        - :meth:`Modules.WithBasis.ParentMethods.submodule`
        - :class:`QuotientModuleWithBasis`
    """
    @staticmethod
    def __classcall_private__(cls, basis, support_order, ambient=None,
                              unitriangular=False, category=None, *args, **opts):
        r"""
        Normalize the input.

        TESTS::

            sage: from sage.modules.with_basis.subquotient import SubmoduleWithBasis
            sage: X = CombinatorialFreeModule(QQ, range(3)); x = X.basis()
            sage: Y1 = SubmoduleWithBasis((x[0]-x[1], x[1]-x[2]), [0,1,2], X)
            sage: Y2 = SubmoduleWithBasis([x[0]-x[1], x[1]-x[2]], (0,1,2), X)
            sage: Y1 is Y2
            True
        """
        basis = Family(basis)
        if ambient is None:
            ambient = basis.an_element().parent()
        Mod = ModulesWithBasis(ambient.category().base_ring())
        default_category = Mod.Subobjects()
        # Submodules of filtered modules are always canonically filtered.
        # We add this to the category if it has not been explicitly passed.
        if category is None and ambient.category().is_subcategory(Mod.Filtered()):
            default_category = default_category.Filtered()
        category = default_category.or_subcategory(category, join=True)
        return super().__classcall__(cls,
                    basis, tuple(support_order), ambient, unitriangular, category,
                    *args, **opts)

    def __init__(self, basis, support_order, ambient, unitriangular, category,
                 *args, **opts):
        r"""
        Initialization.

        TESTS::

            sage: from sage.modules.with_basis.subquotient import SubmoduleWithBasis
            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: ybas = (x[0]-x[1], x[1]-x[2])
            sage: Y = SubmoduleWithBasis(ybas, [0, 1, 2], X)
            sage: Y.print_options(prefix='y')
            sage: Y.basis().list()
            [y[0], y[1]]
            sage: [ y.lift() for y in Y.basis() ]
            [x[0] - x[1], x[1] - x[2]]
            sage: TestSuite(Y).run()
        """
        ring = ambient.base_ring()
        CombinatorialFreeModule.__init__(self, ring, basis.keys(),
                                         category=category.Subobjects(),
                                         *args, **opts)
        self._ambient = ambient
        self._basis = basis
        self._unitriangular = unitriangular
        self._support_order = support_order
        self.lift_on_basis = self._basis.__getitem__
        self.lift.register_as_coercion()

    def ambient(self):
        """
        Return the ambient space of ``self``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3)); x = X.basis()
            sage: Y = X.submodule((x[0]-x[1], x[1]-x[2]))
            sage: Y.ambient() is X
            True
        """
        return self._ambient

    @cached_method
    def _support_key(self, x):
        """
        Return a key corresponding to the index ``x`` for ordering the
        basis of ``self``.

        EXAMPLES::

            sage: A = GradedModulesWithBasis(ZZ).example()
            sage: M = A.submodule(list(A.basis(3)), already_echelonized=True)
            sage: [M._support_key(x) for x in M._support_order]
            [0, 1, 2]
        """
        return self._support_order.index(x)

    @lazy_attribute
    def lift(self):
        r"""
        The lift (embedding) map from ``self`` to the ambient space.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x");             x = X.basis()
            sage: Y = X.submodule((x[0]-x[1], x[1]-x[2]), already_echelonized=True); y = Y.basis()
            sage: Y.lift
            Generic morphism:
              From: Free module generated by {0, 1} over Rational Field
              To:   Free module generated by {0, 1, 2} over Rational Field
            sage: [ Y.lift(u) for u in y ]
            [x[0] - x[1], x[1] - x[2]]
            sage: (y[0] + y[1]).lift()
            x[0] - x[2]
        """
        return self.module_morphism(self.lift_on_basis,
                                    codomain=self._ambient,
                                    triangular="lower",
                                    unitriangular=self._unitriangular,
                                    key=self._support_key,
                                    inverse_on_support="compute")

    @lazy_attribute
    def reduce(self):
        r"""
        The reduce map.

        This map reduces elements of the ambient space modulo this
        submodule.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: Y = X.submodule((x[0]-x[1], x[1]-x[2]), already_echelonized=True)
            sage: Y.reduce
            Generic endomorphism of Free module generated by {0, 1, 2} over Rational Field
            sage: Y.reduce(x[1])
            x[2]
            sage: Y.reduce(2*x[0] + x[1])
            3*x[2]

        TESTS::

            sage: all( Y.reduce(u.lift()) == 0 for u in Y.basis() )
            True
        """
        return self.lift.cokernel_projection()

    @lazy_attribute
    def retract(self):
        r"""
        The retract map from the ambient space.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(3), prefix="x"); x = X.basis()
            sage: Y = X.submodule((x[0]-x[1], x[1]-x[2]), already_echelonized=True)
            sage: Y.print_options(prefix='y')
            sage: Y.retract
            Generic morphism:
              From: Free module generated by {0, 1, 2} over Rational Field
              To:   Free module generated by {0, 1} over Rational Field
            sage: Y.retract(x[0] - x[2])
            y[0] + y[1]

        TESTS::

            sage: all( Y.retract(u.lift()) == u for u in Y.basis() )
            True
        """
        return self.lift.section()

    def is_submodule(self, other):
        r"""
        Return whether ``self`` is a submodule of ``other``.

        INPUT:

        - ``other`` -- another submodule of the same ambient module
          or the ambient module itself

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: G = X.submodule([x[0]-x[2]])
            sage: H = X.submodule([x[0]-x[1], x[2]])
            sage: F.is_submodule(X)
            True
            sage: G.is_submodule(F)
            True
            sage: H.is_submodule(F)
            False
            sage: H.is_submodule(G)
            False

        Infinite dimensional examples::

            sage: X = CombinatorialFreeModule(QQ, ZZ); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: G = X.submodule([x[0]-x[2]])
            sage: H = X.submodule([x[0]-x[1]])
            sage: F.is_submodule(X)
            True
            sage: G.is_submodule(F)
            True
            sage: H.is_submodule(F)
            True
            sage: H.is_submodule(G)
            False

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: Y = CombinatorialFreeModule(QQ, range(6)); y = Y.basis()
            sage: G = Y.submodule([y[0]-y[1], y[1]-y[2], y[2]-y[3]])
            sage: F.is_submodule(G)
            Traceback (most recent call last):
            ...
            ValueError: other (=...) should be a submodule of the same ambient space
        """
        if other is self._ambient:
            return True
        if not (isinstance(self, SubmoduleWithBasis) and self.ambient() is other.ambient()):
            raise ValueError("other (=%s) should be a submodule of the same ambient space" % other)
        if self not in ModulesWithBasis.FiniteDimensional:
            raise NotImplementedError("only implemented for finite dimensional submodules")
        if self.dimension() > other.dimension():  # quick dimension check
            return False
        if not set(self._support_order) <= set(other._support_order):  # quick support check
            return False
        for b in self.basis():
            try:
                other.retract(b.lift())
            except ValueError:
                return False
        return True

    def _common_submodules(self, other):
        """
        Helper method to return a pair of submodules of the same ambient
        free modules to do the corresponding linear algebra.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-3*x[2], x[2]-5*x[3]])
            sage: G = X.submodule([x[0]-x[1], x[1]-2*x[2], x[2]-3*x[3]])
            sage: H = X.submodule([x[0]-x[1], x[1]-2*x[2], x[2]-3*x[3]], support_order=(3,2,1,0))
            sage: F._common_submodules(G)
            (Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [  1   0   0 -15]
             [  0   1   0 -15]
             [  0   0   1  -5],
             Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [ 1  0  0 -6]
             [ 0  1  0 -6]
             [ 0  0  1 -3])
            sage: H._common_submodules(F)
            (Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [   1    0    0 -1/6]
             [   0    1    0 -1/2]
             [   0    0    1   -1],
             Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [    1     0     0 -1/15]
             [    0     1     0  -1/3]
             [    0     0     1    -1])
            sage: G._common_submodules(H)
            (Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [ 1  0  0 -6]
             [ 0  1  0 -6]
             [ 0  0  1 -3],
             Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [ 1  0  0 -6]
             [ 0  1  0 -6]
             [ 0  0  1 -3])
            sage: H._common_submodules(G)
            (Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [   1    0    0 -1/6]
             [   0    1    0 -1/2]
             [   0    0    1   -1],
             Vector space of degree 4 and dimension 3 over Rational Field
             Basis matrix:
             [   1    0    0 -1/6]
             [   0    1    0 -1/2]
             [   0    0    1   -1])
        """
        from sage.modules.free_module import FreeModule
        supp_order = self._support_order
        A = FreeModule(self.base_ring(), len(supp_order))
        U = A.submodule([A([vec[supp] for supp in supp_order]) for vec in self._basis], check=False)
        V = A.submodule([A([vec[supp] for supp in supp_order]) for vec in other._basis], check=False)
        return (U, V)

    def is_equal_subspace(self, other):
        r"""
        Return whether ``self`` is an equal submodule to ``other``.

        .. NOTE::

            This is the mathematical notion of equality (as sets that are
            isomorphic as vector spaces), which is weaker than the `==`
            which takes into account things like the support order.

        INPUT:

        - ``other`` -- another submodule of the same ambient module
          or the ambient module itself

        EXAMPLES::

            sage: R.<z> = LaurentPolynomialRing(QQ)
            sage: X = CombinatorialFreeModule(R, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], z*x[1]-z*x[2], z^2*x[2]-z^2*x[3]])
            sage: G = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: H = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]], support_order=(3,2,1,0))
            sage: F.is_equal_subspace(F)
            True
            sage: F == G
            False
            sage: F.is_equal_subspace(G)
            True
            sage: F.is_equal_subspace(H)
            True
            sage: G == H  # different support orders
            False
            sage: G.is_equal_subspace(H)
            True

        ::

            sage: X = CombinatorialFreeModule(QQ, ZZ); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[3]])
            sage: G = X.submodule([x[0]-x[1], x[2]])
            sage: H = X.submodule([x[0]+x[1], x[1]+3*x[2]])
            sage: Hp = X.submodule([x[0]+x[1], x[1]+3*x[2]], prefix='Hp')
            sage: F.is_equal_subspace(X)
            False
            sage: F.is_equal_subspace(G)
            False
            sage: G.is_equal_subspace(H)
            False
            sage: H == Hp
            False
            sage: H.is_equal_subspace(Hp)
            True
        """
        if self is other:  # trivial case
            return True
        if not isinstance(self, SubmoduleWithBasis) and self.ambient() is other.ambient():
            raise ArithmeticError("other (=%s) should be a submodule of the same ambient space" % other)
        if self.dimension() != other.dimension():  # quick dimension check
            return False
        if self not in ModulesWithBasis.FiniteDimensional:
            raise NotImplementedError("only implemented for finite dimensional submodules")
        if set(self._basis) == set(other._basis):
            return True
        if set(self._support_order) != set(other._support_order):  # different supports
            return False
        U, V = self._common_submodules(other)
        return U == V

    def __add__(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: G = X.submodule([x[0]-x[2]])
            sage: H = X.submodule([x[0]-x[1], x[2]])
            sage: FG = F + G; FG
            Free module generated by {0, 1, 2} over Rational Field
            sage: [FG.lift(b) for b in FG.basis()]
            [B[0] - B[3], B[1] - B[3], B[2] - B[3]]
            sage: FH = F + H; FH
            Free module generated by {0, 1, 2, 3} over Rational Field
            sage: [FH.lift(b) for b in FH.basis()]
            [B[0], B[1], B[2], B[3]]
            sage: GH = G + H; GH
            Free module generated by {0, 1, 2} over Rational Field
            sage: [GH.lift(b) for b in GH.basis()]
            [B[0], B[1], B[2]]

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: Y = CombinatorialFreeModule(QQ, range(5)); y = Y.basis()
            sage: U = Y.submodule([y[0]-y[2]+y[3]])
            sage: F + U
            Traceback (most recent call last):
            ...
            ArithmeticError: both subspaces must have the same ambient space
            sage: F + 3
            Traceback (most recent call last):
            ...
            TypeError: both objects must be submodules
        """
        if not isinstance(other, SubmoduleWithBasis):
            raise TypeError("both objects must be submodules")
        if other.ambient() != self.ambient():
            raise ArithmeticError("both subspaces must have the same ambient space")
        return self.ambient().submodule(set(list(self._basis) + list(other._basis)), check=False)

    subspace_sum = __add__

    def __and__(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: G = X.submodule([x[0]-x[2]])
            sage: H = X.submodule([x[0]-x[1], x[2]])
            sage: FG = F & G; FG
            Free module generated by {0} over Rational Field
            sage: [FG.lift(b) for b in FG.basis()]
            [B[0] - B[2]]
            sage: FH = F & H; FH
            Free module generated by {0} over Rational Field
            sage: [FH.lift(b) for b in FH.basis()]
            [B[0] - B[1]]
            sage: GH = G & H; GH
            Free module generated by {} over Rational Field
            sage: [GH.lift(b) for b in GH.basis()]
            []

            sage: F.intersection(X) is F
            True

        TESTS::

            sage: X = CombinatorialFreeModule(QQ, range(4)); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]])
            sage: Y = CombinatorialFreeModule(QQ, range(5)); y = Y.basis()
            sage: U = Y.submodule([y[0]-y[2]+y[3]])
            sage: F & U
            Traceback (most recent call last):
            ...
            ArithmeticError: both subspaces must have the same ambient space
            sage: F & 3
            Traceback (most recent call last):
            ...
            TypeError: both objects must be submodules
        """
        if other is self._ambient:
            return self
        if not isinstance(other, SubmoduleWithBasis):
            raise TypeError("both objects must be submodules")
        if other.ambient() != self.ambient():
            raise ArithmeticError("both subspaces must have the same ambient space")
        U, V = self._common_submodules(other)
        UV = U & V  # the intersection
        A = self._ambient
        supp = self._support_order
        return A.submodule([A.element_class(A, {supp[i]: c for i, c in vec.iteritems()})
                            for vec in UV.basis()])

    intersection = __and__
    __rand__ = __and__

    def subspace(self, gens, *args, **opts):
        r"""
        The submodule of the ambient space spanned by a finite set
        of generators ``gens`` (as a submodule).

        INPUT:

        - ``gens`` -- a list or family of elements of ``self``

        For additional optional arguments, see
        :meth:`ModulesWithBasis.ParentMethods.submodule`.

        EXAMPLES::

            sage: X = CombinatorialFreeModule(QQ, range(4), prefix='X'); x = X.basis()
            sage: F = X.submodule([x[0]-x[1], x[1]-x[2], x[2]-x[3]], prefix='F'); f = F.basis()
            sage: U = F.submodule([f[0] + 2*f[1] - 5*f[2], f[1] + 2*f[2]]); U
            Free module generated by {0, 1} over Rational Field
            sage: [U.lift(u) for u in U.basis()]
            [F[0] - 9*F[2], F[1] + 2*F[2]]
            sage: V = F.subspace([f[0] + 2*f[1] - 5*f[2], f[1] + 2*f[2]]); V
            Free module generated by {0, 1} over Rational Field
            sage: [V.lift(u) for u in V.basis()]
            [X[0] - 9*X[2] + 8*X[3], X[1] + 2*X[2] - 3*X[3]]
        """
        gens = [self._ambient(g) for g in gens]
        return self._ambient.submodule(gens, *args, **opts)
