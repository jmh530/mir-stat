/++
This module contains data structures useful for specialized linear algebra
operations.

Authors: John Michael Hall

Copyright: 2024 Mir Stat Authors.

Macros:
NDSLICEREF = $(GREF_ALTTEXT mir-algorithm, $(TT $2), $2, mir, ndslice, $1)$(NBSP)
+/

module mir.math.internal.linearAlgebra.types;

static if (is(typeof({ import mir.blas; }))) {

static import cblas;
import mir.blas: Uplo;
import mir.internal.utility: Iota;
import mir.ndslice.slice: IteratorOf, LightConstOfLightScopeOf, LightImmutableOfLightConstOf, Slice, SliceKind;
import mir.qualifier: lightConst, LightConstOf, lightImmutable, LightImmutableOf, lightScope, LightScopeOf;
import mir.ndslice.traits: isMatrix;
import std.meta: staticMap;

alias Diag = cblas.Diag;

alias LabelsOf(T : Slice!(Iterator, N, kind, Labels), Iterator, size_t N, SliceKind kind, Labels...) = Labels;

/++
A view of a matrix that is assumed to be triangular.

A `TriangularView` matrix is implicitly convertible to a two-dimensional
$(NDSLICEREF slice, Slice).

The parameter `uplo` controls whether the matrix is upper or lower triangular.
For instance, with `uplo` set to `Uplo.Upper`, the upper triangular portion of
matrix is assumed to be filled with values, while the lower portion is assumed
to be filled with zeros.

The parameter `diag` controls whether the matrix is unit triangular or not.
With `diag` set to `Diag.Unit`, the matrix is assumed to be unit triangular,
while with it set to `Diag.NonUnit`, the matrix is not assumed to be.

It is the responsibility of the user to ensure that it is actually triangular
(strictly or otherwise), particularly following operations.

Params:
    uplo = controls whether the matrix is assumed to be upper or lower triangular
    diag = controls whether the matrix is assumed to be unit triangular or not
    Iterator = iterator of slice
    kind = kind of slice
    Labels = type of labels of slice
+/
struct TriangularView(Uplo uplo, Diag diag, Iterator, SliceKind kind, Labels...)
{
    import mir.ndslice.slice: Universal;

    ///
    Slice!(Iterator, 2u, kind, Labels) matrix;
    alias matrix this;

    /// Labels count.
    enum size_t L = Labels.length;

    /++
    Returns: View with stripped out reference counted context.
    The lifetime of the result mustn't be longer then the lifetime of the original slice.
    +/
    auto lightScope()() return scope @property
    {
        auto ret = Slice!(LightScopeOf!Iterator, 2u, kind, staticMap!(LightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return TriangularView!(uplo, diag, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// ditto
    auto lightScope()() @trusted return scope const @property
    {
        auto ret = Slice!(LightConstOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightConstOfLightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return TriangularView!(uplo, diag, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// ditto
    auto lightScope()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightImmutableOfLightConstOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return TriangularView!(uplo, diag, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// Returns: Mutable slice over immutable data.
    TriangularView!(uplo, diag, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightImmutable()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!Iterator, 2, kind, staticMap!(LightImmutableOf, Labels))(_structure, .lightImmutable(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightImmutable(_labels[i]);
        return typeof(return)(ret);
    }

    /// Returns: Mutable slice over const data.
    TriangularView!(uplo, diag, LightConstOf!Iterator, kind, staticMap!(LightConstOf, Labels)) lightConst()() return scope const @property @trusted
    {
        auto ret = Slice!(LightConstOf!Iterator, 2u, kind, staticMap!(LightConstOf, Labels))(_structure, .lightConst(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightConst(_labels[i]);
        return typeof(return)(ret);
    }

    /// ditto
    TriangularView!(uplo, diag, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightConst()() return scope immutable @property
    {
        return this.lightImmutable;
    }

    /++
    Transpose `TriangularView` matrix
    +/
    TriangularView!(uplo == Uplo.Upper ? Uplo.Lower : Uplo.Upper, diag, Iterator, Universal, Labels) transposed()() return scope @property
    {
        import mir.ndslice.dynamic: transposed;
        import mir.ndslice.topology: universal;

        return typeof(return)(this.matrix.universal.transposed);
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)(X);

    static assert(is(typeof(Y) == TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)));
}

/// With labels
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: slice;
    import mir.ndslice.slice: Contiguous;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = slice!(double, string, string)(3, 3);
    X[] = a;
    X.label[] = ["A", "B", "C"];
    X.label!1[] = X.label[];

    auto Y = TriangularView!(Uplo.Upper, Diag.NonUnit, double*, Contiguous, string*, string*)(X);

    static assert(is(typeof(Y) == TriangularView!(Uplo.Upper, Diag.NonUnit, double*, Contiguous, string*, string*)));
}

/// Transposed example
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous, Universal;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];
    static immutable b = [[0.010,       0,     0],
                          [0.0030, 0.0225,     0],
                          [0.006,   0.012, 0.040]];

    auto X = mininitRcslice!double(3, 3);
    auto result = mininitRcslice!double(3, 3);
    X[] = a;
    result[] = b;

    auto Y = TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)(X);
    auto Z1 = Y.transposed;
    auto Z2 = Z1.transposed;

    assert(Z1.equal(result));
    assert(Z2.equal(Y));

    static assert(is(typeof(Z1) == TriangularView!(Uplo.Lower, Diag.NonUnit, RCI!double, Universal)));
    static assert(is(typeof(Z2) == TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Universal)));
}

//
version(mir_stat_test_blas)
pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)(X);
    auto Y_LS = Y.lightScope;
    auto Y_C_LS = (cast(const)Y).lightScope;
    auto Y_I_LS = (cast(immutable)Y).lightScope;
    auto Y_LS_C_LC = (cast(const)Y_LS).lightConst;
    auto Y_LS_I_LI = (cast(immutable)Y_LS).lightImmutable;
    auto Y_LS_I_LC = (cast(immutable)Y_LS).lightConst;

    static assert(is(typeof(Y) == TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)));
    static assert(is(typeof(Y_LS) == TriangularView!(Uplo.Upper, Diag.NonUnit, double*, Contiguous)));
    static assert(is(typeof(Y_C_LS) == TriangularView!(Uplo.Upper, Diag.NonUnit, const(double)*, Contiguous)));
    static assert(is(typeof(Y_I_LS) == TriangularView!(Uplo.Upper, Diag.NonUnit, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_C_LC) == TriangularView!(Uplo.Upper, Diag.NonUnit, const(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LI) == TriangularView!(Uplo.Upper, Diag.NonUnit, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LC) == TriangularView!(Uplo.Upper, Diag.NonUnit, immutable(double)*, Contiguous)));
}

/++
Given a matrix `x`, returns a view that is assumed to be triangular.

Params:
    uplo = controls whether the matrix is assumed to be upper or lower triangular
    diag = controls whether the matrix is assumed to be unit triangular or not
+/
template asTriangular(Uplo uplo, Diag diag)
{
    /++
    Params:
        x = m(rows) x n(cols) matrix

    Result:
        m(rows) x n(cols)
    +/
    TriangularView!(uplo, diag, Iterator, kind, Labels) asTriangular(Iterator, SliceKind kind, Labels...)(Slice!(Iterator, 2, kind, Labels) x)
    {
        return typeof(return)(x);
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asTriangular!(Uplo.Upper, Diag.NonUnit);

    static assert(is(typeof(Y) == TriangularView!(Uplo.Upper, Diag.NonUnit, RCI!double, Contiguous)));
}

/++
True if type `T` can be implicitly converted to a $(LREF TriangularView), false otherwise.
+/
template isTriangular(T)
{
    enum bool isTriangular = is(T : TriangularView!(uplo, diag, Iterator, kind, Labels), Uplo uplo, Diag diag, Iterator, SliceKind kind, Labels...);
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asTriangular!(Uplo.Upper, Diag.NonUnit);

    static assert(!isTriangular!(typeof(X)));
    static assert(isTriangular!(typeof(Y)));
}

/++
A view of a square matrix that is assumed to be self-adjoint.

A `SelfAdjointView` matrix is implicitly convertible to a two-dimensional
$(NDSLICEREF slice, Slice).

For a matrix whose iterator is floating point-like, the matrix is also
symmetric. For instance in this case, if `uplo` is set to `Uplo.Upper`, then the
upper triangle of the matrix is assumed to be filled with values, while the
lower triangle can be assumed to mirror the upper.

It is the responsibility of the user to ensure that it is actually self-adjoint,
particularly following operations.

Params:
    uplo = controls whether the upper or lower triangle is stored
    Iterator = iterator
    kind = kind
    Labels = type of labels
+/
struct SelfAdjointView(Uplo uplo, Iterator, SliceKind kind, Labels...)
{
    import mir.ndslice.slice: Universal;

    ///
    Slice!(Iterator, 2u, kind, Labels) matrix;
    alias matrix this;

    /// Labels count.
    enum size_t L = Labels.length;

    /++
    Returns: View with stripped out reference counted context.
    The lifetime of the result mustn't be longer then the lifetime of the original slice.
    +/
    auto lightScope()() return scope @property
    {
        auto ret = Slice!(LightScopeOf!Iterator, 2u, kind, staticMap!(LightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return SelfAdjointView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// ditto
    auto lightScope()() @trusted return scope const @property
    {
        auto ret = Slice!(LightConstOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightConstOfLightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return SelfAdjointView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// ditto
    auto lightScope()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightImmutableOfLightConstOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return SelfAdjointView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret);
    }

    /// Returns: Mutable slice over immutable data.
    SelfAdjointView!(uplo, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightImmutable()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!Iterator, 2, kind, staticMap!(LightImmutableOf, Labels))(_structure, .lightImmutable(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightImmutable(_labels[i]);
        return typeof(return)(ret);
    }

    /// Returns: Mutable slice over const data.
    SelfAdjointView!(uplo, LightConstOf!Iterator, kind, staticMap!(LightConstOf, Labels)) lightConst()() return scope const @property @trusted
    {
        auto ret = Slice!(LightConstOf!Iterator, 2u, kind, staticMap!(LightConstOf, Labels))(_structure, .lightConst(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightConst(_labels[i]);
        return typeof(return)(ret);
    }

    /// ditto
    SelfAdjointView!(uplo, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightConst()() return scope immutable @property
    {
        return this.lightImmutable;
    }

    /++
    Transpose `SelfAdjointView` matrix.
    +/
    SelfAdjointView!(uplo == Uplo.Upper ? Uplo.Lower : Uplo.Upper, Iterator, Universal, Labels) transposed()() return scope @property
    {
        import mir.ndslice.dynamic: transposed;
        import mir.ndslice.topology: universal;

        return typeof(return)(this.matrix.universal.transposed);
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)(X);

    static assert(is(typeof(Y) == SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)));
}

/// With labels
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: slice;
    import mir.ndslice.slice: Contiguous;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = slice!(double, string, string)(3, 3);
    X[] = a;
    X.label[] = ["A", "B", "C"];
    X.label!1[] = X.label[];

    auto Y = SelfAdjointView!(Uplo.Upper, double*, Contiguous, string*, string*)(X);

    static assert(is(typeof(Y) == SelfAdjointView!(Uplo.Upper, double*, Contiguous, string*, string*)));
}

/// Transposed example
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous, Universal;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];
    static immutable b = [[0.010,       0,     0],
                          [0.0030, 0.0225,     0],
                          [0.006,   0.012, 0.040]];

    auto X = mininitRcslice!double(3, 3);
    auto result = mininitRcslice!double(3, 3);
    X[] = a;
    result[] = b;

    auto Y = SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)(X);
    auto Z1 = Y.transposed;
    auto Z2 = Z1.transposed;

    assert(Z1.equal(result));
    assert(Z2.equal(Y));

    static assert(is(typeof(Z1) == SelfAdjointView!(Uplo.Lower, RCI!double, Universal)));
    static assert(is(typeof(Z2) == SelfAdjointView!(Uplo.Upper, RCI!double, Universal)));
}

//
version(mir_stat_test_blas)
pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)(X);
    auto Y_LS = Y.lightScope;
    auto Y_C_LS = (cast(const)Y).lightScope;
    auto Y_I_LS = (cast(immutable)Y).lightScope;
    auto Y_LS_C_LC = (cast(const)Y_LS).lightConst;
    auto Y_LS_I_LI = (cast(immutable)Y_LS).lightImmutable;
    auto Y_LS_I_LC = (cast(immutable)Y_LS).lightConst;

    static assert(is(typeof(Y) == SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)));
    static assert(is(typeof(Y_LS) == SelfAdjointView!(Uplo.Upper, double*, Contiguous)));
    static assert(is(typeof(Y_C_LS) == SelfAdjointView!(Uplo.Upper, const(double)*, Contiguous)));
    static assert(is(typeof(Y_I_LS) == SelfAdjointView!(Uplo.Upper, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_C_LC) == SelfAdjointView!(Uplo.Upper, const(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LI) == SelfAdjointView!(Uplo.Upper, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LC) == SelfAdjointView!(Uplo.Upper, immutable(double)*, Contiguous)));
}

/++
Given a matrix `x`, returns a view that is assumed to be self-adjoint.

Params:
    uplo = controls whether the upper or lower triangle is assumed to be stored
+/
template asSelfAdjoint(Uplo uplo)
{
    /++
    Params:
        x = m(rows) x m(cols) matrix

    Result:
        m(rows) x m(cols)
    +/
    SelfAdjointView!(uplo, Iterator, kind, Labels) asSelfAdjoint(Iterator, SliceKind kind, Labels...)(Slice!(Iterator, 2, kind, Labels) x)
    in
    {
        assert(x.length!0 == x.length!1, "`x` must be a square matrix");
    }
    do
    {
        return typeof(return)(x);
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asSelfAdjoint!(Uplo.Upper);

    static assert(is(typeof(Y) == SelfAdjointView!(Uplo.Upper, RCI!double, Contiguous)));
}

/++
True if type `T` can be implicitly converted to a $(LREF SelfAdjointView), false otherwise.
+/
template isSelfAdjoint(T)
{
    enum bool isSelfAdjoint = is(T : SelfAdjointView!(uplo, Iterator, kind, Labels), Uplo uplo, Iterator, SliceKind kind, Labels...);
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asSelfAdjoint!(Uplo.Upper);

    static assert(isSelfAdjoint!(typeof(Y)));
    static assert(!isSelfAdjoint!(typeof(X)));
}

/++
A view of a square matrix that is assumed to be positive definite.

A `PositiveDefiniteView` matrix is implicitly convertible to a
$(LREF SelfAdjointView).

It is the responsibility of the user to ensure that it is actually positive
definite, particularly following operations.

Params:
    uplo = controls whether the upper or lower triangle is stored
    Iterator = iterator
    kind = kind
    Labels = type of labels
+/
struct PositiveDefiniteView(Uplo uplo, Iterator, SliceKind kind, Labels...)
{
    import mir.ndslice.slice: Universal;

    ///
    SelfAdjointView!(uplo, Iterator, kind, Labels) selfAdjointView;
    alias selfAdjointView this;

    /// Labels count.
    enum size_t L = Labels.length;

    /++
    Returns: View with stripped out reference counted context.
    The lifetime of the result mustn't be longer then the lifetime of the original slice.
    +/
    auto lightScope()() return scope @property
    {
        auto ret = Slice!(LightScopeOf!Iterator, 2u, kind, staticMap!(LightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return PositiveDefiniteView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret.asSelfAdjoint!uplo);
    }

    /// ditto
    auto lightScope()() @trusted return scope const @property
    {
        auto ret = Slice!(LightConstOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightConstOfLightScopeOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return PositiveDefiniteView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret.asSelfAdjoint!uplo);
    }

    /// ditto
    auto lightScope()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!(LightScopeOf!Iterator), 2u, kind, staticMap!(LightImmutableOfLightConstOf, Labels))
            (_structure, .lightScope(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightScope(_labels[i]);
        return PositiveDefiniteView!(uplo, IteratorOf!(typeof(ret)), kind, LabelsOf!(typeof(ret)))(ret.asSelfAdjoint!uplo);
    }

    /// Returns: Mutable slice over immutable data.
    PositiveDefiniteView!(uplo, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightImmutable()() return scope immutable @property
    {
        auto ret = Slice!(LightImmutableOf!Iterator, 2, kind, staticMap!(LightImmutableOf, Labels))(_structure, .lightImmutable(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightImmutable(_labels[i]);
        return typeof(return)(ret.asSelfAdjoint!uplo);
    }

    /// Returns: Mutable slice over const data.
    PositiveDefiniteView!(uplo, LightConstOf!Iterator, kind, staticMap!(LightConstOf, Labels)) lightConst()() return scope const @property @trusted
    {
        auto ret = Slice!(LightConstOf!Iterator, 2u, kind, staticMap!(LightConstOf, Labels))(_structure, .lightConst(_iterator));
        foreach(i; Iota!L)
            ret._labels[i] = .lightConst(_labels[i]);
        return typeof(return)(ret.asSelfAdjoint!uplo);
    }

    /// ditto
    PositiveDefiniteView!(uplo, LightImmutableOf!Iterator, kind, staticMap!(LightImmutableOf, Labels)) lightConst()() return scope immutable @property
    {
        return this.lightImmutable;
    }

    /++
    Transpose `PositiveDefiniteView` matrix.
    +/
    PositiveDefiniteView!(uplo == Uplo.Upper ? Uplo.Lower : Uplo.Upper, Iterator, Universal, Labels) transposed()() return scope @property
    {
        import mir.ndslice.dynamic: transposed;
        import mir.ndslice.topology: universal;

        return typeof(return)(this.selfAdjointView.matrix.universal.transposed.asSelfAdjoint!(uplo == Uplo.Upper ? Uplo.Lower : Uplo.Upper));
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)(X.asSelfAdjoint!(Uplo.Upper));

    static assert(is(typeof(Y) == PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)));
}

/// With labels
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: slice;
    import mir.ndslice.slice: Contiguous;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = slice!(double, string, string)(3, 3);
    X[] = a;
    X.label[] = ["A", "B", "C"];
    X.label!1[] = X.label[];

    auto Y = PositiveDefiniteView!(Uplo.Upper, double*, Contiguous, string*, string*)(X.asSelfAdjoint!(Uplo.Upper));

    static assert(is(typeof(Y) == PositiveDefiniteView!(Uplo.Upper, double*, Contiguous, string*, string*)));
}

/// Transposed example
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous, Universal;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];
    static immutable b = [[0.010,       0,     0],
                          [0.0030, 0.0225,     0],
                          [0.006,   0.012, 0.040]];

    auto X = mininitRcslice!double(3, 3);
    auto result = mininitRcslice!double(3, 3);
    X[] = a;
    result[] = b;

    auto Y = PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)(X.asSelfAdjoint!(Uplo.Upper));
    auto Z1 = Y.transposed;
    auto Z2 = Z1.transposed;

    assert(Z1.equal(result));
    assert(Z2.equal(Y));

    static assert(is(typeof(Z1) == PositiveDefiniteView!(Uplo.Lower, RCI!double, Universal)));
    static assert(is(typeof(Z2) == PositiveDefiniteView!(Uplo.Upper, RCI!double, Universal)));
}

//
version(mir_stat_test_blas)
pure nothrow @nogc
unittest
{
    import cblas: Uplo, Diag;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)(X.asSelfAdjoint!(Uplo.Upper));
    auto Y_LS = Y.lightScope;
    auto Y_C_LS = (cast(const)Y).lightScope;
    auto Y_I_LS = (cast(immutable)Y).lightScope;
    auto Y_LS_C_LC = (cast(const)Y_LS).lightConst;
    auto Y_LS_I_LI = (cast(immutable)Y_LS).lightImmutable;
    auto Y_LS_I_LC = (cast(immutable)Y_LS).lightConst;

    static assert(is(typeof(Y) == PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)));
    static assert(is(typeof(Y_LS) == PositiveDefiniteView!(Uplo.Upper, double*, Contiguous)));
    static assert(is(typeof(Y_C_LS) == PositiveDefiniteView!(Uplo.Upper, const(double)*, Contiguous)));
    static assert(is(typeof(Y_I_LS) == PositiveDefiniteView!(Uplo.Upper, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_C_LC) == PositiveDefiniteView!(Uplo.Upper, const(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LI) == PositiveDefiniteView!(Uplo.Upper, immutable(double)*, Contiguous)));
    static assert(is(typeof(Y_LS_I_LC) == PositiveDefiniteView!(Uplo.Upper, immutable(double)*, Contiguous)));
}

/++
Given a matrix `x`, returns a view that is assumed to be positive definite.

Params:
    uplo = controls whether the upper or lower triangle is assumed to be stored
+/
template asPositiveDefinite(Uplo uplo)
{
    /++
    Params:
        x = m(rows) x m(cols) matrix

    Result:
        m(rows) x m(cols)
    +/
    PositiveDefiniteView!(uplo, Iterator, kind, Labels) asPositiveDefinite(Iterator, SliceKind kind, Labels...)(Slice!(Iterator, 2, kind, Labels) x)
    in
    {
        assert(x.length!0 == x.length!1, "`x` must be a square matrix");
    }
    do
    {
        return typeof(return)(x.asSelfAdjoint!uplo);
    }

    /++
    Params:
        x = m(rows) x m(cols) matrix

    Result:
        m(rows) x m(cols)
    +/
    PositiveDefiniteView!(uplo, Iterator, kind, Labels) asPositiveDefinite(Iterator, SliceKind kind, Labels...)(SelfAdjointView!(uplo, Iterator, kind, Labels) x)
    {
        return typeof(return)(x);
    }
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: Contiguous;
    import mir.rc.array: RCI;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asPositiveDefinite!(Uplo.Upper);
    auto Z = X.asSelfAdjoint!(Uplo.Upper).asPositiveDefinite!(Uplo.Upper);

    static assert(is(typeof(Y) == PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)));
    static assert(is(typeof(Z) == PositiveDefiniteView!(Uplo.Upper, RCI!double, Contiguous)));
}

/++
True if type `T` can be implicitly converted to a $(LREF PositiveDefiniteView), false otherwise.
+/
template isPositiveDefinite(T)
{
    enum bool isPositiveDefinite = is(T : PositiveDefiniteView!(uplo, Iterator, kind, Labels), Uplo uplo, Iterator, SliceKind kind, Labels...);
}

///
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import cblas: Uplo;
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [[0.010, 0.0030, 0.006],
                          [0,     0.0225, 0.012],
                          [0,     0,      0.040]];

    auto X = mininitRcslice!double(3, 3);
    X[] = a;

    auto Y = X.asPositiveDefinite!(Uplo.Upper);

    static assert(isPositiveDefinite!(typeof(Y)));
    static assert(!isPositiveDefinite!(typeof(X)));
}

}
