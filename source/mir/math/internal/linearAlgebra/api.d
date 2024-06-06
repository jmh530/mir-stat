/++
This module provides a number of functions used for linear algebra.

Authors: John Michael Hall

Copyright: 2024 Mir Stat Authors.

+/

module mir.math.internal.linearAlgebra.api;

static if (is(typeof({ import mir.blas; }))) {

static import cblas;
import mir.blas: Uplo;
import mir.internal.utility: isComplex, isFloatingPoint;
import mir.ndslice.slice: Slice, SliceKind;
import mir.rc.array: RCI;
import std.traits: Unqual;

alias Diag = cblas.Diag;

/++
Matrix multiplication. Allocates result to using Mir refcounted arrays.

Additional overloads are provided for symmetric matrices.

This function has multiple overloads that include the following functionality:
- a.mtimes(b) where `a` and `b` are both two-dimensional slices. The result is a
two-dimensional slice.
- a.mtimes(b) where `a` is a two-dimensional slice and `b` is a one-dimensional
slice. The result is a one-dimensional slice. In this case, `b` can be thought
of as a column vector.
- b.mtimes(a) where `a` is a two-dimensional slice and `b` is a one-dimensional
slice. The result is a one-dimensional slice. In this case, `b` can be thought
of as a row vector.

Params:
    a = m(rows) x k(cols) matrix
    b = k(rows) x n(cols) matrix
Result:
    m(rows) x n(cols)
+/
@safe pure nothrow @nogc
Slice!(RCI!T, 2) mtimes(T, SliceKind kindA, SliceKind kindB)(
    Slice!(const(T)*, 2, kindA) a,
    Slice!(const(T)*, 2, kindB) b
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
}
out (c)
{
    assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
    assert(c.length!1 == b.length!1, "The second dimension of the result must match the second dimension of `b`");
}
do
{
    import mir.math.internal.linearAlgebra.kernel: mtimesGeneralKernel;
    import mir.ndslice.allocation: mininitRcslice;

    auto c = mininitRcslice!T(a.length!0, b.length!1);
    mtimesGeneralKernel(a, b, c.lightScope);
    return c;
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 2, kindA) a,
    auto ref const Slice!(RCI!B, 2, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    auto scopeB = b.lightScope.lightConst;
    return mtimes(scopeA, scopeB);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 2, kindA) a,
    Slice!(const(B)*, 2, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    return mtimes(scopeA, b);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    Slice!(const(A)*, 2, kindA) a,
    auto ref const Slice!(RCI!B, 2, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
}
do
{
    auto scopeB = b.lightScope.lightConst;
    return mtimes(a, scopeB);
}

/++
Params:
    a = m(rows) x n(cols) matrix
    b = n(rows) x 1(cols) vector
Result:
    m(rows) x 1(cols)
+/
@safe pure nothrow @nogc
Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
    Slice!(const(T)*, 2, kindA) a,
    Slice!(const(T)*, 1, kindB) b
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
}
out (c)
{
    assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
}
do
{
    import mir.math.internal.linearAlgebra.kernel: mtimesGeneralKernel;
    import mir.ndslice.allocation: mininitRcslice;

    auto c = mininitRcslice!T(a.length!0);
    mtimesGeneralKernel(a, b, c.lightScope);
    return c;
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 2, kindA) a,
    auto ref const Slice!(RCI!B, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    auto scopeB = b.lightScope.lightConst;
    return mtimes(scopeA, scopeB);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 2, kindA) a,
    Slice!(const(B)*, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    return mtimes(scopeA, b);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    Slice!(const(A)*, 2, kindA) a,
    auto ref const Slice!(RCI!B, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
}
do
{
    auto scopeB = b.lightScope.lightConst;
    return mtimes(a, scopeB);
}

/++
Params:
    a = 1(rows) x n(cols) vector
    b = n(rows) x m(cols) matrix
Result:
    1(rows) x m(cols)
+/
@safe pure nothrow @nogc
Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
    Slice!(const(T)*, 1, kindB) a,
    Slice!(const(T)*, 2, kindA) b
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length == b.length!0, "The length of `a` must match the first dimension of `b`");
}
out (c)
{
    assert(c.length!0 == b.length!1, "The first dimension of the result must match the second dimension of `b`");
}
do
{
    import mir.math.internal.linearAlgebra.kernel: mtimesGeneralKernel;
    import mir.ndslice.allocation: mininitRcslice;

    auto c = mininitRcslice!T(b.length!1);
    mtimesGeneralKernel(a, b, c.lightScope);
    return c;
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!B, 1, kindB) a,
    auto ref const Slice!(RCI!A, 2, kindA) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length == b.length!0, "The length of `a` must match the first dimension of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    auto scopeB = b.lightScope.lightConst;
    return mtimes(scopeA, scopeB);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    Slice!(const(A)*, 1, kindA) a,
    auto ref const Slice!(RCI!B, 2, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length == b.length!0, "The length of `a` must match the first dimension of `b`");
}
do
{
    auto scopeB = b.lightScope.lightConst;
    return mtimes(a, scopeB);
}

/// ditto
@safe pure nothrow @nogc
Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 1, kindA) a,
    Slice!(const(B)*, 2, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length == b.length!0, "The length of `a` must match the first dimension of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    return mtimes(scopeA, b);
}

/++
Params:
    a = 1(rows) x n(cols) vector
    b = n(rows) x 1(cols) vector
Result:
    dot product
+/
@safe pure nothrow @nogc
Unqual!A mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    Slice!(const(A)*, 1, kindA) a,
    Slice!(const(B)*, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the length of `b`");
}
do
{
    import mir.math.internal.linearAlgebra.kernel: mtimesGeneralKernel;

    typeof(return) c = void;
    mtimesGeneralKernel(a, b, c);
    return c;
}

/// ditto
@safe pure nothrow @nogc
Unqual!A mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 1, kindA) a,
    auto ref const Slice!(RCI!B, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the length of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    auto scopeB = b.lightScope.lightConst;
    return mtimes(scopeA, scopeB);
}

/// ditto
@safe pure nothrow @nogc
Unqual!A mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    auto ref const Slice!(RCI!A, 1, kindA) a,
    Slice!(const(B)*, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the length of `b`");
}
do
{
    auto scopeA = a.lightScope.lightConst;
    return mtimes(scopeA, b);
}

/// ditto
@safe pure nothrow @nogc
Unqual!A mtimes(A, B, SliceKind kindA, SliceKind kindB)(
    Slice!(const(A)*, 1, kindA) a,
    auto ref const Slice!(RCI!B, 1, kindB) b
)
    if (is(Unqual!A == Unqual!B))
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the length of `b`");
}
do
{
    auto scopeB = b.lightScope.lightConst;
    return mtimes(a, scopeB);
}

/// Matrix-matrix multiplication (real)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.ndslice.dynamic: transposed;
    import mir.ndslice.allocation: mininitRcslice;

    auto a = mininitRcslice!double(3, 5);
    auto b = mininitRcslice!double(5, 4);

    a[] =
    [[-5,  1,  7, 7, -4],
     [-1, -5,  6, 3, -3],
     [-5, -2, -3, 6,  0]];

    b[] =
    [[-5, -3,  3,  1],
     [ 4,  3,  6,  4],
     [-4, -2, -2,  2],
     [-1,  9,  4,  8],
     [ 9,  8,  3, -2]];

    assert(mtimes(a, b) ==
        [[-42,  35,  -7, 77],
         [-69, -21, -42, 21],
         [ 23,  69,   3, 29]]);

    assert(mtimes(b.transposed, a.transposed) ==
        [[-42, -69, 23],
         [ 35, -21, 69],
         [ -7, -42,  3],
         [ 77,  21, 29]]);
}

// test mixed strides
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.ndslice.dynamic: transposed;
    import mir.ndslice.allocation: mininitRcslice;

    auto a = mininitRcslice!double(3, 5);
    auto b = mininitRcslice!double(5, 4);

    a[] =
    [[-5,  1,  7, 7, -4],
     [-1, -5,  6, 3, -3],
     [-5, -2, -3, 6,  0]];

    b[] =
    [[-5, -3,  3,  1],
     [ 4,  3,  6,  4],
     [-4, -2, -2,  2],
     [-1,  9,  4,  8],
     [ 9,  8,  3, -2]];

    auto at = mininitRcslice!double(5, 3);
    at[] =
    [[-5, -1, -5],
     [ 1, -5, -2],
     [ 7,  6, -3],
     [ 7,  3,  6],
     [-4, -3,  0]];
     assert(mtimes(b.transposed, at) ==
        [[-42, -69, 23],
         [ 35, -21, 69],
         [ -7, -42,  3],
         [ 77,  21, 29]]);

    auto bt = mininitRcslice!double(4, 5);
    bt[] =
    [[-5, 4, -4, -1,  9],
     [-3, 3, -2,  9,  8],
     [ 3, 6, -2,  4,  3],
     [ 1, 4,  2,  8, -2]];
     assert(mtimes(bt, a.transposed) ==
        [[-42, -69, 23],
         [ 35, -21, 69],
         [ -7, -42,  3],
         [ 77,  21, 29]]);
}

/// Matrix-matrix multiplication (complex)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: mininitRcslice;
    import mir.complex: Complex;

    auto a = mininitRcslice!(Complex!double)(3, 5);
    auto b = mininitRcslice!(Complex!double)(5, 4);

    a[] =
    [[-5,  1,  7, 7, -4],
     [-1, -5,  6, 3, -3],
     [-5, -2, -3, 6,  0]];

    b[] =
    [[-5, -3,  3,  1],
     [ 4,  3,  6,  4],
     [-4, -2, -2,  2],
     [-1,  9,  4,  8],
     [ 9, 8,  3, -2]];

    assert(mtimes(a, b) ==
        [[-42,  35,  -7, 77],
         [-69, -21, -42, 21],
         [ 23,  69,   3, 29]]);
}

/// Matrix-matrix multiplication, specialization for MxN times Nx1
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2, -3], [-2.0, 2, 3, 10], [0.0, 2, 1, 1]];
    static immutable b = [2.0, 3, 4, 5];
    static immutable c = [14.0, 64, 15];

    auto X = mininitRcslice!double(3, 4);
    auto y = mininitRcslice!double(4);
    auto result = mininitRcslice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    auto Xy = X.mtimes(y);
    assert(Xy.equal(result));
    auto yXT = y.mtimes(X.transposed);
    assert(yXT.equal(result));
}

/// Reference-counted dot product
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [-5.0,  1,  7,  7, -4];
    static immutable b = [ 4.0, -4, -2, 10,  4];

    auto x = mininitRcslice!double(5);
    auto y = mininitRcslice!double(5);

    x[] = a;
    y[] = b;

    assert(mtimes(x, y) == 16);
}

/// Mix slice & RC dot product
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.slice: sliced;

    static immutable a = [-5.0,  1,  7,  7, -4];
    static immutable b = [ 4.0, -4, -2, 10,  4];

    auto x = mininitRcslice!double(5);
    auto y = b.sliced;

    x[] = a;

    assert(mtimes(x, y) == 16);
    assert(mtimes(y, x) == 16);
}

/++
Similar to above, but allows for inputs to be symmetric.

Params:
    uplo = controls whether matrix is upper symmetric or lower symmetric
+/
template mtimes(Uplo uplo = Uplo.Upper)
{
    import mir.math.internal.linearAlgebra.types: SelfAdjointView;

    /+
    Params:
        a = m(rows) x m(cols) symmetric matrix
        b = m(rows) x n(cols) matrix
    Result:
        m(rows) x n(cols)
    +/
    Slice!(RCI!T, 2) mtimes(T, SliceKind kindA, SliceKind kindB)(
        SelfAdjointView!(uplo, const(T)*, kindA) a,
        Slice!(const(T)*, 2, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    out (c)
    {
        assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of the result must match the second dimension of `b`");
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesSymmetricKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0, b.length!1);
        mtimesSymmetricKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const SelfAdjointView!(uplo, RCI!A, kindA) a,
        auto ref const Slice!(RCI!B, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const SelfAdjointView!(uplo, RCI!A, kindA) a,
        Slice!(const(B)*, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        SelfAdjointView!(uplo, const(A)*, kindA) a,
        auto ref const Slice!(RCI!B, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /++
    Params:
        a = m(rows) x m(cols) symmetric matrix
        b = m(rows) x 1(cols) vector
    Result:
        m(rows) x 1(cols)
    +/
    @safe pure nothrow @nogc
    Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
        SelfAdjointView!(uplo, const(T)*, kindA) a,
        Slice!(const(T)*, 1, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    out (c)
    {
        assert(c.length == a.length);
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesSymmetricKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0);
        mtimesSymmetricKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const SelfAdjointView!(uplo, RCI!A, kindA) a,
        auto ref const Slice!(RCI!B, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const SelfAdjointView!(uplo, RCI!A, kindA) a,
        Slice!(const(B)*, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        SelfAdjointView!(uplo, const(A)*, kindA) a,
        auto ref const Slice!(RCI!B, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /+
    Params:
        a = m(rows) x n(cols) matrix
        b = n(rows) x n(cols) symmetric matrix
    Result:
        m(rows) x n(cols)
    +/
    Slice!(RCI!T, 2) mtimes(T, SliceKind kindA, SliceKind kindB)(
        Slice!(const(T)*, 2, kindA) a,
        SelfAdjointView!(uplo, const(T)*, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    out (c)
    {
        assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of the result must match the second dimension of `b`");
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesSymmetricRightKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0, b.length!1);
        mtimesSymmetricRightKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 2, kindA) a,
        auto ref const SelfAdjointView!(uplo, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 2, kindA) a,
        SelfAdjointView!(uplo, const(B)*, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        Slice!(const(A)*, 2, kindA) a,
        auto ref const SelfAdjointView!(uplo, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /++
    Params:
        a = 1(rows) x m(cols) vector
        b = m(rows) x m(cols) symmetric matrix
    Result:
        m(rows) x 1(cols)
    +/
    @safe pure nothrow @nogc
    Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
        Slice!(const(T)*, 1, kindA) a,
        SelfAdjointView!(uplo, const(T)*, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    out (c)
    {
        assert(c.length == a.length);
    }
    do
    {
        return mtimes(b, a);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 1, kindA) a,
        auto ref const SelfAdjointView!(uplo, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 1, kindA) a,
        SelfAdjointView!(uplo, const(B)*, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        Slice!(const(A)*, 1, kindA) a,
        auto ref const SelfAdjointView!(uplo, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }
}

/// Symmetric Matrix-Matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asSelfAdjoint;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable b = [[2.0, 3], [4.0, 3], [0.0, -5]];
    static immutable c = [[26.0, 14], [18.0, 6], [16.0, 10]];

    auto X = mininitRcslice!double(3, 3);
    auto Y = mininitRcslice!double(3, 2);
    auto result = mininitRcslice!double(3, 2);

    X[] = a;
    Y[] = b;
    result[] = c;

    auto XY = X.asSelfAdjoint!(Uplo.Upper).mtimes(Y);
    assert(XY.equal(result));
    auto YTX = Y.transposed.mtimes(X.asSelfAdjoint!(Uplo.Upper));
    assert(YTX.equal(result.transposed));
}

/// Symmetric Matrix, specialization for vectors
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asSelfAdjoint;
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 28, 17];

    auto X = mininitRcslice!double(3, 3);
    auto y = mininitRcslice!double(3);
    auto result = mininitRcslice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    auto Xy = X.asSelfAdjoint!(Uplo.Upper).mtimes(y);
    assert(Xy.equal(result));
    auto yX = y.mtimes(X.asSelfAdjoint!(Uplo.Upper));
    assert(yX.equal(result));
}

/// Symmetric Matrix, specialization for vectors (GC version)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asSelfAdjoint;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 28, 17];

    auto X = uninitSlice!double(3, 3);
    auto y = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    auto Xy = X.asSelfAdjoint!(Uplo.Upper).mtimes(y);
    assert(Xy.equal(result));
    auto yX = y.mtimes(X.asSelfAdjoint!(Uplo.Upper));
    assert(yX.equal(result));
}

/++
Similar to above, but allows for inputs to be triangular.

Params:
    uplo = controls whether `a` is upper triangular or lower triangular
    diag = controls whether `a` is is non-unit triangular or unit triangular
+/
template mtimes(Uplo uplo = Uplo.Upper, Diag diag = Diag.NonUnit)
{
    import mir.math.internal.linearAlgebra.types: TriangularView;

    /+
    Params:
        a = m(rows) x m(cols) triangular matrix
        b = m(rows) x n(cols) matrix
    Result:
        m(rows) x n(cols)
    +/
    Slice!(RCI!T, 2) mtimes(T, SliceKind kindA, SliceKind kindB)(
        TriangularView!(uplo, diag, const(T)*, kindA) a,
        Slice!(const(T)*, 2, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    out (c)
    {
        assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of the result must match the second dimension of `b`");
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesTriangularKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0, b.length!1);
        mtimesTriangularKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const TriangularView!(uplo, diag, RCI!A, kindA) a,
        auto ref const Slice!(RCI!B, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const TriangularView!(uplo, diag, RCI!A, kindA) a,
        Slice!(const(B)*, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        TriangularView!(uplo, diag, const(A)*, kindA) a,
        auto ref const Slice!(RCI!B, 2, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /++
    Params:
        a = m(rows) x m(cols) triangular matrix
        b = m(rows) x 1(cols) vector
    Result:
        m(rows) x 1(cols)
    +/
    @safe pure nothrow @nogc
    Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
        TriangularView!(uplo, diag, const(T)*, kindA) a,
        Slice!(const(T)*, 1, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    out (c)
    {
        assert(c.length == a.length);
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesTriangularKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0);
        mtimesTriangularKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const TriangularView!(uplo, diag, RCI!A, kindA) a,
        auto ref const Slice!(RCI!B, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const TriangularView!(uplo, diag, RCI!A, kindA) a,
        Slice!(const(B)*, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        TriangularView!(uplo, diag, const(A)*, kindA) a,
        auto ref const Slice!(RCI!B, 1, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /+
    Params:
        a = m(rows) x n(cols) matrix
        b = n(rows) x n(cols) triangular matrix
    Result:
        m(rows) x n(cols)
    +/
    Slice!(RCI!T, 2) mtimes(T, SliceKind kindA, SliceKind kindB)(
        Slice!(const(T)*, 2, kindA) a,
        TriangularView!(uplo, diag, const(T)*, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    out (c)
    {
        assert(c.length!0 == a.length!0, "The first dimension of the result must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of the result must match the second dimension of `b`");
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesTriangularRightKernel;
        import mir.ndslice.allocation: mininitRcslice;
        import std.stdio: writeln;
        debug writeln("in this function");
        auto c = mininitRcslice!T(a.length!0, b.length!1);
        mtimesTriangularRightKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 2, kindA) a,
        auto ref const TriangularView!(uplo, diag, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 2, kindA) a,
        TriangularView!(uplo, diag, const(B)*, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 2) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        Slice!(const(A)*, 2, kindA) a,
        auto ref const TriangularView!(uplo, diag, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }

    /++
    Params:
        a = 1(rows) x m(cols) vector
        b = m(rows) x m(cols) triangular matrix
    Result:
        m(rows) x 1(cols)
    +/
    @safe pure nothrow @nogc
    Slice!(RCI!T, 1) mtimes(T, SliceKind kindA, SliceKind kindB)(
        Slice!(const(T)*, 1, kindA) a,
        TriangularView!(uplo, diag, const(T)*, kindB) b
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    out (c)
    {
        assert(c.length == a.length);
    }
    do
    {
        import mir.math.internal.linearAlgebra.kernel: mtimesTriangularRightKernel;
        import mir.ndslice.allocation: mininitRcslice;

        auto c = mininitRcslice!T(a.length!0);
        mtimesTriangularRightKernel(a, b, c.lightScope);
        return c;
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 1, kindA) a,
        auto ref const TriangularView!(uplo, diag, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        auto scopeB = b.lightScope.lightConst;
        return mtimes(scopeA, scopeB);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        auto ref const Slice!(RCI!A, 1, kindA) a,
        TriangularView!(uplo, diag, const(B)*, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeA = a.lightScope.lightConst;
        return mtimes(scopeA, b);
    }

    /// ditto
    @safe pure nothrow @nogc
    Slice!(RCI!(Unqual!A), 1) mtimes(A, B, SliceKind kindA, SliceKind kindB)(
        Slice!(const(A)*, 1, kindA) a,
        auto ref const TriangularView!(uplo, diag, RCI!B, kindB) b
    )
        if (is(Unqual!A == Unqual!B))
    in
    {
        assert(a.length == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
    }
    do
    {
        auto scopeB = b.lightScope.lightConst;
        return mtimes(a, scopeB);
    }
}

/// Triangular Matrix-Matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asTriangular;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable b = [[2.0, 3], [4.0, 3], [0.0, -5]];
    static immutable c = [[26.0, 14], [8.0, -9], [0.0, -5]];

    auto X = mininitRcslice!double(3, 3);
    auto Y = mininitRcslice!double(3, 2);
    auto result = mininitRcslice!double(3, 2);

    X[] = a;
    Y[] = b;
    result[] = c;

    auto XY = X.asTriangular!(Uplo.Upper, Diag.NonUnit).mtimes(Y);
    assert(XY.equal(result));

    auto Z = X.asTriangular!(Uplo.Upper, Diag.NonUnit).transposed;

    auto YTX = Y.transposed.mtimes(X.asTriangular!(Uplo.Upper, Diag.NonUnit).transposed);
    assert(YTX.equal(result.transposed));
}

/// Triangular Matrix, specialization for vectors
version(mir_stat_test_blas)
@safe pure nothrow @nogc
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asTriangular;
    import mir.ndslice.allocation: mininitRcslice;

    static immutable a = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 18, 4];

    auto X = mininitRcslice!double(3, 3);
    auto y = mininitRcslice!double(3);
    auto result = mininitRcslice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    auto Xy = X.asTriangular!(Uplo.Upper, Diag.NonUnit).mtimes(y);
    assert(Xy.equal(result));
    auto yX = y.mtimes(X.asTriangular!(Uplo.Upper, Diag.NonUnit).transposed);
    assert(yX.equal(result));
}

/// Triangular Matrix, specialization for vectors (GC version)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.math.internal.linearAlgebra.types: asTriangular;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 18, 4];

    auto X = uninitSlice!double(3, 3);
    auto y = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    auto Xy = X.asTriangular!(Uplo.Upper, Diag.NonUnit).mtimes(y);
    assert(Xy.equal(result));
    auto yX = y.mtimes(X.asTriangular!(Uplo.Upper, Diag.NonUnit).transposed);
    assert(yX.equal(result));
}

}
