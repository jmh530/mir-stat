/++
This module contains kernels used in linear algebra. For the purposes of this
module, each function takes `Slice`s and has a separate name, rather than
separate types and overloads.

Authors: John Michael Hall

Copyright: 2024 Mir Stat Authors.

+/

module mir.math.internal.linearAlgebra.kernel;

static if (is(typeof({ import mir.blas; }))) {

static import cblas;
import mir.blas: Uplo;
import mir.internal.utility: isComplex, isFloatingPoint;
import mir.ndslice.slice: Slice, SliceKind;

alias Diag = cblas.Diag;

/++
General matrix multiplication.

Params:
    a = m(rows) x n(cols) matrix
    b = n(rows) x p(cols) matrix
    c = m(rows) x p(cols) matrix
+/
@safe pure nothrow @nogc
void mtimesGeneralKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
    Slice!(const(T)*, 2, kindA) a,
    Slice!(const(T)*, 2, kindB) b,
    Slice!(T*, 2, kindC) c
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
    assert(c.length!0 == a.length!0, "The first dimension of `c` must match the first dimension of `a`");
    assert(c.length!1 == b.length!1, "The second dimension of `c` must match the second dimension of `b`");
}
do
{
    import mir.blas: gemm;

    gemm(cast(T)1, a, b, cast(T)0, c);
}

/+
Params:
    a = m(rows) x n(cols) matrix
    b = n-dimensional vector
    c = m-dimensional vector
+/
@safe pure nothrow @nogc
void mtimesGeneralKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
    Slice!(const(T)*, 2, kindA) a,
    Slice!(const(T)*, 1, kindB) b,
    Slice!(T*, 1, kindC) c
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
    assert(c.length!0 == a.length!0, "The length of `c` must match the first dimension of `a`");
}
do
{
    import mir.blas: gemv;

    gemv(cast(T)1, a, b, cast(T)0, c);
}

/+
Params:
    a = m-dimensional vector
    b = m(rows) x n(cols) matrix
    c = n-dimensional vector
+/
@safe pure nothrow @nogc
void mtimesGeneralKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
    Slice!(const(T)*, 1, kindB) a,
    Slice!(const(T)*, 2, kindA) b,
    Slice!(T*, 1, kindC) c
)
    if (isFloatingPoint!T || isComplex!T)
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the first dimension of `b`");
    assert(c.length!0 == b.length!1, "The length of `c` must match the second dimension  of `b`");
}
do
{
    import mir.ndslice.dynamic: transposed;
    import mir.ndslice.topology: universal;

    mtimesGeneralKernel(b.universal.transposed, a, c);
}

/++
Params:
    a = 1(rows) x m(cols) vector
    b = m(rows) x 1(cols) vector
    c = scaler
+/
@safe pure nothrow @nogc
void mtimesGeneralKernel(T, SliceKind kindA, SliceKind kindB)(
    Slice!(const(T)*, 1, kindA) a,
    Slice!(const(T)*, 1, kindB) b,
    out T c
)
    if (isFloatingPoint!T)
in
{
    assert(a.length!0 == b.length!0, "The length of `a` must match the length of `b`");
}
do
{
    import mir.blas: dot;

    c = dot(a, b);
}

/// General matrix multiplication (real)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[-5.0,  1,  7, 7, -4],
                          [-1.0, -5,  6, 3, -3],
                          [-5.0, -2, -3, 6,  0]];
    static immutable b = [[-5.0, -3,  3,  1],
                          [ 4.0,  3,  6,  4],
                          [-4.0, -2, -2,  2],
                          [-1.0,  9,  4,  8],
                          [ 9.0, 8,  3, -2]];
    static immutable c = [[-42.0,  35,  -7, 77],
                          [-69.0, -21, -42, 21],
                          [ 23.0,  69,   3, 29]];

    auto X = uninitSlice!double(3, 5);
    auto Y = uninitSlice!double(5, 4);
    auto XY = uninitSlice!double(3, 4);
    auto YtXt = uninitSlice!double(4, 3);
    auto result = uninitSlice!double(3, 4);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesGeneralKernel(Y, XY);
    assert(XY.equal(result));
    Y.transposed.mtimesGeneralKernel(X.transposed, YtXt);
    assert(YtXt.equal(result.transposed));
}

// test mixed strides
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[-5.0,  1,  7, 7, -4],
                          [-1.0, -5,  6, 3, -3],
                          [-5.0, -2, -3, 6,  0]];
    static immutable at = [[-5.0, -1, -5],
                           [ 1.0, -5, -2],
                           [ 7.0,  6, -3],
                           [ 7.0,  3,  6],
                           [-4.0, -3,  0]];
    static immutable b = [[-5.0, -3,  3,  1],
                          [ 4.0,  3,  6,  4],
                          [-4.0, -2, -2,  2],
                          [-1.0,  9,  4,  8],
                          [ 9.0, 8,  3, -2]];
    static immutable bt = [[-5.0, 4, -4, -1,  9],
                           [-3.0, 3, -2,  9,  8],
                           [ 3.0, 6, -2,  4,  3],
                           [ 1.0, 4,  2,  8, -2]];
    static immutable c = [[-42.0, -69, 23],
                          [ 35.0, -21, 69],
                          [ -7.0, -42,  3],
                          [ 77.0,  21, 29]];

    auto X = uninitSlice!double(3, 5);
    auto Xt = uninitSlice!double(5, 3);
    auto Y = uninitSlice!double(5, 4);
    auto Yt = uninitSlice!double(4, 5);
    auto YtXt1 = uninitSlice!double(4, 3);
    auto YtXt2 = uninitSlice!double(4, 3);
    auto result = uninitSlice!double(4, 3);

    X[] = a;
    Xt[] = at;
    Y[] = b;
    Yt[] = bt;
    result[] = c;

    Y.transposed.mtimesGeneralKernel(Xt, YtXt1);
    assert(YtXt1.equal(result));
    Yt.mtimesGeneralKernel(X.transposed, YtXt2);
    assert(YtXt2.equal(result));
}

/// General matrix multiplication (complex)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.complex: Complex;
    import mir.ndslice.allocation: uninitSlice;

    alias C = Complex!double;

    static immutable a = [[-5.0,  1,  7, 7, -4],
                          [-1.0, -5,  6, 3, -3],
                          [-5.0, -2, -3, 6,  0]];
    static immutable b = [[-5.0, -3,  3,  1],
                          [ 4.0,  3,  6,  4],
                          [-4.0, -2, -2,  2],
                          [-1.0,  9,  4,  8],
                          [ 9.0, 8,  3, -2]];
    static immutable c = [[-42.0,  35,  -7, 77],
                          [-69.0, -21, -42, 21],
                          [ 23.0,  69,   3, 29]];

    auto X = uninitSlice!C(3, 5);
    auto Y = uninitSlice!C(5, 4);
    auto XY = uninitSlice!C(3, 4);
    auto result = uninitSlice!C(3, 4);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesGeneralKernel(Y, XY);
    assert(XY.equal(result));
}

/// Matrix-vector multiplication, specialization for MxN times Nx1 & 1xM times MxN
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2, -3], [-2.0, 2, 3, 10], [0.0, 2, 1, 1]];
    static immutable b = [2.0, 3, 4, 5];
    static immutable c = [14.0, 64, 15];

    auto X = uninitSlice!double(3, 4);
    auto y = uninitSlice!double(4);
    auto Xy = uninitSlice!double(3);
    auto yXT = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    X.mtimesGeneralKernel(y, Xy);
    assert(Xy.equal(result));
    y.mtimesGeneralKernel(X.transposed, yXT);
    assert(yXT.equal(result));
}

/// Dot product
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [-5.0,  1,  7,  7, -4];
    static immutable b = [ 4.0, -4, -2, 10,  4];

    auto x = uninitSlice!double(5);
    auto y = uninitSlice!double(5);
    double z;

    x[] = a;
    y[] = b;

    x.mtimesGeneralKernel(y, z);
    assert(z == 16);
}

/++
Symmetric matrix multiplication.

Similar to `mtimesGeneralKernel`, but `a` is assumed to be a symmetric matrix.

Params:
    uplo = controls whether `a` is upper symmetric or lower symmetric
+/
template mtimesSymmetricKernel(Uplo uplo = Uplo.Upper)
{
    /+
    Params:
        a = m(rows) x m(cols) symmetric matrix
        b = m(rows) x n(cols) matrix
        c = m(rows) x n(cols) matrix
    +/
    void mtimesSymmetricKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 2, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` assumed to be a square matrix");
        assert(c.length!0 == a.length!0, "The first dimension of `c` must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of `c` must match the second dimension of `b`");
    }

    do
    {
        import mir.blas: symm;
        static import cblas;

        symm(cblas.Side.Left, uplo, cast(T)1, a, b, cast(T)0, c);
    }

    /+
    Params:
        a = m(rows) x m(cols) symmetric matrix
        b = m-dimensional vector
        c = m-dimensional vector
    +/
    @safe pure nothrow @nogc
    void mtimesSymmetricKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 1, kindB) b,
        Slice!(T*, 1, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
        assert(c.length!0 == a.length!0, "The length of `c` must match the first dimension of `a`");
    }
    do
    {
        import mir.blas: symv;

        symv(uplo, cast(T)1, a, b, cast(T)0, c);
    }
}

/// Symmetric matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable b = [[2.0, 3], [4.0, 3], [0.0, -5]];
    static immutable c = [[26.0, 14], [18.0, 6], [16.0, 10]];

    auto X = uninitSlice!double(3, 3);
    auto Y = uninitSlice!double(3, 2);
    auto XY = uninitSlice!double(3, 2);
    auto result = uninitSlice!double(3, 2);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesSymmetricKernel(Y, XY);
    assert(XY.equal(result));
}

/// Symmetric matrix multiplication, specialization for MxM times Mx1
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 28, 17];

    auto X = uninitSlice!double(3, 3);
    auto y = uninitSlice!double(3);
    auto Xy = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    X.mtimesSymmetricKernel(y, Xy);
    assert(Xy.equal(result));
}

/++
Symmetric matrix multiplication.

Similar to `mtimesGeneralKernel`, but `b` is assumed to be a symmetric matrix.

Params:
    uplo = controls whether `b` is upper symmetric or lower symmetric
+/
template mtimesSymmetricRightKernel(Uplo uplo = Uplo.Upper)
{
    /+
    Params:
        a = m(rows) x n(cols) matrix
        b = n(rows) x n(cols) symmetric matrix
        c = m(rows) x n(cols) matrix
    +/
    void mtimesSymmetricRightKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 2, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` assumed to be a square matrix");
        assert(c.length!0 == a.length!0, "The first dimension of `c` must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of `c` must match the second dimension of `b`");
    }

    do
    {
        import mir.blas: symm;
        static import cblas;

        symm(cblas.Side.Right, uplo, cast(T)1, b, a, cast(T)0, c);
    }

    /+
    Params:
        a = m-dimensional vector
        b = m(rows) x m(cols) symmetric matrix
        c = m-dimensional vector
    +/
    @safe pure nothrow @nogc
    void mtimesSymmetricRightKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 1, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 1, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!0 == b.length!0, "The length of `a` must match the second dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
        assert(c.length!0 == a.length!0, "The length of `c` must match the length of `a`");
    }
    do
    {
        import mir.ndslice.dynamic: transposed;
        import mir.ndslice.topology: universal;

        mtimesSymmetricKernel(b.universal.transposed, a, c);
    }
}

/// Symmetric matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[2.0, 4, 0], [3.0, 3, -5]];
    static immutable b = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable c = [[26.0, 18, 16], [14.0, 6, 10]];

    auto X = uninitSlice!double(2, 3);
    auto Y = uninitSlice!double(3, 3);
    auto XY = uninitSlice!double(2, 3);
    auto result = uninitSlice!double(2, 3);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesSymmetricRightKernel(Y, XY);
    assert(XY.equal(result));
}

/// Symmetric matrix multiplication, specialization for 1xM times MxM
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [2.0, 3, 4];
    static immutable b = [[3.0, 5, 2], [5.0, 2, 3], [2.0, 3, 1]];
    static immutable c = [29, 28, 17];

    auto x = uninitSlice!double(3);
    auto Y = uninitSlice!double(3, 3);
    auto xY = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    x[] = a;
    Y[] = b;
    result[] = c;

    x.mtimesSymmetricRightKernel(Y, xY);
    assert(xY.equal(result));
}

/++
Triangular matrix multiplication.

Similar to `mtimesGeneralKernel`, but `a` is assumed to be a triangular matrix.

Params:
    uplo = controls whether `a` is upper triangular or lower triangular
    diag = controls whether `a` is is non-unit triangular or unit triangular
+/
@safe pure nothrow @nogc
template mtimesTriangularKernel(Uplo uplo = Uplo.Upper, Diag diag = Diag.NonUnit)
{
    /+
    Params:
        a = m(rows) x m(cols) triangular matrix
        b = m(rows) x n(cols) matrix
        c = m(rows) x n(cols) matrix
    +/
    void mtimesTriangularKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 2, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
        assert(c.length!0 == a.length!0, "The first dimension of `c` must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of `c` must match the second dimension of `b`");
    }
    do
    {
        import mir.blas: Side, trmm;

        foreach (size_t i; 0 .. b.length!0) {
            foreach (size_t j; 0 .. b.length!1) {
                c[i, j] = b[i, j]; //note: purposefully filling with `b` here since trmm overwrites it
            }
        }
        trmm(Side.Left, uplo, diag, cast(T)1, a, c);
    }

    /++
    Params:
        a = m(rows) x m(cols) triangular matrix
        b = m-dimensional vector
        c = m-dimensional vector
    +/
    void mtimesTriangularKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 1, kindB) b,
        Slice!(T*, 1, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the length of `b`");
        assert(a.length!0 == a.length!1, "`a` must be a square matrix");
        assert(c.length == a.length, "The length of `c` must match the first dimension of `a`");
    }
    do
    {
        import mir.blas: trmv;

        foreach (size_t i; 0 .. b.length) {
            c[i] = b[i];
        }
        trmv(uplo, diag, a, c);
    }
}

/// Triangular matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;
    import mir.ndslice.dynamic: transposed;

    static immutable a = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable b = [[2.0, 3], [4.0, 3], [0.0, -5]];
    static immutable c = [[26.0, 14], [8.0, -9], [0.0, -5]];

    auto X = uninitSlice!double(3, 3);
    auto Y = uninitSlice!double(3, 2);
    auto XY = uninitSlice!double(3, 2);
    auto result = uninitSlice!double(3, 2);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesTriangularKernel(Y, XY);
    assert(XY.equal(result));
}

/// Triangular matrix multiplication, specialization for MxM times Mx1
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable b = [2.0, 3, 4];
    static immutable c = [29, 18, 4];

    auto X = uninitSlice!double(3, 3);
    auto y = uninitSlice!double(3);
    auto Xy = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    X[] = a;
    y[] = b;
    result[] = c;

    X.mtimesTriangularKernel(y, Xy);
    assert(Xy.equal(result));
}

/++
Triangular (right) matrix multiplication.

Similar to `mtimesGeneralKernel`, but `b` is assumed to be a triangular matrix.

Params:
    uplo = controls whether `b` is upper triangular or lower triangular
    diag = controls whether `b` is is non-unit triangular or unit triangular
+/
@safe pure nothrow @nogc
template mtimesTriangularRightKernel(Uplo uplo = Uplo.Upper, Diag diag = Diag.NonUnit)
{
    /+
    Params:
        a = m(rows) x n(cols) matrix
        b = n(rows) x n(cols) triangular matrix
        c = m(rows) x n(cols) matrix
    +/
    void mtimesTriangularRightKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 2, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 2, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!1 == b.length!0, "The second dimension of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
        assert(c.length!0 == a.length!0, "The first dimension of `c` must match the first dimension of `a`");
        assert(c.length!1 == b.length!1, "The second dimension of `c` must match the second dimension of `b`");
    }
    do
    {
        import mir.blas: Side, trmm;

        foreach (size_t i; 0 .. a.length!0) {
            foreach (size_t j; 0 .. a.length!1) {
                c[i, j] = a[i, j]; //note: purposefully filling with `a` here since trmm overwrites it
            }
        }
        trmm(Side.Right, uplo, diag, cast(T)1, b, c);
    }

    /++
    Params:
        a = m-dimensional vector
        b = m(rows) x m(cols) triangular matrix
        c = m-dimensional vector
    +/
    void mtimesTriangularRightKernel(T, SliceKind kindA, SliceKind kindB, SliceKind kindC)(
        Slice!(const(T)*, 1, kindA) a,
        Slice!(const(T)*, 2, kindB) b,
        Slice!(T*, 1, kindC) c
    )
        if (isFloatingPoint!T)
    in
    {
        assert(a.length!0 == b.length!0, "The length of `a` must match the first dimension of `b`");
        assert(b.length!0 == b.length!1, "`b` must be a square matrix");
        assert(c.length == b.length!1, "The length of `c` must match the second dimension of `b`");
    }
    do
    {
        import mir.ndslice.dynamic: transposed;
        import mir.ndslice.topology: universal;

        static if (uplo == Uplo.Upper) {
            .mtimesTriangularKernel!(Uplo.Lower, diag)(b.universal.transposed, a, c); //NOTE: switching Uplo.Lower to Uplo.Upper because of the transpose
        } else {
            .mtimesTriangularKernel!(Uplo.Upper, diag)(b.universal.transposed, a, c); //reverse
        }
    }
}

/// Triangular matrix multiplication
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[2.0, 4, 0], [3.0, 3, -5]];
    static immutable b = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable c = [[6.0, 18, 16], [9.0, 21, 10]];

    auto X = uninitSlice!double(2, 3);
    auto Y = uninitSlice!double(3, 3);
    auto XY = uninitSlice!double(2, 3);
    auto result = uninitSlice!double(2, 3);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesTriangularRightKernel(Y, XY);
    assert(XY.equal(result));
}

/// Triangular matrix multiplication, specialization for 1xM times MxM
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [2.0, 3, 4];
    static immutable b = [[3.0, 5, 2], [0.0, 2, 3], [0.0, 0, 1]];
    static immutable c = [6.0, 16, 17];

    auto x = uninitSlice!double(3);
    auto Y = uninitSlice!double(3, 3);
    auto xY = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    x[] = a;
    Y[] = b;
    result[] = c;

    x.mtimesTriangularRightKernel(Y, xY);
    assert(xY.equal(result));
}

// Triangular matrix multiplication (lower)
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [[2.0, 4, 0], [3.0, 3, -5]];
    static immutable b = [[3.0, 0, 0], [5.0, 2, 0], [2.0, 3, 1]];
    static immutable c = [[26.0, 8, 0], [14.0, -9, -5]];

    auto X = uninitSlice!double(2, 3);
    auto Y = uninitSlice!double(3, 3);
    auto XY = uninitSlice!double(2, 3);
    auto result = uninitSlice!double(2, 3);

    X[] = a;
    Y[] = b;
    result[] = c;

    X.mtimesTriangularRightKernel!(Uplo.Lower)(Y, XY);
    assert(XY.equal(result));
}

// Triangular matrix multiplication (lower), specialization for 1xM times MxM
version(mir_stat_test_blas)
@safe pure nothrow
unittest
{
    import mir.algorithm.iteration: equal;
    import mir.ndslice.allocation: uninitSlice;

    static immutable a = [2.0, 3, 4];
    static immutable b = [[3.0, 0, 0], [5.0, 2, 0], [2.0, 3, 1]];
    static immutable c = [29.0, 18, 4];

    auto x = uninitSlice!double(3);
    auto Y = uninitSlice!double(3, 3);
    auto xY = uninitSlice!double(3);
    auto result = uninitSlice!double(3);

    x[] = a;
    Y[] = b;
    result[] = c;

    x.mtimesTriangularRightKernel!(Uplo.Lower)(Y, xY);
    assert(xY.equal(result));
}

}
