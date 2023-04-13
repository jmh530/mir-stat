/++
This module contains algorithms for multivariate descriptive statistics.

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall

Copyright: 2023 Mir Stat Authors.

Macros:
SUBREF = $(REF_ALTTEXT $(TT $2), $2, mir, stat, $1)$(NBSP)
MATHREF = $(GREF_ALTTEXT mir-algorithm, $(TT $2), $2, mir, math, $1)$(NBSP)
MATHREF_ALT = $(GREF_ALTTEXT mir-algorithm, $(B $(TT $2)), $2, mir, math, $1)$(NBSP)
NDSLICEREF = $(GREF_ALTTEXT mir-algorithm, $(TT $2), $2, mir, ndslice, $1)$(NBSP)
T2=$(TR $(TDNW $(LREF $1)) $(TD $+))
T3=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $3))
T4=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $3) $(TD $4))

+/

// TODO: Add handling for N-dimensional slices that produce the covariance matrix

module mir.stat.descriptive.multivariate;

import mir.internal.utility: isFloatingPoint;
import mir.math.sum: Summation, Summator;
import std.traits: isMutable;

/++
Covariance algorithms.

See Also:
    $(WEB en.wikipedia.org/wiki/Algorithms_for_calculating_variance, Algorithms for calculating variance).
+/
enum CovarianceAlgo
{
    /++
    Performs Welford's online algorithm for updating covariance. While it only
    iterates each input once, it can be slower for smaller inputs. However, it
    is also more accurate. Can also `put` another CovarianceAccumulator of the
    same type, which uses the parallel algorithm from Chan et al.
    +/
    online,
    
    /++
    Calculates covariance using E(x*y) - E(x)*E(y) (alowing for adjustments for 
    population/sample variance). This algorithm can be numerically unstable.
    +/
    naive,

    /++
    Calculates covariance using a two-pass algorithm whereby the inputs are first 
    centered and then the sum of products is calculated from that. May be faster
    than `online` and generally more accurate than the `naive` algorithm.
    +/
    twoPass,

    /++
    Calculates covariance assuming the mean of the inputs is zero. 
    +/
    assumeZeroMean,
    onlineOld,
    online2
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isMutable!T && covarianceAlgo == CovarianceAlgo.naive)
{
    import mir.math.sum: Summator;
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, Slice, SliceKind;
    import mir.primitives: isInputRange, front, empty, popFront;

    ///
    private size_t _count;
    ///
    Summator!(T, summation) summatorLeft;
    ///
    Summator!(T, summation) summatorRight;
    ///
    Summator!(T, summation) sumOfProducts;

    ///
    this(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX &&
            isInputRange!RangeY)
    {
        import core.lifetime: move;
        this.put(x.move, y.move);
    }

    ///
    void put(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
        Slice!(IteratorX, 1, kindX) x,
        Slice!(IteratorY, 1, kindY) y
    )
    in
    {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
    }
    do
    {
        import mir.ndslice.topology: zip, map;

        _count += x.length;
        summatorLeft.put(x);
        summatorRight.put(y);
        sumOfProducts.put(x.zip(y).map!"a * b");
    }

    ///
    void put(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this.put(x.toSlice, y.toSlice);
    }

    ///
    void put(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && !isConvertibleToSlice!RangeX &&
            isInputRange!RangeY && !isConvertibleToSlice!RangeY)
    {
        do
        {
            assert(!(!x.empty && y.empty) && !(x.empty && !y.empty),
                   "x and y must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(x.front, y.front);
            x.popFront;
            y.popFront;
        } while(!x.empty || !y.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of x and y sould be caught by above assert
    }

    ///
    void put()(T x, T y)
    {
        _count++;
        summatorLeft.put(x);
        summatorRight.put(y);
        sumOfProducts.put(x * y);
    }

    ///
    void put()(CovarianceAccumulator!(T, covarianceAlgo, summation) v)
    {
        _count += v.count;
        summatorLeft.put(v.summatorLeft.sum);
        summatorRight.put(v.summatorRight.sum);
        sumOfProducts.put(v.sumOfProducts.sum);
    }

const:

    ///
    size_t count() @property
    {
        return _count;
    }

    ///
    F meanLeft(F = T)() const @property
    {
        return cast(F) summatorLeft.sum / count;
    }

    ///
    F meanRight(F = T)() const @property
    {
        return cast(F) summatorRight.sum / count;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) sumOfProducts.sum / (count + isPopulation - 1) - 
            (cast(F) summatorLeft.sum * cast(F) summatorRight.sum) * (F(1) / (count * (count + isPopulation - 1)));
    }
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == 82.25 / 12 - (29.25 * 36) / (12 * 12);
    v.covariance(false).shouldApprox == 82.25 / 11 - (29.25 * 36) / (12 * 12) * (12.0 / 11);

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == 94.25 / 13 - (33.25 * 39) / (13 * 13);
    v.covariance(false).shouldApprox == 94.25 / 12 - (33.25 * 39) / (13 * 13) * (13.0 / 12);
}

// Check dynamic array
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25];

    CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == 82.25 / 12 - (29.25 * 36) / (12 * 12);
    v.covariance(false).shouldApprox == 82.25 / 11 - (29.25 * 36) / (12 * 12) * (12.0 / 11);

    v.meanLeft.shouldApprox == 2.4375;
    v.meanRight.shouldApprox == 3;

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == 94.25 / 13 - (33.25 * 39) / (13 * 13);
    v.covariance(false).shouldApprox == 94.25 / 12 - (33.25 * 39) / (13 * 13) * (13.0 / 12);
}

// rcslice test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto x = mininitRcslice!double(12);
    auto y = mininitRcslice!double(12);
    x[] = a;
    y[] = b;
    auto v = CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Check adding CovarianceAccumultors
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: sum, Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x1 = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25].sliced;
    auto y1 = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5].sliced;
    auto x2 = [  2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y2 = [ 9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive) v1;
    v1.put(x1, y1);
    CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive) v2;
    v2.put(x2, y2);
    v1.put(v2);

    v1.covariance(true).shouldApprox == -5.5 / 12;
    v1.covariance(false).shouldApprox == -5.5 / 11;
}

// Test input range
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;
    import std.range: iota;

    auto x = iota(0, 5);
    auto y = iota(-3, 2);
    CovarianceAccumulator!(double, CovarianceAlgo.naive, Summation.naive) v;
    v.put(x, y);
    v.covariance(true).should == 10.0 / 5;
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isFloatingPoint!T && isMutable!T && covarianceAlgo == CovarianceAlgo.online2)
{
    import mir.math.sum: Summator;
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, Slice, SliceKind;
    import mir.primitives: isInputRange, front, empty, popFront;

    ///
    private size_t _count;
    ///
    private T _meanLeft = 0;
    ///
    private T _meanRight = 0;
    ///
    Summator!(T, summation) centeredSumOfProducts;

    ///
    this(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && isInputRange!RangeY)
    {
        import core.lifetime: move;
        _meanLeft = x.front;
        _meanRight = y.front;
        _count++;
        x.popFront;
        y.popFront;
        this.put(x.move, y.move);
    }

    ///
    this()(T x, T y)
    {
        _meanLeft = x;
        _meanRight = y;
        _count++;
    }

    ///
    void put(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
        Slice!(IteratorX, 1, kindX) x,
        Slice!(IteratorY, 1, kindY) y
    )
    in
    {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
    }
    do
    {
        import mir.ndslice.topology: zip;

        if (count == 0) {
            _meanLeft = x[0];
            _meanRight = y[0];
            _count = 1;
            if (x.length > 1) {
                foreach(e; x[1 .. $].zip(y[1 .. $])) {
                    this.put(e[0], e[1]);
                }
            }
        } else {
            foreach(e; x.zip(y)) {
                this.put(e[0], e[1]);
            }
        }
    }

    ///
    void put(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this.put(x.toSlice, y.toSlice);
    }

    ///
    void put(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && !isConvertibleToSlice!RangeX &&
            isInputRange!RangeY && !isConvertibleToSlice!RangeY)
    {
        do
        {
            assert(!(!x.empty && y.empty) && !(x.empty && !y.empty),
                   "x and y must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(x.front, y.front);
            x.popFront;
            y.popFront;
        } while(!x.empty || !y.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of x and y sould be caught by above assert
    }

    ///
    void put()(T x, T y)
    {
        T delta = x - meanLeft;
        _count++;
        _meanLeft += (x - _meanLeft) / _count;
        _meanRight += (y - _meanRight) / _count;
        centeredSumOfProducts.put(delta * (y - meanRight));
    }

    ///
    void put()(CovarianceAccumulator!(T, covarianceAlgo, summation) v)
    {
        size_t oldCount = count;
        T deltaLeft = v.meanLeft - meanLeft;
        T deltaRight = v.meanRight - meanRight;
        _count += v.count;
        _meanLeft = (_meanLeft * oldCount + v.count * v.meanLeft) / _count;
        _meanRight = (_meanRight * oldCount + v.count * v.meanRight) / _count;
        centeredSumOfProducts.put(v.centeredSumOfProducts.sum + deltaLeft * deltaRight * v.count * oldCount / count);
    }

const:

    ///
    size_t count() @property
    {
        return _count;
    }

    ///
    F meanLeft(F = T)() const @property
    {
        return _meanLeft;
    }
    ///
    F meanRight(F = T)() const @property
    {
        return _meanRight;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) centeredSumOfProducts.sum / (count + isPopulation - 1);
    }
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isFloatingPoint!T && isMutable!T && covarianceAlgo == CovarianceAlgo.online)
{
    import mir.math.sum: Summator;
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, Slice, SliceKind;
    import mir.primitives: isInputRange, front, empty, popFront;

    private size_t _count;
    private T _meanLeft = 0;
    private T _meanRight = 0;
    private T _covariance = 0;

    ///
    this(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && isInputRange!RangeY)
    {
        import core.lifetime: move;
        _meanLeft = x.front;
        _meanRight = y.front;
        _count++;
        x.popFront;
        y.popFront;
        this.put(x.move, y.move);
    }

    ///
    this()(T x, T y)
    {
        _meanLeft = x;
        _meanRight = y;
        _count++;
    }

    ///
    void put(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
        Slice!(IteratorX, 1, kindX) x,
        Slice!(IteratorY, 1, kindY) y
    )
    in
    {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
    }
    do
    {
        import mir.ndslice.topology: zip;

        if (count == 0) {
            _meanLeft = x[0];
            _meanRight = y[0];
            _count = 1;
            if (x.length > 1) {
                foreach(e; x[1 .. $].zip(y[1 .. $])) {
                    this.put(e[0], e[1]);
                }
            }
        } else {
            foreach(e; x.zip(y)) {
                this.put(e[0], e[1]);
            }
        }
    }

    ///
    void put(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this.put(x.toSlice, y.toSlice);
    }

    ///
    void put(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && !isConvertibleToSlice!RangeX &&
            isInputRange!RangeY && !isConvertibleToSlice!RangeY)
    {
        do
        {
            assert(!(!x.empty && y.empty) && !(x.empty && !y.empty),
                   "x and y must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(x.front, y.front);
            x.popFront;
            y.popFront;
        } while(!x.empty || !y.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of x and y sould be caught by above assert
    }

    ///
    void put()(T x, T y)
    {
        size_t oldCount = _count;
        T delta = x - _meanLeft;
        _count++;
        _meanLeft += (x - _meanLeft) / _count;
        _meanRight += (y - _meanRight) / _count;
        _covariance += (delta * (y - _meanRight) - _covariance) / _count;
    }

    ///
    void put()(CovarianceAccumulator!(T, covarianceAlgo, summation) v)
    {
        size_t oldCount = count;
        T deltaLeft = v.meanLeft - meanLeft;
        T deltaRight = v.meanRight - meanRight;
        _count += v.count;
        _meanLeft = (_meanLeft * oldCount + v.count * v.meanLeft) / _count;
        _meanRight = (_meanRight * oldCount + v.count * v.meanRight) / _count;
        _covariance = _covariance * oldCount / count + v._covariance * v.count / count + deltaLeft * deltaRight * v.count * oldCount / (count * count);
    }

const:

    ///
    size_t count() @property
    {
        return _count;
    }

    ///
    F meanLeft(F = T)() const @property
    {
        return _meanLeft;
    }

    ///
    F meanRight(F = T)() const @property
    {
        return _meanRight;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) _covariance * count / (count + isPopulation - 1);
    }
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isFloatingPoint!T && isMutable!T && covarianceAlgo == CovarianceAlgo.onlineOld)
{
    import mir.math.sum: Summator;
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, Slice, SliceKind;
    import mir.primitives: isInputRange, front, empty, popFront;

    ///
    private size_t _count;
    ///
    Summator!(T, summation) summatorLeft;
    ///
    Summator!(T, summation) summatorRight;
    ///
    Summator!(T, summation) centeredSumOfProducts;

    ///
    this(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && isInputRange!RangeY)
    {
        import core.lifetime: move;
        this.put(x.move, y.move);
    }

    ///
    this()(T x, T y)
    {
        this.put(x, y);
    }

    ///
    void put(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
        Slice!(IteratorX, 1, kindX) x,
        Slice!(IteratorY, 1, kindY) y
    )
    in
    {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
    }
    do
    {
        import mir.ndslice.topology: zip;

        foreach(e; x.zip(y)) {
            this.put(e[0], e[1]);
        }
    }

    ///
    void put(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this.put(x.toSlice, y.toSlice);
    }

    ///
    void put(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && !isConvertibleToSlice!RangeX &&
            isInputRange!RangeY && !isConvertibleToSlice!RangeY)
    {
        do
        {
            assert(!(!x.empty && y.empty) && !(x.empty && !y.empty),
                   "x and y must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(x.front, y.front);
            x.popFront;
            y.popFront;
        } while(!x.empty || !y.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of x and y sould be caught by above assert
    }

    ///
    void put()(T x, T y)
    {
        T delta = x;
        if (count > 0) {
            delta -= meanLeft;
        }
        _count++;
        summatorLeft.put(x);
        summatorRight.put(y);
        centeredSumOfProducts.put(delta * (y - meanRight));
    }

    ///
    void put()(CovarianceAccumulator!(T, covarianceAlgo, summation) v)
    {
        size_t oldCount = count;
        T deltaLeft = v.meanLeft;
        T deltaRight = v.meanRight;
        if (oldCount > 0) {
            deltaLeft -= meanLeft;
            deltaRight -= meanRight;
        }
        _count += v.count;
        summatorLeft.put!T(v.summatorLeft.sum);
        summatorRight.put(v.summatorRight.sum);
        centeredSumOfProducts.put(v.centeredSumOfProducts.sum + deltaLeft * deltaRight * v.count * oldCount / count);
    }

const:

    ///
    size_t count() @property
    {
        return _count;
    }

    ///
    F meanLeft(F = T)() const @property
    {
        return summatorLeft.sum / count;
    }
    ///
    F meanRight(F = T)() const @property
    {
        return summatorRight.sum / count;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) centeredSumOfProducts.sum / (count + isPopulation - 1);
    }
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == -5.5 / 13;
    v.covariance(false).shouldApprox == -5.5 / 12;
}

// Check dynamic array
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25];

    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;

    v.meanLeft.shouldApprox == 2.4375;
    v.meanRight.shouldApprox == 3;

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == -5.5 / 13;
    v.covariance(false).shouldApprox == -5.5 / 12;
}

// rcslice test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto x = mininitRcslice!double(12);
    auto y = mininitRcslice!double(12);
    x[] = a;
    y[] = b;
    auto v = CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Check adding CovarianceAccumultors
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: sum, Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x1 = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25].sliced;
    auto y1 = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5].sliced;
    auto x2 = [  2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y2 = [ 9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v1;
    v1.put(x1, y1);
    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v2;
    v2.put(x2, y2);
    v1.put(v2);

    v1.covariance(true).shouldApprox == -5.5 / 12;
    v1.covariance(false).shouldApprox == -5.5 / 11;
}

// Initializing with one point
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive)(4.0, 3.0);
    v.covariance(true).should == 0;
}

// Test input range
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;
    import std.range: iota;

    auto x = iota(0, 5);
    auto y = iota(-3, 2);
    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v;
    v.put(x, y);
    v.covariance(true).should == 10 / 5;
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isMutable!T && covarianceAlgo == CovarianceAlgo.twoPass)
{
    import mir.functional: naryFun;
    import mir.math.stat: MeanAccumulator;
    import mir.math.sum: Summator;
    import mir.ndslice.slice: isConvertibleToSlice, isSlice, Slice, SliceKind;

    ///
    MeanAccumulator!(T, summation) meanAccumulatorLeft;
    ///
    MeanAccumulator!(T, summation) meanAccumulatorRight;
    ///
    Summator!(T, summation) centeredSumOfProducts;

    ///
    this(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
         Slice!(IteratorX, 1, kindX) x, Slice!(IteratorY, 1, kindY) y)
     in
     {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
     }
     do
    {
        import core.lifetime: move;
        import mir.ndslice.internal: LeftOp;
        import mir.ndslice.topology: map, vmap, zip;

        meanAccumulatorLeft.put(x.lightScope);
        meanAccumulatorRight.put(y.lightScope);
        centeredSumOfProducts.put(x.move.vmap(LeftOp!("-", T)(meanAccumulatorLeft.mean)).zip(y.move.vmap(LeftOp!("-", T)(meanAccumulatorRight.mean))).map!(naryFun!"a * b"));
    }

    ///
    this(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this(x.toSlice, y.toSlice);
    }

    ///
    this()(T x, T y)
    {
        meanAccumulatorLeft.put(x);
        meanAccumulatorRight.put(y);
        centeredSumOfProducts.put(cast(T) 0);
    }

const:

    ///
    size_t count() @property
    {
        assert(meanAccumulatorLeft.count == meanAccumulatorRight.count);
        return meanAccumulatorLeft.count;
    }

    ///
    F meanLeft(F = T)() const @property
    {
        return meanAccumulatorLeft.mean!F;
    }
    ///
    F meanRight(F = T)() const @property
    {
        return meanAccumulatorRight.mean!F;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) centeredSumOfProducts.sum / (count + isPopulation - 1);
    }
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Check dynamic array
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25];

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;

    v.meanLeft.shouldApprox == 2.4375;
    v.meanRight.shouldApprox ==3;
}

// rcslice test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto x = mininitRcslice!double(12);
    auto y = mininitRcslice!double(12);
    x[] = a;
    y[] = b;
    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Check Vmap
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;
    auto x = a + 1;
    auto y = b - 1;

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Initializing with one point
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(4.0, 3.0);
    v.centeredSumOfProducts.sum.should == 0;
}

// withAsSlice test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array: RCArray;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];

    auto x = RCArray!double(12);
    foreach(i, ref e; x)
        e = a[i];
    auto y = RCArray!double(12);
    foreach(i, ref e; y)
        e = b[i];

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.twoPass, Summation.naive)(x, y);
    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

///
struct CovarianceAccumulator(T, CovarianceAlgo covarianceAlgo, Summation summation)
    if (isMutable!T && covarianceAlgo == CovarianceAlgo.assumeZeroMean)
{
    import mir.math.sum: Summator;
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice, isConvertibleToSlice, isSlice;
    import mir.primitives: isInputRange, front, empty, popFront;

    private size_t _count;

    ///
    Summator!(T, summation) centeredSumOfProducts;

    ///
    this(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && isInputRange!RangeY)
    {
        this.put(x, y);
    }

    ///
    this()(T x, T y)
    {
        this.put(x, y);
    }

    ///
    void put(IteratorX, IteratorY, SliceKind kindX, SliceKind kindY)(
        Slice!(IteratorX, 1, kindX) x,
        Slice!(IteratorY, 1, kindY) y
    )
    in
    {
        assert(x.length == y.length,
               "CovarianceAcumulator.put: both vectors must have the same length");
    }
    do
    {
        import mir.ndslice.topology: zip, map;

        _count += x.length;
        centeredSumOfProducts.put(x.zip(y).map!"a * b");
    }

    ///
    void put(SliceLikeX, SliceLikeY)(SliceLikeX x, SliceLikeY y)
        if (isConvertibleToSlice!SliceLikeX && !isSlice!SliceLikeX &&
            isConvertibleToSlice!SliceLikeY && !isSlice!SliceLikeY)
    {
        import mir.ndslice.slice: toSlice;
        this.put(x.toSlice, y.toSlice);
    }

    ///
    void put(RangeX, RangeY)(RangeX x, RangeY y)
        if (isInputRange!RangeX && !isConvertibleToSlice!RangeX &&
            isInputRange!RangeY && !isConvertibleToSlice!RangeY)
    {
        do
        {
            assert(!(!x.empty && y.empty) && !(x.empty && !y.empty),
                   "x and y must both be empty at the same time, one cannot be empty while the other has remaining items");
            this.put(x.front, y.front);
            x.popFront;
            y.popFront;
        } while(!x.empty || !y.empty); // Using an || instead of && so that the loop does not end early. mis-matched lengths of x and y sould be caught by above assert
    }

    ///
    void put()(T x, T y)
    {
        _count++;
        centeredSumOfProducts.put(x * y);
    }

    ///
    void put()(CovarianceAccumulator!(T, covarianceAlgo, summation) v)
    {
        _count += v.count;
        centeredSumOfProducts.put(v.centeredSumOfProducts.sum);
    }

const:

    ///
    size_t count() @property
    {
        return _count;
    }
    
    ///
    F meanLeft(F = T)() const @property
    {
        return 0;
    }

    ///
    F meanRight(F = T)() const @property
    {
        return 0;
    }

    ///
    F covariance(F = T)(bool isPopulation) @property
    {
        return cast(F) centeredSumOfProducts.sum / (count + isPopulation - 1);
    }
}

///
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.stat: center;
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;
    auto x = a.center;
    auto y = b.center;

    CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == 6.5 / 13;
    v.covariance(false).shouldApprox == 6.5 / 12;
}

// Check dynamic array
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.stat: center, mean;
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: should, shouldApprox;

    auto a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    auto b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto aMean = a.mean;
    auto bMean = b.mean;
    auto x = a.dup;
    auto y = b.dup;
    for (size_t i; i < a.length; i++) {
        x[i] -= aMean;
        y[i] -= bMean;
    }

    CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive) v;
    v.put(x, y);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;

    v.put(4.0, 3.0);
    v.covariance(true).shouldApprox == 6.5 / 13;
    v.covariance(false).shouldApprox == 6.5 / 12;
    
    v.meanLeft.should == 0;
    v.meanRight.should == 0;
}

// rcslice test
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.stat: center;
    import mir.math.sum: Summation;
    import mir.ndslice.allocation: mininitRcslice;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto x = mininitRcslice!double(12);
    auto y = mininitRcslice!double(12);
    x[] = a;
    y[] = b;
    auto v = CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive)(x.center, y.center);

    v.covariance(true).shouldApprox == -5.5 / 12;
    v.covariance(false).shouldApprox == -5.5 / 11;
}

// Check adding CovarianceAccumultors
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: sum, Summation;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto a1 = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25].sliced;
    auto b1 = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5].sliced;
    auto a2 = [  2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto b2 = [ 9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;
    auto meanA = (a1.sum + a2.sum) / 12;
    auto meanB = (b1.sum + b2.sum) / 12;
    auto x1 = a1 - meanA;
    auto y1 = b1 - meanB;
    auto x2 = a2 - meanA;
    auto y2 = b2 - meanB;

    CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive) v1;
    v1.put(x1, y1);
    CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive) v2;
    v2.put(x2, y2);
    v1.put(v2);

    v1.covariance(true).shouldApprox == -5.5 / 12;
    v1.covariance(false).shouldApprox == -5.5 / 11;
}

// Initializing with one point
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;

    auto v = CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive)(4.0, 3.0);
    v.centeredSumOfProducts.sum.should == 12;
}

// Test input range
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.sum: Summation;
    import mir.test: should;
    import std.range: iota;

    auto x = iota(0, 5);
    auto y = iota(-3, 2);
    auto v = CovarianceAccumulator!(double, CovarianceAlgo.assumeZeroMean, Summation.naive)(x, y);
    v.centeredSumOfProducts.sum.should == 0;
}

/++
Calculates the covariance of the inputs.

If `x` and `y` are both slices or convertible to slices, then they must be 
one-dimensional.

By default, if `F` is not floating point type, then the result will have a
`double` type if `F` is implicitly convertible to a floating point type.

Params:
    F = controls type of output
    covarianceAlgo = algorithm for calculating covariance (default: CovarianceAlgo.online)
    summation = algorithm for calculating sums (default: Summation.appropriate)
Returns:
    The covariance of the inputs
+/
template cov(F, 
             CovarianceAlgo covarianceAlgo = CovarianceAlgo.online, 
             Summation summation = Summation.appropriate)
    if (isFloatingPoint!F)
{
    import mir.math.common: fmamath;
    import mir.math.stat: meanType;
    import mir.math.sum: ResolveSummationType;
    import std.traits: CommonType;
    /++
    Params:
        x = range, must be finite iterable
        y = range, must be finite iterable
        isPopulation = true if population covariance, false if sample covariance (default)
    +/
    @fmamath F cov(RangeX, RangeY)(RangeX x, RangeY y, bool isPopulation = false)
    {
        import core.lifetime: move;

        auto covarianceAccumulator = CovarianceAccumulator!(F, covarianceAlgo, ResolveSummationType!(summation, RangeX, F))(x.move, y.move);
        return covarianceAccumulator.covariance(isPopulation);
    }
}

/// ditto
template cov(
    CovarianceAlgo covarianceAlgo = CovarianceAlgo.online, 
    Summation summation = Summation.appropriate)
{
    import mir.math.common: fmamath;
    import mir.math.stat: meanType;
    import std.traits: CommonType;
    /++
    Params:
        x = range, must be finite iterable
        y = range, must be finite iterable
        isPopulation = true if population covariance, false if sample covariance (default)
    +/
    @fmamath CommonType!(meanType!RangeX, meanType!RangeY) cov(RangeX, RangeY)(RangeX x, RangeY y, bool isPopulation = false)
    {
        import core.lifetime: move;

        alias F = typeof(return);
        return .cov!(F, covarianceAlgo, summation)(x.move, y.move, isPopulation);
    }
}

/// ditto
template cov(F, string covarianceAlgo, string summation = "appropriate")
{
    mixin("alias cov = .cov!(F, CovarianceAlgo." ~ covarianceAlgo ~ ", Summation." ~ summation ~ ");");
}

/// ditto
template cov(string covarianceAlgo, string summation = "appropriate")
{
    mixin("alias cov = .cov!(CovarianceAlgo." ~ covarianceAlgo ~ ", Summation." ~ summation ~ ");");
}

/// Covariance of vectors
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                2.0,   7.5,   5.0,  1.0,  1.5,  0.0].sliced;
    auto y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    x.cov(y, true).shouldApprox == -5.5 / 12;
    x.cov(y).shouldApprox == -5.5 / 11;
}

/// Can also set algorithm type
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.math.common: approxEqual;
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto a = [0.0, 1.0, 1.5, 2.0, 3.5, 4.25,
              2.0, 7.5, 5.0, 1.0, 1.5, 0.0].sliced;
    auto b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
               9.25, -0.75,   2.5, 1.25,   -1, 2.25].sliced;

    auto x = a + 10.0 ^^ 9;
    auto y = b + 10.0 ^^ 9;

    x.cov(y).shouldApprox == -5.5 / 11;

    // The naive algorithm is numerically unstable in this case
    assert(!x.cov!"naive"(y).approxEqual(-5.5 / 11));

    // But the two-pass algorithm provides a consistent answer
    x.cov!"twoPass"(y).shouldApprox == -5.5 / 11;

    // And the assumeZeroMean algorithm is way off
    assert(!x.cov!"assumeZeroMean"(y).approxEqual(-5.5 / 11));
}

/// Can also set algorithm type
version(mir_stat_test)
@safe
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;
    import std.array: array;
    import std.random: uniform;
    import std.range : generate, takeExactly;

    auto a = generate!(() => uniform(0.0, 1)).takeExactly(1_000_000).array.sliced;
    auto b = generate!(() => uniform(0.0, 1)).takeExactly(1_000_000).array.sliced;
    for (size_t i; i < a.length; i++) {
        b[i] += a[i];
    }
    auto x = a + 10.0 ^^ 12;
    auto y = b + 10.0 ^^ 12;
    x.cov(y).shouldApprox(0.0001) == a.cov(b);
    x.cov!"twoPass"(y).shouldApprox(0.0001) == a.cov!"twoPass"(b);
    /*
    import std.stdio: writeln;
    writeln(x.cov!"twoPass"(y));
    CovarianceAccumulator!(double, CovarianceAlgo.online, Summation.naive) v1;
    CovarianceAccumulator!(double, CovarianceAlgo.onlineOld, Summation.naive) v2;
    v1.put(x, y);
    v2.put(x, y);
    writeln(v1.meanLeft);
    writeln(v2.meanLeft);
    writeln(v1.meanRight);
    writeln(v2.meanRight);
    writeln(v1.covariance(false));
    writeln(v2.covariance(false));
    writeln(x.cov(y));
    writeln(a.cov(b));
    writeln(a.cov!"twoPass"(b));
    */
}

/// Can also set algorithm or output type
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.topology: repeat;
    import mir.test: shouldApprox;

    //Set population covariance, covariance algorithm, sum algorithm or output type

    auto a = [1.0, 1e100, 1, -1e100].sliced;
    auto b = [1.0e100, 1, 1, -1e100].sliced;
    auto x = a * 10_000;
    auto y = b * 10_000;

    /++
    Due to Floating Point precision, when centering `x`, subtracting the mean 
    from the second and fourth numbers has no effect (for `y` the same is true
    for the first and fourth). Further, after centering and multiplying `x` and
    `y`, the third numbers in the slice has precision too low to be included in
    the centered sum of the products. 
    +/
    x.cov(y).shouldApprox == 1.0e208 / 3;
    x.cov(y, true).shouldApprox == 1.0e208 / 4;

    x.cov!("online")(y).shouldApprox == 1.0e208 / 3;
    x.cov!("online", "kbn")(y).shouldApprox == 1.0e208 / 3;
    x.cov!("online", "kb2")(y).shouldApprox == 1.0e208 / 3;
    x.cov!("online", "precise")(y).shouldApprox == 1.0e208 / 3;
    x.cov!(double, "online", "precise")(y).shouldApprox == 1.0e208 / 3;

    auto z = uint.max.repeat(3);
    z.cov!float(z).shouldApprox == 0.0;
    static assert(is(typeof(z.cov!float(z)) == float));
}

/++
For integral slices, pass output type as template parameter to ensure output
type is correct.
+/
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.test: shouldApprox;

    auto x = [0, 1, 1, 2, 4, 4,
              2, 7, 5, 1, 2, 0].sliced;
    auto y = [6, 3, 7, 1, 1, 1,
              9, 5, 3, 1, 3, 7].sliced;

    x.cov(y).shouldApprox == -18.583333 / 11;
    static assert(is(typeof(x.cov(y)) == double));

    x.cov!float(y).shouldApprox == -18.583333 / 11;
    static assert(is(typeof(x.cov!float(y)) == float));
}

// make sure works with dynamic array
version(mir_stat_test)
@safe pure nothrow
unittest
{
    import mir.test: shouldApprox;

    double[] x = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                    2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    double[] y = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                   9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    x.cov(y).shouldApprox == -5.5 / 11;
}

/// Works with @nogc
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: mininitRcslice;
    import mir.test: shouldApprox;

    static immutable a = [  0.0,   1.0,   1.5,  2.0,  3.5, 4.25,
                            2.0,   7.5,   5.0,  1.0,  1.5,  0.0];
    static immutable b = [-0.75,   6.0, -0.25, 8.25, 5.75,  3.5,
                           9.25, -0.75,   2.5, 1.25,   -1, 2.25];
    auto x = mininitRcslice!double(12);
    auto y = mininitRcslice!double(12);
    x[] = a;
    y[] = b;

    x.cov(y, true).shouldApprox == -5.5 / 12;
    x.cov(y).shouldApprox == -5.5 / 11;
}

version(mir_stat_test_cov_performance)
unittest
{
    import mir.math.sum: Summation;
    import mir.ndslice.slice: sliced;
    import std.random: uniform;
    import std.array : array;
    import std.datetime.stopwatch;
    import std.range : generate, takeExactly;
    import std.stdio: writeln;

    size_t n = 10000;
    double aLow = -4.0;
    double aHigh = 3;
    double bLow = -3;
    double bHigh = 2;
    double[] a1 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b1 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x1 = a1.sliced;
    auto y1 = b1.sliced;
    double[] a2 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b2 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x2 = a2.sliced;
    auto y2 = b2.sliced;
    double[] a3 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b3 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x3 = a3.sliced;
    auto y3 = b3.sliced;
    double[] a4 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b4 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x4 = a4.sliced;
    auto y4 = b4.sliced;
    double[] a5 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b5 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x5 = a5.sliced;
    auto y5 = b5.sliced;
    double[] a6 = generate!(() => uniform(aLow, aHigh)).takeExactly(n).array;
    double[] b6 = generate!(() => uniform(bLow, bHigh)).takeExactly(n).array;
    auto x6 = a6.sliced;
    auto y6 = b6.sliced;

    void f1()
    {
        x1.cov!"online"(y1);
    }
    void f2()
    {
        x2.cov!"naive"(y2);
    }
    void f3()
    {
        x3.cov!"twoPass"(y3);
    }
    void f4()
    {
        x4.cov!"assumeZeroMean"(y4);
    }
    void f5()
    {
        x5.cov!"onlineOld"(y5);
    }
    void f6()
    {
        x6.cov!"online2"(y6);
    }

    auto r = benchmark!(f1, f2, f3, f4, f5, f6)(10_000);
    for (size_t i; i < r.length; i++) {
        writeln(r[i]);
    }
}

//TODO: compile with dub test --build=mir_stat_test_cov_performance --compiler=ldc2