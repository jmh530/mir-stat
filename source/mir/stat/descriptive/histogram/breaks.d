/++
This module contains algorithms for calculating breaks for histograms.

License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

Authors: John Michael Hall, Ilya Yaroshenko

Copyright: 2020 Mir Stat Authors.

Macros:
SUBREF = $(REF_ALTTEXT $(TT $2), $2, mir, stat, $1)$(NBSP)
MATHREF = $(REF_ALTTEXT $(TT $2), $2, mir, math, $1)$(NBSP)
NDSLICEREF = $(REF_ALTTEXT $(TT $2), $2, mir, ndslice, $1)$(NBSP)
T2=$(TR $(TDNW $(LREF $1)) $(TD $+))
T4=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $3) $(TD $4))
+/

module mir.stat.descriptive.histogram.breaks;

import mir.math.stat: VarianceAlgo;
import mir.primitives: hasShape;
import mir.stat.descriptive.univariate: QuantileAlgo;
import std.traits: isIntegral;

/++
Computes the number of breaks for a histogram using the Sturges' formula.

Calculates the ceiling of the base 2 logarithm of the number of elements of the
input and then adds 1.

Sturges' formula implicitly assumes an approximately normal distribution.

Params:
    CountType = the type that is used to count in histogram bins
Returns:
    The number of breaks

See_also: 
    $(LREF scott),
    $(LREF freedmanDiaconis),
    $(WEB en.wikipedia.org/wiki/Histogram, Histogram)
+/
template sturges(CountType)
    if(isIntegral!CountType)
{
    import mir.primitives: hasShape;

    /++
    Params:
        x = input
    +/
    CountType sturges(T)(T x)
        if (hasShape!T)
    {
        size_t n = x.elementCount;

        return sturges(n);
    }

    /++
    Params:
        n = count
    +/
    CountType sturges(size_t n)
    {
        assert(n > 0, "sturges: elementCount must be greater than zero");

        import mir.math.common: ceil, log2;
        import mir.primitives: elementCount;

        return cast(CountType) (ceil(log2(cast(double) n)) + 1);
    }
}

/++
Params:
    x = input
+/
size_t sturges(T)(T x)
    if (hasShape!T)
{
    return .sturges!size_t(x);
}

///
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7].sliced;

    auto k = x.sturges!size_t;
    assert(k == 4);
    static assert(is(typeof(k) == size_t));

    auto l = x.sturges;
    assert(l == 4);
}

/++
Calculates the number of breaks for a histogram given a step size `h`.

The formula calculates the ceiling of the result from first calculating the
difference between the maximum value of the input and the minimum value of the
input and then dividing that by the step size.

Params:
    CountType = the type that is used to count in histogram bins
Returns:
    The number of breaks

See_also:
    $(LREF scott),
    $(LREF struges),
    $(LREF freedmanDiaconis),
    $(WEB en.wikipedia.org/wiki/Histogram, Histogram)
+/
template binsFromWidth(CountType)
    if(isIntegral!CountType)
{
    import mir.internal.utility: isFloatingPoint;
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice;
    import std.traits: Unqual;

    /++
    Params:
        slice = slice
    +/
    CountType binsFromWidth(Iterator, size_t N, SliceKind kind, H)(
        Slice!(Iterator, N, kind) slice, H h)
    {
        import mir.algorithm.iteration: minmaxIndex;
        import mir.math.common: ceil;

        auto indexes = slice.minmaxIndex;

        return cast(CountType) ceil((slice[indexes[1]] - slice[indexes[0]]) / h);
    }

    /++
    Params:
        array = array
    +/
    CountType binsFromWidth(T, H)(T[] array, H h)
        if (isFloatingPoint!(Unqual!T))
    {
        import mir.ndslice.slice: sliced;

        return binsFromWidth(array.sliced, h);
    }

    /++
    Params:
        withAsSlice = withAsSlice
    +/
    CountType binsFromWidth(T, H)(T withAsSlice, H h)
        if (hasAsSlice!T)
    {
        return binsFromWidth(withAsSlice.asSlice, h);
    }
}

/// binsFromWidth example
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11].sliced;

    auto k = x.binsFromWidth!size_t(5.496295);
    assert(k == 3);
    static assert(is(typeof(k) == size_t));
}

// withAsSlice test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.rc.array: RCArray;
    import mir.algorithm.iteration: minmaxPos, minPos, maxPos, minmaxIndex, minIndex, maxIndex;

    static immutable a = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

    auto x = RCArray!double(12);
    foreach(i, ref e; x)
        e = a[i];
    auto y = x.asSlice;

    auto k = x.binsFromWidth!size_t(5.496295);
    assert(k == 3);
}

// dynamic array test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

    auto k = x.binsFromWidth!size_t(5.496295);
    assert(k == 3);
}

/++
Computes the number of breaks for a histogram using the Scott's normal
reference rule.

Calculates the break width using 3.49 times the sample standard deviation
divided by the cubed root of the number of elements of the input. The number of
breaks is then calculated using `binsFromWidth`.

It is considered optimal for normally distributed random samples.

Params:
    CountType = the type that is used to count in histogram bins
    varianceAlgo = Algorithm used to calculate variance
Returns:
    The number of breaks

See_also: 
    $(LREF struges),
    $(LREF freedmanDiaconis),
    $(LREF binsFromWidth),
    $(WEB en.wikipedia.org/wiki/Histogram, Histogram)
+/
template scott(CountType, VarianceAlgo varianceAlgo = VarianceAlgo.online)
    if(isIntegral!CountType)
{
    import mir.internal.utility: isFloatingPoint;
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice;
    import std.traits: isIterable, Unqual;

    /++
    Params:
        r = range
    +/
    CountType scott(Range)(Range r)
        if (isIterable!Range)
    {
        import core.lifetime: move;
        import mir.math.common: pow, sqrt;
        import mir.math.stat: meanType, VarianceAccumulator;
        import mir.math.sum: ResolveSummationType, Summation;
        import mir.primitives: elementCount;

        alias F = meanType!(Range);
        auto varianceAccumulator = VarianceAccumulator!(F, varianceAlgo, ResolveSummationType!(Summation.appropriate, Range, F))(r);
        auto h = 3.49 * varianceAccumulator.variance!F(false).sqrt / (pow(varianceAccumulator.count, 1.0 / 3));
        assert(h > 0, "scott: bin width must be greater than zero");

        return binsFromWidth!CountType(r, h);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    varianceAlgo = Algorithm used to calculate variance
+/
template scott(VarianceAlgo varianceAlgo = VarianceAlgo.online)
{
    import mir.internal.utility: isFloatingPoint;
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice;
    import std.traits: isIterable, Unqual;

    /++
    Params:
        r = range
    +/
    size_t scott(Range)(Range r)
        if (isIterable!Range)
    {
        return .scott!(size_t, varianceAlgo)(r);
    }
}

/// Scott example
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11].sliced;

    auto k = x.scott!size_t;
    assert(k == 3);
    static assert(is(typeof(k) == size_t));

    auto l = x.scott;
    assert(l == 3);
}

// withAsSlice test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.rc.array: RCArray;

    static immutable a = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

    auto x = RCArray!double(12);
    foreach(i, ref e; x)
        e = a[i];

    auto k = x.scott!size_t;
    assert(k == 3);
}

// dynamic array test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

    auto k = x.scott!size_t;
    assert(k == 3);
}

/++
Computes the number of breaks for a histogram using Freedman-Diaconis' choice.

Calculates the break width using two times the interquartile range of the input
divided by the cubed root of the number of elements of the input. The number of
breaks is then calculated using `binsFromWidth`.

When the interquartile range is zero, this function makes the same adjustments
as R's `nclass.FD` function, which steadily widens the range used in the
interquartile range to a maximum width of [1/512, 511/512] (with an adjustment
to the multiplier as well).

It is less sensitive to the standard deviation of outliers than Scott's rule.

Params:
    CountType = the type that is used to count in histogram bins
    quantileAlgo = the quantile algorithm used (default: QuantileAlgo.type7)
    allowModifySlice = controls whether the input is modified in place (default: false)
Returns:
    The number of breaks

See_also: 
    $(LREF struges),
    $(LREF freedmanDiaconis),
    $(LREF binsFromWidth),
    $(WEB en.wikipedia.org/wiki/Histogram, Histogram)
+/
template freedmanDiaconis(CountType,
                          QuantileAlgo quantileAlgo = QuantileAlgo.type7,
                          bool allowModifySlice = false)
    if(isIntegral!CountType)
{
    import mir.internal.utility: isFloatingPoint;
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice;
    import std.traits: Unqual;

    /++
    Params:
        slice = slice
    +/
    CountType freedmanDiaconis(Iterator, size_t N, SliceKind kind)(
            Slice!(Iterator, N, kind) slice)
    {
        import core.lifetime: move;
        import mir.math.common: pow;
        import mir.ndslice.topology: flattened;
        import mir.primitives: elementCount;
        import mir.stat.descriptive.univariate: interquartileRange, quantileType;

        alias F = quantileType!(Slice!(Iterator), quantileAlgo);
        // This allows the interquartileRange to be calculated over the slice
        // multiple times using allowModifySlice=true. This helps reduce some
        // work
        static if (!allowModifySlice) {
            import mir.ndslice.allocation: rcslice;
            import mir.ndslice.topology: as;

            auto view = slice.lightScope;
            auto val = view.as!(Unqual!(slice.DeepElement)).rcslice;
            auto temp = val.lightScope.flattened;
        } else {
            auto temp = slice.flattened;
        }
        auto iqr = temp.interquartileRange!(F, quantileAlgo, true);
        typeof(iqr) h;
        if (iqr == 0) {
            ushort divisor = 8;
            F q;
            while (iqr == 0 && divisor <= 512)
            {
                q = 1.0 / divisor;
                iqr = temp.interquartileRange!(F, quantileAlgo, true)(q, 1 - q);
                divisor = cast(ushort) (divisor << 1);
            }
            assert(iqr > 0, "freedmanDiaconis: interquartile range must be larger than zero");
            // R makes this adjustment to the multiplier, instead of times 2 for
            // this case
            h = iqr / (1 - 2 * q);
        } else {
            h = 2 * iqr;
        }
        h /= pow(temp.elementCount, 1.0 / 3);
        return binsFromWidth!CountType(temp.move, h);
    }

    /++
    Params:
        array = array
    +/
    CountType freedmanDiaconis(T)(T[] array)
        if (isFloatingPoint!(Unqual!T))
    {
        import mir.ndslice.slice: sliced;

        return freedmanDiaconis(array.sliced);
    }

    /++
    Params:
        withAsSlice = withAsSlice
    +/
    CountType freedmanDiaconis(T)(T withAsSlice)
        if (hasAsSlice!T)
    {
        return freedmanDiaconis(withAsSlice.asSlice);
    }
}

template freedmanDiaconis(QuantileAlgo quantileAlgo = QuantileAlgo.type7,
                          bool allowModifySlice = false)
{
    import mir.ndslice.slice: Slice, SliceKind, hasAsSlice;

    /++
    Params:
        slice = slice
    +/
    size_t freedmanDiaconis(Iterator, size_t N, SliceKind kind)(
            Slice!(Iterator, N, kind) slice)
    {
        return .freedmanDiaconis!(size_t, quantileAlgo, allowModifySlice)(slice);
    }

    /++
    Params:
        array = array
    +/
    size_t freedmanDiaconis(T)(T[] array)
        if (isFloatingPoint!(Unqual!T))
    {
        import mir.ndslice.slice: sliced;

        return .freedmanDiaconis!(size_t, quantileAlgo, allowModifySlice)(array.sliced);
    }

    /++
    Params:
        withAsSlice = withAsSlice
    +/
    size_t freedmanDiaconis(T)(T withAsSlice)
        if (hasAsSlice!T)
    {
        return .freedmanDiaconis!(size_t, quantileAlgo, allowModifySlice)(withAsSlice.asSlice);
    }
}

/// freedmanDiaconis Example
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7].sliced;

    auto k = x.freedmanDiaconis!size_t;
    assert(k == ((7.0 - 0.0) / 3.5));
    static assert(is(typeof(k) == size_t));

    auto l = x.freedmanDiaconis;
    assert(l == ((7.0 - 0.0) / 3.5));
}

// Test example vs. R in extreme case
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.math.common: ceil, pow;
    import mir.stat.descriptive: quantile;
    import mir.ndslice.allocation: slice;

    auto x = slice!double([500], 0.0);
    x[0] = 1;

    auto k = x.freedmanDiaconis!size_t;
    assert(k == ceil(((1.0 - 0.0) / ((x.quantile(1.0 - 1.0 / 512) / (1.0 - 2.0 / 512)) / pow(500.0, 1.0 / 3)))));
}

// withAsSlice test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.rc.array: RCArray;

    static immutable a = [0.0, 1, 2, 3, 4, 5, 6, 7];

    auto x = RCArray!double(8);
    foreach(i, ref e; x)
        e = a[i];

    auto k = x.freedmanDiaconis!size_t;
    assert(k == ((7.0 - 0.0) / 3.5));
}

// dynamic array test
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;

    auto x = [0.0, 1, 2, 3, 4, 5, 6, 7];

    auto k = x.freedmanDiaconis!size_t;
    assert(k == ((7.0 - 0.0) / 3.5));
}
