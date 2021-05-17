/++
This module contains a histogram accumulator.

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

module mir.stat.descriptive.histogram.accumulator;

import std.meta: allSatisfy;
import mir.stat.descriptive.histogram.traits: isAxis;

struct DenseStorage(Storage)
{
    import std.traits: isNumeric;

    Storage storage; 
    
    void put(size_t i = 0)(size_t x)
        if (is(Storage : T[], T) ||
            is(Storage : T[N], T, size_t N))
    {
        storage[i]++;
    }
    
    void put(size_t i = 0)(size_t x)
        if (isNumeric!Storage)
    {
        storage++;
    }

    void put(size_t i = 0)()
        if (is(Storage : T[], T) ||
            is(Storage : T[N], T, size_t N))
    {
        storage[i]++;
    }

    void put(size_t i = 0)()
        if (isNumeric!Storage)
    {
        storage++;
    }
}

/++
Accumulator used to generate histogram.

If the `Axis` has an `options` member, the histogram may optionally allow
for overflow and underflow members.

Params:
    Axis = the type of the axis used to create the histogram bins

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis),
    $(LREF FrequencyAccumulator)
+/
struct HistogramAccumulator(Storage, Axis...)
    if (Axis.length > 0 &&
        allSatisfy!(isAxis, Axis))
{
    import std.meta: allSatisfy, staticMap;
    import std.traits: hasMember, isIterable, isSomeString;
    import mir.primitives: hasShape, DeepElementType;
    import mir.stat.descriptive.histogram.traits: includeOverflow, includeUnderflow,
        BinTypeOf, isCategoryAxis;
// need to check that storage dimension matches axis
// isRandomAccessRange for storage does not work for single value
//isRandomAccessRange!Storage &&
//        __traits(compiles, {alias deepElementType = DeepElementType!(Storage);})&&

//1) Need to handle template constraints properly
//2) over/underflow is not currently handling multiple axes properly, just checking the first axis
//need to do any allow over/underflow, then just put them all
//3) Need to be able to put another dense storage
private:
    static if (includeOverflow!Axis)
    {
        static if (N == 1) {
            ///
            DenseStorage!CountType overflowStorage;
        } else {
            ///
            DenseStorage!(CountType[N]) overflowStorage;
        }
    }

    static if (includeUnderflow!Axis)
    {
        static if (N == 1) {
            ///
            DenseStorage!CountType underflowStorage;
        } else {
            ///
            DenseStorage!(CountType[N]) underflowStorage;
        }
    }

public:

    ///
    Axis axis;

    ///
    Storage counts;

    //
    enum N = Axis.length;

    ///
    alias CountType = DeepElementType!Storage;
    
    static if (includeOverflow!Axis)
    {
        ///
        alias OverflowType = typeof(overflowStorage.storage);
    }

    static if (includeUnderflow!Axis)
    {
        ///
        alias UnderflowType = typeof(underflowStorage.storage);
    }

    ///
    this(Storage x, Axis y)
    {
        counts = x;
        axis = y;
    }

    ///
    void put(Range)(Range r)
        if (N == 1 && 
            isIterable!Range && 
            !(isCategoryAxis!(Axis[0]) && isSomeString!Range))
    {
        foreach(x; r)
        {
            put(x);
        }
    }

    ///
    void put(T...)(T x)
    {
        import mir.ndslice.topology: iota;

        static foreach(i; iota(N))
        {
            putSingleImpl!(T[i], i)(x[i]);
        }
    }

    private
    void putSingleImpl(T, size_t i)(T x)
        if (is(T == BinTypeOf!(Axis[i])) || 
            (isCategoryAxis!(Axis[i]) && isSomeString!(T)))
    {
        static if (includeOverflow!(Axis[i]) && includeUnderflow!(Axis[i])) {
            if (axis[i].isOverflow(x)) {
                overflowStorage.put!i();
            } else if (axis[i].isUnderflow(x)) {
                underflowStorage.put!i();
            } else {
                putStorage!(T, i)(x);
            }
        } else static if (!includeOverflow!(Axis[i]) && includeUnderflow!(Axis[i])) {
            if (axis[i].isUnderflow(x)) {
                underflowStorage.put!i();
            } else {
                putStorage!(T, i)(x);
            }
        } else static if (includeOverflow!(Axis[i]) && !includeUnderflow!(Axis[i])) {
            if (axis[i].isOverflow(x)) {
                overflowStorage.put!i();
            } else {
                putStorage!(T, i)(x);
            }
        } else {
            putStorage!(T, i)(x);
        }
    }

    private
    void putStorage(T, size_t i)(T x)
    {
        counts[axis[i].index(x)]++;
    }

    ///
    void put(HistogramAccumulator!(Storage, Axis) h)
    {
        import mir.stat.descriptive.histogram.traits: hasAxisOptions;

        assert(axis == h.axis);
        counts[] += h.counts[];
        static if (hasAxisOptions!(Axis[0])) {
            static if (Axis[0].options.enableOverflow && hasMember!(typeof(h), "overflowStorage")) {
                static if (N == 1) {
                    overflowStorage.storage += h.overflowStorage.storage;
                } else {
                    overflowStorage.storage[] += h.overflowStorage.storage[];
                }
            }
            static if (Axis[0].options.enableUnderflow && hasMember!(typeof(h), "underflowStorage")) {
                static if (N == 1) {
                    underflowStorage.storage += h.underflowStorage.storage;
                } else {
                    underflowStorage.storage[] += h.underflowStorage.storage[];
                }
            }
        }
    }

    static if (includeOverflow!Axis)
    {
        ///
        OverflowType overflow()()
        {
            return overflowStorage.storage;
        }
    }

    static if (includeUnderflow!Axis)
    {
        ///
        UnderflowType underflow()()
        {
            return underflowStorage.storage;
        }
    }
}

// Check IntegralAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions())(5, 2.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(integralAxis))(counts, integralAxis);
    h.put([2.0, 2.5, 3.0, 3.5]);
    assert(counts == [2, 2, 0, 0, 0]);
    h.put(4.0);
    assert(counts == [2, 2, 1, 0, 0]);
}

// Check over/underflow IntegralAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions(false, true, true))(5, 2.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(integralAxis))(counts, integralAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
    h.put(3.0);
    assert(counts == [0, 1, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check over/underflow IntegralAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis, EnableUnderflow;

    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions(EnableUnderflow(true)))(5, 2.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(integralAxis))(counts, integralAxis);
    assert(h.underflow == 0);
    h.put(1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.underflow == 1);
    h.put(3.0);
    assert(counts == [0, 1, 0, 0, 0]);
    assert(h.underflow == 1);
}

// Check EnumAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, EnumAxis;

    enum Foo {
        A,
        B
    }
    EnumAxis!(size_t, Foo) enumAxis;
    size_t[] counts = [0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(enumAxis))(counts, enumAxis);
    h.put([Foo.A, Foo.B, Foo.B, Foo.B]);
    assert(counts == [1, 3]);
    h.put(Foo.A);
    assert(counts == [2, 3]);
}

// Check CategoryAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, CategoryAxis;

    enum Foo {
        A,
        B
    }
    CategoryAxis!(size_t, Foo, AxisOptions(false, true)) categoryAxis;
    size_t[] counts = [0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(categoryAxis))(counts, categoryAxis);
    h.put([Foo.A, Foo.B, Foo.B, Foo.B]);
    assert(counts == [1, 3]);
    h.put(Foo.A);
    assert(counts == [2, 3]);

    // Check strings
    h.put("B");
    assert(counts == [2, 4]);
    assert(h.overflow == 0);
    h.put("C");
    assert(h.overflow == 1);
    h.put(["C", "D"]);
    assert(h.overflow == 3);
    h.put(["CD"]);
    assert(h.overflow == 4);
}

// Check RegularAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, RegularAxis;

    auto regularAxis = RegularAxis!(size_t, double, AxisOptions())(5, 2.0, 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(regularAxis))(counts, regularAxis);
    h.put([2.0, 2.5, 3.0, 11.5]);
    assert(counts == [3, 0, 0, 0, 1]);
    h.put(7.0);
    assert(counts == [3, 0, 1, 0, 1]);
}

// Check over/underflow RegularAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, RegularAxis;

    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(false, true, true))(5, 2.0, 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(regularAxis))(counts, regularAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check RegularAxis, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, RegularAxis;

    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(true))(5, 2.0, 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(regularAxis))(counts, regularAxis);
    h.put([2.5, 3.0, 3.5, 12.0]);
    assert(counts == [3, 0, 0, 0, 1]);
    h.put(7.0);
    assert(counts == [3, 0, 1, 0, 1]);
}

// Check over/underflow RegularAxis, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, RegularAxis;

    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(true, true, true))(5, 2.0, 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(regularAxis))(counts, regularAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check TransformAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.math.common: log10;
    import mir.stat.descriptive.histogram.axis: AxisOptions, TransformAxis, inverseTransformMapping;

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseTransformMapping!log10, AxisOptions())(5, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(transformAxis))(counts, transformAxis);
    h.put([10.0 ^^ 2.0, 10.0 ^^ 2.5, 10.0 ^^ 3.0, 10.0 ^^ 11.5]);
    assert(counts == [3, 0, 0, 0, 1]);
    h.put(10.0 ^^ 7.0);
    assert(counts == [3, 0, 1, 0, 1]);
}

// Check over/underflow TransformAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.math.common: log10;
    import mir.stat.descriptive.histogram.axis: AxisOptions, TransformAxis, inverseTransformMapping;

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseTransformMapping!log10, AxisOptions(false, true, true))(5, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(transformAxis))(counts, transformAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(10.0 ^^ 13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(10.0 ^^ 1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check Transform, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.math.common: log10;
    import mir.stat.descriptive.histogram.axis: AxisOptions, TransformAxis, inverseTransformMapping;

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseTransformMapping!log10, AxisOptions(true))(5, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(transformAxis))(counts, transformAxis);
    h.put([10.0 ^^ 2.5, 10.0 ^^ 3.0, 10.0 ^^ 3.5, 10.0 ^^ 12.0]);
    assert(counts == [3, 0, 0, 0, 1]);
    h.put(10.0 ^^ 7.0);
    assert(counts == [3, 0, 1, 0, 1]);
}

// Check over/underflow TransformAxis, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.math.common: log10;
    import mir.stat.descriptive.histogram.axis: AxisOptions, TransformAxis, inverseTransformMapping;

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseTransformMapping!log10, AxisOptions(true, true, true))(5, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(transformAxis))(counts, transformAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(10.0 ^^ 13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(10.0 ^^ 1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check VariableAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.stat.descriptive.histogram.axis: AxisOptions, VariableAxis;

    auto axisSlice = [2.0, 3, 4, 5, 6, 7].sliced;
    auto variableAxis = VariableAxis!(size_t, double*, AxisOptions())(axisSlice);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(variableAxis))(counts, variableAxis);
    h.put([2.0, 2.5, 3.0, 3.5]);
    assert(counts == [2, 2, 0, 0, 0]);
    h.put(4.0);
    assert(counts == [2, 2, 1, 0, 0]);
}

// Check over/underflow VariableAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.stat.descriptive.histogram.axis: AxisOptions, VariableAxis;

    auto axisSlice = [2.0, 3, 4, 5, 6, 7].sliced;
    auto variableAxis = VariableAxis!(size_t, double*, AxisOptions(false, true, true))(axisSlice);
    size_t[] counts = [0, 0, 0, 0, 0];

    auto h = HistogramAccumulator!(size_t[], typeof(variableAxis))(counts, variableAxis);
    assert(h.overflow == 0);
    assert(h.underflow == 0);
    h.put(13.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 0);
    h.put(1.0);
    assert(counts == [0, 0, 0, 0, 0]);
    assert(h.overflow == 1);
    assert(h.underflow == 1);
}

// Check put HistogramAccumulator
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    auto integralAxis1 = IntegralAxis!(size_t, double, AxisOptions())(5, 2.0);
    size_t[] counts1 = [0, 0, 0, 0, 0];
    auto integralAxis2 = IntegralAxis!(size_t, double, AxisOptions())(5, 2.0);
    size_t[] counts2 = [0, 0, 0, 0, 0];

    auto h1 = HistogramAccumulator!(size_t[], typeof(integralAxis1))(counts1, integralAxis1);
    h1.put([2.0, 2.5, 3.0, 3.5, 4.0]);
    assert(counts1 == [2, 2, 1, 0, 0]);
    auto h2 = HistogramAccumulator!(size_t[], typeof(integralAxis2))(counts2, integralAxis2);
    h2.put([4.0, 5.0, 5.5, 6.0, 6.5]);
    assert(counts2 == [0, 0, 1, 2, 2]);
    h2.put(h1);
    assert(counts2 == [2, 2, 2, 2, 2]);
}

// Check put HistogramAccumulator with over/underflow
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis, 
        EnableOverflow, EnableUnderflow;

    auto integralAxis1 = IntegralAxis!(size_t, double, AxisOptions(EnableOverflow(true), EnableUnderflow(true)))(5, 2.0);
    size_t[] counts1 = [0, 0, 0, 0, 0];
    auto integralAxis2 = IntegralAxis!(size_t, double, AxisOptions(EnableOverflow(true), EnableUnderflow(true)))(5, 2.0);
    size_t[] counts2 = [0, 0, 0, 0, 0];

    auto h1 = HistogramAccumulator!(size_t[], typeof(integralAxis1))(counts1, integralAxis1);
    h1.put(-1.0);
    auto h2 = HistogramAccumulator!(size_t[], typeof(integralAxis2))(counts2, integralAxis2);
    h2.put(9.0);
    h2.put(h1);
    assert(h2.overflow == 1);
    assert(h2.underflow == 1);
}

// Check custom CircleAxis
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.functional: RefTuple, refTuple;
    import mir.stat.descriptive.histogram.axis: AxisOptions;
    
    struct Point
    {
        double x;
        double y;
    }

    static struct CircleAxis(AxisOptions axisOptions)
    {
        alias CountType = size_t;
        alias BinType = Point;

        enum CountType N_bin = 1;

        CountType index(BinType x)
        {
            if (!isOverflow(x)) {
                return 0;
            } else {
                assert(0, "index: input may not overflow");
            }
        }

        bool isOverflow()(BinType x) const
        {
            return x.x * x.x + x.y + x.y < 1.0;
        }
    }

    auto circleAxis = CircleAxis!(AxisOptions())();
    size_t[1] count = 0;

    auto h = HistogramAccumulator!(size_t[1], typeof(circleAxis))(count, circleAxis);
    auto p = Point(0.25, 0.5);
    h.put(p);
}
