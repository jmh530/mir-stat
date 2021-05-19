/++
This module contains an API for creating reference-counted histograms.

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

module mir.stat.descriptive.histogram.api.rc;

import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;


import std.meta: allSatisfy;
import std.range.primitives: isRandomAccessRange;
import mir.primitives: DeepElementType;
import mir.rc.array: RCI;
import mir.ndslice.slice: Slice, SliceKind;
import mir.stat.descriptive.histogram.axis: AxisOptions;
import mir.stat.descriptive.histogram.traits: isAxis;

/++
Params:
    slice = slice
    axis = axis
+/
HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
    rchistogramImplBasic(Iterator, size_t N, SliceKind kind, Axis)(
               Slice!(Iterator, N, kind) x, Axis axis)
    if (isAxis!Axis)
{
    import mir.ndslice.allocation: mininitRcslice;

    auto counts = mininitRcslice!(Axis.CountType)(axis.N_bin);
    foreach(ref e; counts) {
        e = 0;
    }
    auto h = HistogramAccumulator!(typeof(counts), Axis)(counts, axis);
    h.put(x);
    return h;
}

// Check rchistogramImplBasic
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.rc.array: RCI;
    import mir.ndslice.slice: sliced;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions())(5, 2.0);
    auto x = [2.0, 2.5, 3.0, 3.5].sliced;

    auto h = rchistogramImplBasic(x, integralAxis);
    assert(h.counts == [2, 2, 0, 0, 0]);
    static assert(is(typeof(h.counts) == Slice!(RCI!(integralAxis.CountType))));
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    axisOptions = options
+/
private
template rchistogramImpl(CountType, BinType, alias Axis, AxisOptions axisOptions)
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: CategoryAxis, IntegralAxis, RegularAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), IntegralAxis!(CountType, BinType, axisOptions))
        rchistogramImpl(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low)
        if (__traits(isSame, Axis, IntegralAxis))
    {
        import core.lifetime: move;

        auto integralAxis = IntegralAxis!(CountType, BinType, axisOptions)(N_bin, low);
        return .rchistogramImplBasic(slice.move, integralAxis);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), RegularAxis!(CountType, BinType, axisOptions))
        rchistogramImpl(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, RegularAxis))
    {
        import core.lifetime: move;

        auto regularAxis = RegularAxis!(CountType, BinType, axisOptions)(N_bin, low, high);
        return .rchistogramImplBasic(slice.move, regularAxis);
    }

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), CategoryAxis!(CountType, BinType, axisOptions))
        rchistogramImpl(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, CategoryAxis))
    {
        import core.lifetime: move;

        CategoryAxis!(CountType, BinType, axisOptions) categoryAxis;
        return .rchistogramImplBasic(slice.move, categoryAxis);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
+/
private
template rchistogramImpl(CountType, BinType, alias Axis)
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: EnumAxis;

    ///
    HistogramAccumulator!(Slice!(RCI!(CountType)), EnumAxis!(CountType, BinType))
        rchistogramImpl(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, EnumAxis))
    {
        import core.lifetime: move;

        EnumAxis!(CountType, BinType) enumAxis;
        return .rchistogramImplBasic(slice.move, enumAxis);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    transform = function to transform axis
    inverseTransform = function to undo transform
    axisOptions = options
+/
private
template rchistogramImpl(CountType, BinType, alias Axis, alias transform, alias inverseTransform, AxisOptions axisOptions)
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: TransformAxis;

    /++
    Params:
        slice = slice
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions))
        rchistogramImpl(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, TransformAxis))
    {
        import core.lifetime: move;

        auto transformAxis = TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
        return .rchistogramImplBasic(slice.move, transformAxis);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Iterator = iterator used in slice
    Axis = type of axis
    axisOptions = options
+/
private
template rchistogramImpl(CountType, Iterator, alias Axis, AxisOptions axisOptions)
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: VariableAxis;
    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/

    HistogramAccumulator!(Slice!(RCI!(CountType)), VariableAxis!(CountType, Iterator, axisOptions))
        rchistogramImpl(size_t N, SliceKind kindA, SliceKind kindB)(
                   Slice!(Iterator, N, kindA) dataSlice,
                   Slice!(Iterator, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis))
    {
        import core.lifetime: move;

        auto variableAxis = VariableAxis!(CountType, Iterator, axisOptions)(axisSlice.move);
        return .rchistogramImplBasic(dataSlice.move, variableAxis.move);
    }
}

/++
Computes a histogram of the inputs.

If the `Axis` has an `options` member, the histogram may optionally allow
for overflow and underflow members.

Params:
    slice = slice
    axis = axis

See_also:
    $(LREF HistogramAccumulator),
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
    rchistogram(Iterator, size_t N, SliceKind kind, Axis)(
               Slice!(Iterator, N, kind) slice, Axis axis)
    if (isAxis!Axis)
{
    import core.lifetime: move;
    return .rchistogramImplBasic(slice.move, axis);
}

/++
Params:
    Axis = type of axis
+/
template rchistogram(Axis)
    if (isAxis!Axis)
{
    import mir.stat.descriptive.histogram.axis: IntegralAxis, RegularAxis,
        TransformAxis, EnumAxis, CategoryAxis, VariableAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        axisOptions = options
    +/
    HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
        rchistogram(Iterator, size_t N, SliceKind kind, CountType, BinType)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low)
        if (is(Axis == IntegralAxis))
    {
        import core.lifetime: move;

        auto integralAxis = Axis(N_bin, low);
        return .rchistogramImplBasic(slice.move, integralAxis);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
        axisOptions = options
    +/
    HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
        rchistogram(Iterator, size_t N, SliceKind kind, CountType, BinType, AxisOptions)(
                    Slice!(Iterator, N, kind) slice,
                    CountType N_bin,
                    BinType low,
                    BinType high)
        if (is(Axis == RegularAxis) || is(Axis == TransformAxis))
    {
        import core.lifetime: move;

        auto axis = Axis(N_bin, low, high);
        return .rchistogramImplBasic(slice.move, axis);
    }

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
        rchistogram(Iterator, size_t N, SliceKind kind)(
                    Slice!(Iterator, N, kind) slice)
        if (is(Axis == EnumAxis) || is(Axis == CategoryAxis))
    {
        import core.lifetime: move;

        auto axis = Axis();
        return .rchistogramImplBasic(slice.move, axis);
    }

    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/
    HistogramAccumulator!(Slice!(RCI!(Axis.CountType)), Axis)
        rchistogram(size_t N, SliceKind kindA, SliceKind kindB)(
                    Slice!(Iterator, N, kindA) dataSlice,
                    Slice!(Iterator, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis))
    {
        import core.lifetime: move;

        auto axis = Axis(axisSlice.move);
        return .rchistogramImplBasic(dataSlice.move, axis.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, BinType, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: IntegralAxis, RegularAxis, CategoryAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), IntegralAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low)
        if (__traits(isSame, Axis, IntegralAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), RegularAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, RegularAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low, high);
    }
    
    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), CategoryAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, CategoryAxis) && is(BinType == enum))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    transform = function to transform axis
    axisOptions = options
+/
template rchistogram(CountType, BinType, alias Axis, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: TransformAxis;

    /++
    Params:
        slice = slice
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, TransformAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, transform, inverseTransform, axisOptions)(slice.move, N_bin, low, high);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Iterator = iterator used in slice
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, Iterator, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: VariableAxis;

    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), VariableAxis!(CountType, Iterator, axisOptions))
        rchistogram(size_t N, SliceKind kindA, SliceKind kindB)(
                   Slice!(Iterator, N, kindA) dataSlice,
                   Slice!(Iterator, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, Iterator, Axis, axisOptions)(dataSlice.move, axisSlice.move);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(BinType, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: IntegralAxis, RegularAxis, CategoryAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), IntegralAxis!(DefaultCountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   DefaultCountType N_bin,
                   BinType low)
        if (__traits(isSame, Axis, IntegralAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, BinType, Axis, axisOptions)(slice.move, N_bin, low);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), RegularAxis!(DefaultCountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   DefaultCountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, RegularAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, BinType, Axis, axisOptions)(slice.move, N_bin, low, high);
    }
    
    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), CategoryAxis!(DefaultCountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, CategoryAxis) && is(BinType == enum))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, BinType, Axis, axisOptions)(slice.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
    transform = function to transform axis
    inverseTransform = function to undo transform
    axisOptions = options
+/
template rchistogram(BinType, alias Axis, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: TransformAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;
    

    /++
    Params:
        slice = slice
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), TransformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice,
                   DefaultCountType N_bin,
                   BinType low,
                   BinType high)
        if (__traits(isSame, Axis, TransformAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, BinType, Axis, transform, inverseTransform, axisOptions)(slice.move, N_bin, low, high);
    }
}

/++
Params:
    Iterator = iterator used in slice
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(Iterator, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: VariableAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), VariableAxis!(DefaultCountType, Iterator, axisOptions))
        rchistogram(size_t N, SliceKind kindA, SliceKind kindB)(
                   Slice!(Iterator, N, kindA) dataSlice,
                   Slice!(Iterator, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, Iterator, Axis, axisOptions)(dataSlice.move, axisSlice.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: IntegralAxis, RegularAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), IntegralAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, BinType)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low)
        if (__traits(isSame, Axis, IntegralAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), RegularAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, BinType)(
            Slice!(Iterator, N, kind) slice,
            CountType N_bin,
            BinType low,
            BinType high)
        if (__traits(isSame, Axis, RegularAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low, high);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    transform = function to transform axis
    inverseTransform = function to undo transform
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, alias Axis, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: TransformAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, BinType)(
            Slice!(Iterator, N, kind) slice,
            CountType N_bin,
            BinType low,
            BinType high)
        if (__traits(isSame, Axis, TransformAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, transform, inverseTransform, axisOptions)(slice.move, N_bin, low, high);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: VariableAxis;

    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), VariableAxis!(CountType, IteratorB, axisOptions))
        rchistogram(IteratorA, size_t N, SliceKind kindA, IteratorB, SliceKind kindB)(
                   Slice!(IteratorA, N, kindA) dataSlice,
                   Slice!(IteratorB, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis) && 
            is(DeepElementType!(Slice!(IteratorA, N, kindA)) : DeepElementType!(Slice!(IteratorB, 1, kindB))))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, IteratorB, Axis, axisOptions)(dataSlice.move, axisSlice.move);
    }
}

/++
Params:
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: CategoryAxis, IntegralAxis, RegularAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), IntegralAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, CountType, BinType)(
                   Slice!(Iterator, N, kind) slice,
                   CountType N_bin,
                   BinType low)
        if (__traits(isSame, Axis, IntegralAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low);
    }

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), RegularAxis!(CountType, BinType, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, CountType, BinType)(
            Slice!(Iterator, N, kind) slice,
            CountType N_bin,
            BinType low,
            BinType high)
        if (__traits(isSame, Axis, RegularAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, axisOptions)(slice.move, N_bin, low, high);
    }
    
    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), CategoryAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, CategoryAxis) && is(DeepElementType!(typeof(slice)) == enum))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, DeepElementType!(typeof(slice)), Axis, axisOptions)(slice.move);
    }
}

/++
Params:
    Axis = type of axis
+/
template rchistogram(alias Axis)
    if (__traits(isTemplate, Axis))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: EnumAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), EnumAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind))))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, EnumAxis) && is(DeepElementType!(typeof(slice)) == enum))
    {
        import core.lifetime: move;
        import mir.stat.descriptive.histogram.traits: DefaultCountType;

        return .rchistogramImpl!(DefaultCountType, DeepElementType!(typeof(slice)), Axis)(slice.move);
    }
}

/++
Params:
    Axis = type of axis
    transform = function to transform axis
    inverseTransform = function to undo transform
    axisOptions = options
+/
template rchistogram(alias Axis, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.stat.descriptive.histogram.axis: TransformAxis;

    /++
    Params:
        slice = slice
        N_bin = number of bins
        low = the value of the smallest bin
        high = the value of the largest bin
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind, CountType, BinType)(
            Slice!(Iterator, N, kind) slice,
            CountType N_bin,
            BinType low,
            BinType high)
        if (__traits(isSame, Axis, TransformAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, BinType, Axis, transform, inverseTransform, axisOptions)(slice.move, N_bin, low, high);
    }
}

/++
Params:
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: VariableAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        dataSlice = slice of data
        axisSlice = slice of axis breaks
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), VariableAxis!(DefaultCountType, IteratorB, axisOptions))
        rchistogram(IteratorA, size_t N, SliceKind kindA, IteratorB, SliceKind kindB)(
                   Slice!(IteratorA, N, kindA) dataSlice,
                   Slice!(IteratorB, 1, kindB) axisSlice)
        if (__traits(isSame, Axis, VariableAxis) && 
            is(DeepElementType!(Slice!(IteratorA, N, kindA)) : DeepElementType!(Slice!(IteratorB, 1, kindB))))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, IteratorB, Axis, axisOptions)(dataSlice.move, axisSlice.move);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    Axis = type of axis
+/
template rchistogram(BinType, alias Axis)
    if (__traits(isTemplate, Axis) && is(BinType == enum))
{
    import mir.stat.descriptive.histogram.axis: EnumAxis;
    import mir.stat.descriptive.histogram.traits: DefaultCountType;

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(DefaultCountType)), EnumAxis!(DefaultCountType, BinType))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, EnumAxis))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(DefaultCountType, BinType, Axis)(slice.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Axis = type of axis
+/
template rchistogram(CountType, alias Axis)
    if (__traits(isTemplate, Axis) && !is(CountType == enum))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: EnumAxis;

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), EnumAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind))))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, EnumAxis) && is(DeepElementType!(typeof(slice)) == enum))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, DeepElementType!(typeof(slice)), Axis)(slice.move);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    Axis = type of axis
    axisOptions = options
+/
template rchistogram(CountType, alias Axis, AxisOptions axisOptions = AxisOptions())
    if (__traits(isTemplate, Axis) && !is(CountType == enum))
{
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: CategoryAxis;

    /++
    Params:
        slice = slice
    +/
    HistogramAccumulator!(Slice!(RCI!(CountType)), CategoryAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), axisOptions))
        rchistogram(Iterator, size_t N, SliceKind kind)(
                   Slice!(Iterator, N, kind) slice)
        if (__traits(isSame, Axis, CategoryAxis) && is(DeepElementType!(typeof(slice)) == enum))
    {
        import core.lifetime: move;

        return .rchistogramImpl!(CountType, DeepElementType!(typeof(slice)), Axis, axisOptions)(slice.move);
    }
}

/// Integral Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: IntegralAxis, integralAxis;
    import mir.stat.descriptive.histogram.breaks: sturges;

    static immutable a = [0.0, 0.5, 1, 1.5, 2];
    static immutable b = [2, 2, 1];
    static immutable c = [2, 2, 1, 0];

    auto x = rcslice!double(a);
    auto result1 = rcslice!size_t(b);
    auto result2 = rcslice!size_t(c);

    auto h1 = x.rchistogram!IntegralAxis(3u, 0.0);
    assert(h1.counts == result1);
    static assert(is(h1.CountType == uint));

    // Pass axis directly
    auto iAxis = integralAxis(3u, 0.0);
    auto h2 = x.rchistogram(iAxis);
    assert(h2.counts == result1);

    // Use function to calculate N_bin
    auto iAxis2 = x.integralAxis!sturges(0.0);
    auto h3 = x.rchistogram(iAxis2);
    assert(h3.counts == result2);
}

/// Regular Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: RegularAxis, regularAxis;
    import mir.stat.descriptive.histogram.breaks: sturges;

    static immutable a = [0.0, 1, 4, 5, 6, 9, 10, 13, 14];
    static immutable b = [3, 3, 3];
    static immutable c = [2, 2, 1, 2, 2];

    auto x = rcslice!double(a);
    auto result1 = rcslice!size_t(b);
    auto result2 = rcslice!size_t(c);

    auto h1 = x.rchistogram!RegularAxis(3u, 0.0, 15.0);
    assert(h1.counts == result1);
    static assert(is(h1.CountType == uint));

    // Pass axis directly
    auto regularAxis2 = regularAxis(3u, 0.0, 15.0);
    auto h2 = x.rchistogram(regularAxis2);
    assert(h2.counts == result1);

    // Use function to calculate N_bin
    auto regularAxis3 = x.regularAxis!sturges(0.0, 15.0);
    auto h3 = rchistogram(x, regularAxis3);
    assert(h3.counts == result2);
}

/// Transform Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: TransformAxis, transformAxis, inverseTransformMapping;
    import mir.stat.descriptive.histogram.breaks: sturges;

    
    static immutable a = [10.0 ^^ 2.0, 10.0 ^^ 2.5, 10.0 ^^ 5.0, 10.0 ^^ 11.5];
    static immutable b = [2, 1, 0, 1];
    static immutable c = [3, 0, 1];
    
    auto x = rcslice!double(a);
    auto result1 = rcslice!size_t(b);
    auto result2 = rcslice!size_t(c);

    auto h1 = x.rchistogram!(TransformAxis, log10, inverseTransformMapping!log10)(4u, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(h1.counts == result1);
    static assert(is(h1.CountType == uint));

    // Pass axis directly
    auto regularAxis2 = transformAxis!(log10, inverseTransformMapping!log10)(4u, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    auto h2 = x.rchistogram(regularAxis2);
    assert(h2.counts == result1);
    
    // Use function to calculate N_bin
    auto regularAxis3 = x.transformAxis!(log10, inverseTransformMapping!log10, sturges)(10.0 ^^ 2.0, 10.0 ^^ 12.0);
    auto h3 = x.rchistogram(regularAxis3);
    assert(h3.counts == result2);

    // Can also supply lambda
    auto h4 = x.rchistogram!(TransformAxis, a => log10(a), a => (10.0 ^^ a))(4u, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(h4.counts == result1);

    // Or string lambda
    auto h5 = x.rchistogram!(TransformAxis, "log10(a)", "10.0 ^^ a")(4u, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(h5.counts == result1);
}

/// Enum Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: EnumAxis, enumAxis;

    enum Foo {
        A,
        B
    }

    static immutable a = [Foo.A, Foo.B, Foo.A, Foo.A, Foo.B];
    static immutable b = [3, 2];

    auto x = rcslice!Foo(a);
    auto result = rcslice!size_t(b);

    auto h1 = x.rchistogram!EnumAxis;
    assert(h1.counts == result);
    static assert(is(h1.CountType == size_t));

    // Pass axis directly
    auto eAxis = enumAxis!Foo;
    auto h2 = x.rchistogram(eAxis);
    assert(h2.counts == result);
}

/// Category Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: CategoryAxis, categoryAxis;

    enum Foo {
        A,
        B
    }

    static immutable a = [Foo.A, Foo.B, Foo.A, Foo.A, Foo.B];
    static immutable b = ["A", "B", "A", "A", "B"];
    static immutable c = [3, 2];

    auto x = rcslice!Foo(a);
    auto y = rcslice!string(b);
    auto result = rcslice!size_t(c);

    auto h1 = x.rchistogram!CategoryAxis;
    assert(h1.counts == result);
    static assert(is(h1.CountType == size_t));

    // Can handle string inputs
    auto h2 = y.rchistogram!(Foo, CategoryAxis);
    assert(h2.counts == result);

    // Pass axis directly
    auto cAxis = categoryAxis!Foo();
    auto h3 = x.rchistogram(cAxis);
    assert(h3.counts == result);
}

/// Variable Axis example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.primitives: DeepElementType;
    import mir.stat.descriptive.histogram.axis: VariableAxis, variableAxis;

    static immutable a = [0.0, 0.5, 1, 1.5, 2];
    static immutable b = [0.0, 1, 2, 3];
    static immutable c = [2, 2, 1];

    auto x = rcslice!double(a);
    auto breaks = rcslice!double(b);
    auto result = rcslice!size_t(c);

    auto h1 = x.rchistogram!VariableAxis(breaks);
    assert(h1.counts == result);
    static assert(is(h1.CountType == size_t));

    // Pass axis directly
    auto vAxis = variableAxis(breaks);
    auto h2 = x.rchistogram(vAxis);
    assert(h2.counts == result);
}
