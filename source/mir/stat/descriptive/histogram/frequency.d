/++
This module contains algorithms for frequency statistics.

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall, Ilya Yaroshenko

Copyright: 2026 Mir Stat Authors.

Macros:
SUBREF = $(REF_ALTTEXT $(TT $2), $2, mir, stat, $1)$(NBSP)
MATHREF = $(REF_ALTTEXT $(TT $2), $2, mir, math, $1)$(NBSP)
NDSLICEREF = $(REF_ALTTEXT $(TT $2), $2, mir, ndslice, $1)$(NBSP)
T2=$(TR $(TDNW $(LREF $1)) $(TD $+))
T4=$(TR $(TDNW $(LREF $1)) $(TD $2) $(TD $3) $(TD $4))
+/


module mir.stat.descriptive.histogram.frequency;

import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
import mir.stat.descriptive.histogram.traits: isAxis;

/++
One-dimensional histogram wrapper that also maintains the total recorded count.
The total includes ordinary bins and enabled underflow and overflow counters.
Storage supplied to the constructor may already contain ordinary bin counts;
flow counters start at zero. The total is calculated once at construction.

Counts are exposed read-only. Arrays and slices retain their usual aliasing
semantics: callers must not modify backing storage through external aliases or
independently mutate copies of this wrapper that share storage. Counter types
must be large enough for both bin counts and the total.

Relative-frequency output is not implemented yet.

Params:
    Storage = storage for ordinary bin counts, as in $(LREF HistogramAccumulator)
    FrequencyType = reserved output type for future relative-frequency access
    AxisType = type of the single histogram axis

See_also:
    $(LREF HistogramAccumulator),
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF TransformAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
struct FrequencyAccumulator(Storage, FrequencyType, AxisType)
    if (isAxis!AxisType)
{
    import std.traits: isIterable, isSomeString;
    import mir.stat.descriptive.histogram.traits: includeOverflow, includeUnderflow,
        BinTypeOf, isCategoryAxis;

    private alias Histogram = HistogramAccumulator!(Storage, AxisType);
    private Histogram histogramAccumulator;
    private CountType total;

    /// Type used for bin counts and the running total.
    alias CountType = Histogram.CountType;

    /++
    Params:
        counts = ordinary bin counts; length must match the axis
        axis = axis used to classify observations
    +/
    this(Storage counts, AxisType axis)
    {
        assert(counts.length == axis.N_bin,
            "FrequencyAccumulator: storage length must match the axis");
        histogramAccumulator = Histogram(counts, axis);
        total = 0;
        foreach (value; counts)
            total += value;
    }

    /// Total recorded observations, including flow bins.
    CountType count() const @property
    {
        return total;
    }

    /// Read-only access to ordinary bin counts.
    ref const(Storage) counts() const @property
    {
        return histogramAccumulator.counts;
    }

    /// Read-only access to the axis.
    ref const(AxisType) axis() const @property
    {
        return histogramAccumulator.axis[0];
    }

    static if (includeOverflow!AxisType)
    {
        /// Recorded overflow observations.
        CountType overflow()() { return histogramAccumulator.overflow; }
    }

    static if (includeUnderflow!AxisType)
    {
        /// Recorded underflow observations.
        CountType underflow()() { return histogramAccumulator.underflow; }
    }

    /// Record an iterable of observations, preserving the total after each put.
    void put(Range)(Range r)
        if (isIterable!Range &&
            !(isCategoryAxis!AxisType && isSomeString!Range))
    {
        foreach (x; r)
            put(x);
    }

    /// Record one observation and increment the total after successful insertion.
    void put(T)(T x)
        if (is(T == BinTypeOf!AxisType) ||
            (isCategoryAxis!AxisType && isSomeString!T))
    {
        histogramAccumulator.put(x);
        total++;
    }

    /// Merge another accumulator with a compatible axis.
    void put(ref FrequencyAccumulator f)
    {
        auto addedCount = f.count;
        histogramAccumulator.put(f.histogramAccumulator);
        total += addedCount;
    }
}

/// Collect observations while maintaining a total alongside the bin counts.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    auto f = FrequencyAccumulator!(uint[], double, Axis)([0u, 0u, 0u], Axis(3, 0.0));
    assert(f.count == 0);

    // Bins are [0, 1), [1, 2), and [2, 3).
    f.put(0.5);
    f.put([1.0, 1.5, 2.5]);
    assert(f.counts == [1u, 2u, 1u]);
    assert(f.count == 4);
}

/// The total includes observations recorded in enabled flow bins.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis,
        EnableOverflow, EnableUnderflow;

    alias Axis = IntegralAxis!(uint, double,
        AxisOptions(EnableOverflow(true), EnableUnderflow(true)));
    auto f = FrequencyAccumulator!(uint[], double, Axis)([0u, 0u], Axis(2, 0.0));
    f.put([-1.0, 0.5, 1.5, 3.0]);
    assert(f.counts == [1u, 1u]);
    assert(f.underflow == 1);
    assert(f.overflow == 1);
    assert(f.count == 4);
}

/// Initialize from existing reference-counted storage and merge another total.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    auto counts = rcslice!uint([2u, 1u]);
    alias F = FrequencyAccumulator!(typeof(counts), double, Axis);
    auto f = F(counts, Axis(2, 0.0));
    assert(f.count == 3);

    // Each accumulator has separate backing storage.
    auto other = F(rcslice!uint([0u, 0u]), Axis(2, 0.0));
    other.put([0.5, 1.5]);
    f.put(other);
    assert(f.counts == [3u, 2u]);
    assert(f.count == 5);
    assert(other.count == 2);
}

// Construction, insertion, and merging with array and reference-counted storage.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.ndslice.slice: sliced;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions(false, true, true));

    void check(Storage)(Storage counts, Storage otherCounts)
    {
        alias F = FrequencyAccumulator!(Storage, double, Axis);
        auto axis = Axis(3, 0.0);
        auto f = F(counts, axis);
        static assert(is(f.CountType == uint));
        static assert(!__traits(compiles, f.counts[0] = 9));
        static assert(!__traits(compiles, f.count = 9));
        assert(f.count == 3);
        assert(f.counts == [1u, 2u, 0u]);
        assert(f.axis.N_bin == 3);
        assert(f.underflow == 0 && f.overflow == 0);

        void checkTotal()
        {
            uint sum = f.underflow + f.overflow;
            foreach (value; f.counts)
                sum += value;
            assert(f.count == sum);
        }
        checkTotal();
        f.put(0.5);
        assert(f.count == 4);
        checkTotal();
        f.put([-1.0, 1.5, 3.0]);
        assert(f.count == 7);
        checkTotal();
        f.put([2.0, 2.5].sliced);
        assert(f.counts == [2u, 3u, 2u]);
        assert(f.underflow == 1 && f.overflow == 1);
        checkTotal();
        f.put((double[]).init);
        assert(f.count == 9);
        checkTotal();

        auto other = F(otherCounts, axis);
        assert(other.count == 0);
        other.put([-2.0, 0.5, 4.0, 5.0]);
        f.put(other);
        assert(f.counts == [3u, 3u, 2u]);
        assert(f.underflow == 2 && f.overflow == 3);
        assert(f.count == 13);
        assert(other.count == 4);
        checkTotal();
    }

    check([1u, 2u, 0u], [0u, 0u, 0u]);
    check(rcslice!uint([1u, 2u, 0u]), rcslice!uint([0u, 0u, 0u]));
}

// No-flow axes and category strings use the same counting paths.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis,
        CategoryAxis, EnableOverflow;
    import std.range: only;

    alias Axis = IntegralAxis!(size_t, double, AxisOptions());
    auto f = FrequencyAccumulator!(size_t[], double, Axis)([0UL, 0], Axis(2, 0.0));
    static assert(!__traits(hasMember, typeof(f), "overflow"));
    static assert(!__traits(hasMember, typeof(f), "underflow"));
    f.put(only(0.5, 1.5, 1.75));
    assert(f.count == 3 && f.counts == [1, 2]);

    enum Label { first, second }
    alias Categories = CategoryAxis!(uint, Label, AxisOptions(EnableOverflow(true)));
    auto c = FrequencyAccumulator!(uint[], double, Categories)([0u, 0u], Categories());
    c.put("first");
    c.put(["second", "unknown"]);
    c.put(Label.second);
    assert(c.counts == [1u, 2u]);
    assert(c.overflow == 1 && c.count == 4);
}

// Rejected input must not inflate the total; a partially accepted range stays consistent.
version(mir_stat_test_hist)
unittest
{
    import core.exception: AssertError;
    import std.exception: assertThrown;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    alias F = FrequencyAccumulator!(uint[], double, Axis);
    auto f = F([0u, 0u], Axis(2, 0.0));
    assertThrown!AssertError(f.put(3.0));
    assert(f.count == 0 && f.counts == [0u, 0u]);
    assertThrown!AssertError(f.put([0.5, 3.0]));
    assert(f.count == 1 && f.counts == [1u, 0u]);

    auto incompatible = F([1u, 0u], Axis(2, 1.0));
    assertThrown!AssertError(f.put(incompatible));
    assert(f.count == 1 && f.counts == [1u, 0u]);
    assertThrown!AssertError(F([0u], Axis(2, 0.0)));
}

// Custom axes can provide flow predicates without an AxisOptions member.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    static struct Axis
    {
        alias CountType = uint;
        alias BinType = int;
        enum N_bin = 1;
        uint index(int x) const { assert(x == 0); return 0; }
        bool isUnderflow(int x) const { return x < 0; }
        bool isOverflow(int x) const { return x > 0; }
    }
    alias F = FrequencyAccumulator!(uint[], double, Axis);
    auto f = F([0u], Axis());
    auto other = F([0u], Axis());
    f.put([0, -1]);
    other.put([0, 1, 2]);
    f.put(other);
    assert(f.counts == [2u]);
    assert(f.underflow == 1 && f.overflow == 2);
    assert(f.count == 5);

    f.put(f);
    assert(f.counts == [4u]);
    assert(f.underflow == 2 && f.overflow == 4);
    assert(f.count == 10);
}
