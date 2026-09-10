/++
Read-only views of histogram bins and counts.

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall

Copyright: 2026 Mir Stat Authors.
+/
module mir.stat.descriptive.histogram.view;

import mir.ndslice.slice: isSlice;
import mir.primitives: DeepElementType;
import mir.stat.descriptive.histogram.traits: isAxis;
import std.traits: isDynamicArray, isNumeric;

/++
A bin description and the count read when the element was accessed.

The index is the original ordinary-bin index, including when the view is sliced.
Elements are returned by value; assigning to count does not update the histogram.

Params:
    BinDescription = type returned by the axis's const bin accessor
    Count = histogram count type
+/
struct HistogramBin(BinDescription, Count)
{
    /// Original ordinary-bin index.
    size_t index;
    /// Axis-specific interval or category description.
    BinDescription bin;
    /// Count at the time this element was read.
    Count count;
}

package template supportsBinView(Storage, Axis)
{
    static if (isDynamicArray!Storage)
        private enum supportedStorage = true;
    else static if (isSlice!Storage)
        private enum supportedStorage = Storage.N == 1;
    else
        private enum supportedStorage = false;

    static if (supportedStorage && isAxis!Axis)
        enum supportsBinView = __traits(compiles, {
            Storage storage;
            const Axis axis;
            auto description = axis.bin(size_t.init);
            auto count = storage[size_t.init];
            static assert(isNumeric!(DeepElementType!Storage));
        });
    else
        enum supportsBinView = false;
}

/++
Read-only random-access range of ordinary bins and counts.

Usually obtained from a histogram's bins accessor. Includes zero-count bins;
underflow and overflow are excluded. Numeric bin descriptions expose low and
high, while enum and category descriptions expose slot. Interval closure and
circular behavior remain defined by the axis options.

The view copies the axis and storage handles without allocating a result array.
Dynamic arrays and one-dimensional Mir slices are supported. Shared count
updates are visible on later reads; previously returned counts are values.
Replacing the source histogram's axis or storage does not redirect the view.
Traversal and slicing change only the view's position.

Reference-counted handles retain their allocations. Borrowed storage, including
borrowed variable-axis breaks, must outlive the view and any bin descriptions
that refer to it. Keep storage shape and shared axis boundaries unchanged.
Custom axes must provide const bin access whose result does not allow mutation
of shared boundaries.

Params:
    Storage = dynamic array or one-dimensional Mir slice of numeric counts
    Axis = axis with const runtime bin-description access
+/
struct HistogramBinView(Storage, Axis)
    if (supportsBinView!(Storage, Axis))
{
    private Storage _counts;
    private Axis _axis;
    private size_t _begin;
    private size_t _end;

    /// Type of each returned value.
    alias Element = HistogramBin!(
        typeof((const Axis).init.bin(size_t.init)), DeepElementType!Storage);

    /// Construct a view with one count per ordinary bin.
    this(Storage counts, Axis axis)
    {
        assert(counts.length == axis.N_bin,
            "HistogramBinView: count length must match axis bin count");
        _counts = counts;
        _axis = axis;
        _begin = 0;
        _end = counts.length;
    }

    /// Number of remaining ordinary bins.
    size_t length() const @property { return _end - _begin; }

    /// Whether all bins in this range have been consumed.
    bool empty() const @property { return _begin == _end; }

    /// First remaining element, returned by value.
    Element front() const @property
    {
        assert(!empty, "HistogramBinView.front: empty range");
        return this[0];
    }

    /// Last remaining element, returned by value.
    Element back() const @property
    {
        assert(!empty, "HistogramBinView.back: empty range");
        return this[length - 1];
    }

    /// Advance past the first remaining bin.
    void popFront()
    {
        assert(!empty, "HistogramBinView.popFront: empty range");
        ++_begin;
    }

    /// Remove the last remaining bin from this range.
    void popBack()
    {
        assert(!empty, "HistogramBinView.popBack: empty range");
        --_end;
    }

    /// Copy the traversal position, sharing the backing buffers.
    auto save() @property { return this; }

    /// Read an element relative to the current range.
    Element opIndex(size_t index) const
    {
        assert(index < length, "HistogramBinView: index is out of range");
        auto originalIndex = _begin + index;
        return Element(originalIndex, _axis.bin(originalIndex), _counts[originalIndex]);
    }

    /// Return a subrange; element indices still refer to the original histogram.
    auto opSlice(size_t begin, size_t end)
    {
        assert(begin <= end && end <= length,
            "HistogramBinView: slice is out of range");
        auto result = this;
        result._begin = _begin + begin;
        result._end = _begin + end;
        return result;
    }

    /// Copy the full remaining range.
    auto opSlice() { return this; }

    /// Support $ in index and slice expressions.
    size_t opDollar() const { return length; }
}

/// Iterate over numeric bins alongside their counts.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    auto h = HistogramAccumulator!(uint[], Axis)([0u, 0u, 0u], Axis(3, 0.0));
    h.put([0.5, 1.0, 1.5, 2.5]);

    uint total;
    foreach (entry; h.bins)
    {
        assert(entry.bin.low == entry.index);
        assert(entry.bin.high == entry.index + 1);
        total += entry.count;
    }
    assert(total == 4);
}

/// Index and slice a view without losing the original bin indices.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    auto h = HistogramAccumulator!(uint[], Axis)([1u, 2u, 1u], Axis(3, 0.0));
    auto bins = h.bins;
    assert(bins[1].count == 2);
    assert(bins[1].bin.low == 1.0);

    auto middle = bins[1 .. $];
    assert(middle.length == 2);
    assert(middle.front.index == 1);
    assert(middle.back.index == 2);
}

/// Category bins expose a slot instead of interval boundaries.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, CategoryAxis;

    enum Label { first, second }
    alias Axis = CategoryAxis!(uint, Label, AxisOptions());
    auto h = HistogramAccumulator!(uint[], Axis)([0u, 0u], Axis());
    h.put([Label.first, Label.second, Label.second]);

    auto bins = h.bins;
    assert(bins[0].bin.slot == Label.first);
    assert(bins[0].count == 1);
    assert(bins[1].bin.slot == Label.second);
    assert(bins[1].count == 2);
}

/// Saved views have independent positions and share subsequent count updates.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    auto counts = rcslice!uint([1u, 2u]);
    auto h = HistogramAccumulator!(typeof(counts), Axis)(counts, Axis(2, 0.0));
    auto bins = h.bins;
    auto saved = bins.save;
    auto previous = bins.front;
    bins.popFront();
    h.put(0.5);

    assert(bins.front.index == 1);
    assert(saved.front.index == 0);
    assert(saved.front.count == 2);
    assert(previous.count == 1);

    previous.count = 100;
    assert(h.counts[0] == 2);
}

// Range semantics and sharing for both supported owning storage forms.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;
    import std.range.primitives: isRandomAccessRange, hasLength, hasSlicing,
        hasAssignableElements, isInfinite;

    alias Axis = IntegralAxis!(uint, double, AxisOptions(false, true, true));
    void check(Storage)(Storage counts, Storage replacement)
    {
        auto h = HistogramAccumulator!(Storage, Axis)(counts, Axis(3, 0.0));
        auto bins = h.bins;
        alias View = typeof(bins);
        static assert(isRandomAccessRange!View);
        static assert(hasLength!View && hasSlicing!View);
        static assert(!hasAssignableElements!View && !isInfinite!View);
        static assert(is(typeof(bins.front.count) == uint));
        bins[0].count = 10u; // Assigning to a returned temporary cannot update storage.
        assert(h.counts[0] == 0);

        h.put([-1.0, 0.5, 2.5, 4.0]);
        assert(h.underflow == 1 && h.overflow == 1);
        assert(bins.length == 3);
        assert(bins[0].count == 1 && bins[1].count == 0 && bins[2].count == 1);

        auto copy = bins.save;
        bins.popFront();
        bins.popBack();
        assert(bins.length == 1 && bins.front.index == 1);
        assert(bins.front == bins.back);
        auto sub = copy[1 .. 3][1 .. 2];
        assert(sub.front.index == 2);
        assert(copy[].length == 3);
        assert(copy[3 .. 3].empty);
        bins.popFront();
        assert(bins.empty && bins.length == 0);

        h.counts = replacement;
        h.axis[0] = Axis(3, 10.0);
        h.put(10.5);
        assert(copy.front.count == 1 && copy.front.bin.low == 0.0);
        assert(h.bins.front.count == 6 && h.bins.front.bin.low == 10.0);
    }
    check([0u, 0u, 0u], [5u, 0u, 0u]);
    check(rcslice!uint([0u, 0u, 0u]), rcslice!uint([5u, 0u, 0u]));
}

// All built-in axes preserve their existing bin descriptions.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.slice: sliced;
    import mir.ndslice.allocation: rcslice;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis,
        RegularAxis, TransformAxis, EnumAxis, CategoryAxis, VariableAxis;
    import mir.math.common: approxEqual;

    void checkNumeric(Axis)(Axis axis)
    {
        uint[] counts = [3u, 0u, 7u];
        auto h = HistogramAccumulator!(uint[], Axis)(counts, axis);
        const expectedAxis = axis;
        auto bins = h.bins;
        foreach (i; 0 .. bins.length)
        {
            auto expected = expectedAxis.bin(i);
            assert(bins[i].bin.low.approxEqual(expected.low));
            assert(bins[i].bin.high.approxEqual(expected.high));
            assert(bins[i].count == counts[i]);
        }
    }
    checkNumeric(IntegralAxis!(uint, double, AxisOptions())(3, 0.0));
    checkNumeric(RegularAxis!(uint, double, AxisOptions(true))(3, 0.0, 6.0));
    checkNumeric(RegularAxis!(uint, double, AxisOptions(false, false, false, true))(3, 0.0, 6.0));
    checkNumeric(TransformAxis!(uint, double, "a * 2", "a / 2", AxisOptions())(3, 0.0, 6.0));

    auto breaks = [0.0, 1.0, 3.0, 6.0].sliced;
    checkNumeric(VariableAxis!(uint, double*, AxisOptions())(breaks));

    enum Label { first, second }
    alias Enum = EnumAxis!(uint, Label);
    alias Category = CategoryAxis!(uint, Label, AxisOptions());
    auto enums = HistogramAccumulator!(uint[], Enum)([2u, 4u], Enum()).bins;
    auto categories = HistogramAccumulator!(uint[], Category)([2u, 4u], Category()).bins;
    assert(enums[0].bin.slot == Label.first && enums[1].bin.slot == Label.second);
    assert(categories[0].bin.slot == Label.first && categories[1].count == 4);

    // The view, then the returned bin, retains reference-counted break storage.
    auto makeView()
    {
        auto ownedBreaks = rcslice!double([0.0, 2.0, 5.0]);
        auto axis = VariableAxis!(uint, typeof(ownedBreaks._iterator), AxisOptions())(ownedBreaks);
        auto counts = rcslice!uint([2u, 3u]);
        return HistogramAccumulator!(typeof(counts), typeof(axis))(counts, axis).bins;
    }
    auto owned = makeView();
    assert(owned[1].bin.low == 2.0 && owned[1].count == 3);
    auto description = owned[1].bin;
    owned = typeof(owned).init;
    assert(description.low == 2.0 && description.high == 5.0);
    static assert(!__traits(compiles, description.low = 100.0));
}

// Invalid access and incompatible input are rejected.
version(mir_stat_test_hist)
unittest
{
    import core.exception: AssertError;
    import std.exception: assertThrown;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    alias H = HistogramAccumulator!(uint[], Axis);
    assertThrown!AssertError(H([0u], Axis(2, 0.0)).bins);
    auto bins = H([0u, 0u], Axis(2, 0.0)).bins;
    assertThrown!AssertError(bins[2]);
    assertThrown!AssertError(bins[size_t.max]);
    assertThrown!AssertError(bins[0 .. 3]);
    assertThrown!AssertError(bins[1 .. 0]);
    auto empty = bins[0 .. 0];
    assertThrown!AssertError(empty.front);
    assertThrown!AssertError(empty.back);
    assertThrown!AssertError(empty.popFront());
    assertThrown!AssertError(empty.popBack());

    alias Multi = HistogramAccumulator!(size_t[][], Axis, Axis);
    static assert(!__traits(compiles, Multi.init.bins()));

    struct CountingAxis
    {
        alias CountType = uint;
        alias BinType = double;
        uint N_bin() const { return 2; }
        uint index(double x) const { return cast(uint) x; }
    }
    static assert(isAxis!CountingAxis);
    alias CountingOnly = HistogramAccumulator!(uint[], CountingAxis);
    static assert(!__traits(compiles, CountingOnly.init.bins()));
}
