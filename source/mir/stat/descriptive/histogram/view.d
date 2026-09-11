/++
Read-only views of histogram bins and counts.

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall

Copyright: 2026 Mir Stat Authors.
+/
module mir.stat.descriptive.histogram.view;

import mir.ndslice.slice: isSlice;
import mir.primitives: DeepElementType;
import mir.qualifier: lightConst;
import mir.stat.descriptive.histogram.traits: isAxis;
import std.traits: isDynamicArray, isNumeric, Unqual;

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
            const Storage storage;
            const Axis axis;
            auto readOnlyCounts = lightConst(storage);
            auto readOnlyAxis = lightConst(axis);
            auto description = (const typeof(readOnlyAxis)).init.bin(size_t.init);
            auto count = (const typeof(readOnlyCounts)).init[size_t.init];
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
of shared boundaries. They must also support mir.qualifier.lightConst: value-only
axes can be copied from const, while axes containing mutable references should
provide a lightConst accessor that preserves ownership and makes those references
read-only.

A const histogram can create a view with a mutable traversal position. A const
view can be indexed, saved, and sliced; save and slicing return independent
mutable cursors over the same read-only data. Updates through an existing mutable
histogram remain visible. Constness does not make the backing data immutable.

Params:
    Storage = dynamic array or one-dimensional Mir slice of numeric counts
    Axis = axis with const runtime bin-description access
+/
struct HistogramBinView(Storage, Axis)
    if (supportsBinView!(Storage, Axis))
{
    private alias ReadOnlyStorage = typeof(lightConst((const Storage).init));
    private alias ReadOnlyAxis = typeof(lightConst((const Axis).init));
    private ReadOnlyStorage _counts;
    private ReadOnlyAxis _axis;
    private size_t _begin;
    private size_t _end;

    /// Type of each returned value.
    alias Element = HistogramBin!(
        typeof((const ReadOnlyAxis).init.bin(size_t.init)),
        Unqual!(DeepElementType!ReadOnlyStorage));

    /// Construct a view with one count per ordinary bin.
    this(const Storage counts, const Axis axis)
    {
        assert(counts.length == axis.N_bin,
            "HistogramBinView: count length must match axis bin count");
        _counts = lightConst(counts);
        _axis = lightConst(axis);
        _begin = 0;
        _end = counts.length;
    }

    private this(ReadOnlyStorage counts, ReadOnlyAxis axis, size_t begin, size_t end)
    {
        _counts = counts;
        _axis = axis;
        _begin = begin;
        _end = end;
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
    auto save() const @property
    {
        return HistogramBinView(lightConst(_counts), lightConst(_axis), _begin, _end);
    }

    /// Read an element relative to the current range.
    Element opIndex(size_t index) const
    {
        assert(index < length, "HistogramBinView: index is out of range");
        auto originalIndex = _begin + index;
        return Element(originalIndex, _axis.bin(originalIndex), _counts[originalIndex]);
    }

    /// Return a subrange; element indices still refer to the original histogram.
    auto opSlice(size_t begin, size_t end) const
    {
        assert(begin <= end && end <= length,
            "HistogramBinView: slice is out of range");
        return HistogramBinView(lightConst(_counts), lightConst(_axis),
            _begin + begin, _begin + end);
    }

    /// Copy the full remaining range.
    auto opSlice() const { return save; }

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

/// Const access preserves live counts while allowing independent traversal.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    alias H = HistogramAccumulator!(uint[], Axis);
    auto h = H([1u, 2u, 3u], Axis(3, 0.0));

    // A reporting function needs only const access to the histogram.
    auto readBins(ref const H histogram) { return histogram.bins; }
    const fixed = readBins(h);

    // save copies the traversal position into a mutable cursor.
    // Advancing that cursor leaves the const view at the first bin.
    auto cursor = fixed.save;
    cursor.popFront();
    assert(cursor.front.index == 1);
    assert(fixed.front.index == 0);

    // The const view shares the count buffer; it does not freeze the data.
    // Adding 0.5 through h increments the first bin from 1 to 2.
    h.put(0.5);
    assert(fixed.front.count == 2);

    // Slicing creates another mutable cursor, here covering bins 1 and 2.
    // popBack removes bin 2 from this cursor's range, without changing
    // the histogram's bins or the range covered by fixed.
    auto subset = fixed[1 .. $];
    subset.popBack();
    assert(subset.length == 1);
    assert(fixed.length == 3);
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

// Const sources produce mutable cursors over read-only storage handles.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;
    import std.algorithm: map, equal, find;
    import std.range: retro, take;
    import std.range.primitives: isRandomAccessRange;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    void check(Storage)(Storage counts)
    {
        alias H = HistogramAccumulator!(Storage, Axis);
        auto h = H(counts, Axis(3, 0.0));
        auto readBins(ref const H source) { return source.bins; }
        const fixed = readBins(h);
        auto cursor = fixed.save;
        static assert(is(typeof(cursor) == typeof(h.bins())));
        static assert(isRandomAccessRange!(typeof(cursor)));
        static assert(!__traits(compiles, fixed.popFront()));
        static assert(!__traits(compiles, fixed.popBack()));
        static assert(!__traits(compiles, cursor._counts[0] = 10u));
        static assert(!__traits(compiles, fixed._counts[0] = 10u));

        auto previous = fixed.front;
        h.put(0.5);
        assert(fixed.front.count == 2 && previous.count == 1);
        previous.count = 100;
        assert(h.counts[0] == 2);

        assert(cursor.map!(e => e.count).equal([2u, 2u, 3u]));
        assert(cursor.retro.map!(e => e.count).equal([3u, 2u, 2u]));
        assert(cursor.take(2).map!(e => e.count).equal([2u, 2u]));
        assert(cursor.find!(e => e.count == 3).front.index == 2);
        cursor.popFront();
        const advanced = cursor;
        auto saved = advanced.save;
        auto full = advanced[];
        auto tail = advanced[1 .. $];
        assert(saved.front.index == 1 && full.front.index == 1);
        assert(tail.front.index == 2);
        saved.popFront();
        full.popBack();
        assert(saved.front.index == 2 && full.back.index == 1);
        assert(advanced.length == 2 && fixed.length == 3);
        assert(advanced[0 .. 0].empty);
    }
    check([1u, 2u, 3u]);
    check(rcslice!uint([1u, 2u, 3u]));
}

// Both count and break ownership survive a const source and saved/sliced views.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.ndslice.allocation: rcslice;
    import mir.ndslice.slice: Slice;
    import mir.rc.array: RCI;
    import mir.stat.descriptive.histogram.accumulator: HistogramAccumulator;
    import mir.stat.descriptive.histogram.axis: AxisOptions, VariableAxis, Bin;

    auto makeView()
    {
        auto breaks = rcslice!double([0.0, 1.0, 3.0, 6.0]);
        alias Axis = VariableAxis!(uint, RCI!double, AxisOptions());
        auto counts = rcslice!uint([1u, 2u, 3u]);
        const h = HistogramAccumulator!(typeof(counts), Axis)(counts, Axis(breaks));
        const fixed = h.bins;
        return fixed.save[1 .. $];
    }
    auto view = makeView();
    assert(view.front.count == 2 && view.back.count == 3);
    assert(view.front.bin.low == 1.0 && view.back.bin.high == 6.0);
    static assert(is(typeof(view._counts) == Slice!(RCI!(const uint))));
    static assert(is(typeof(view._axis) ==
        VariableAxis!(uint, RCI!(const double), AxisOptions())));
    auto bin = view.front.bin;
    static assert(is(typeof(bin) == Bin!(Slice!(RCI!(const double)))));
    view = typeof(view).init;
    assert(bin.low == 1.0 && bin.high == 3.0);
}

// Additional storage forms keep const data readable and traversal independent.
version(mir_stat_test_hist)
unittest
{
    import mir.ndslice.slice: Slice, SliceKind;
    import mir.stat.descriptive.histogram.axis: AxisOptions, IntegralAxis;
    import std.algorithm: map, equal;

    alias Axis = IntegralAxis!(uint, double, AxisOptions());
    uint[] backing = [1u, 99u, 2u, 99u, 3u, 99u];
    auto strided = Slice!(uint*, 1, SliceKind.universal)([3], [2], backing.ptr);
    const view = HistogramBinView!(typeof(strided), Axis)(strided, Axis(3, 0.0));
    assert(view.save.map!(e => e.count).equal([1u, 2u, 3u]));
    backing[2] = 4;
    assert(view[1].count == 4);

    const(uint)[] counts = [1u, 2u, 3u];
    const readOnly = HistogramBinView!(typeof(counts), Axis)(counts, Axis(3, 0.0));
    auto cursor = readOnly.save;
    assert(cursor.map!(e => e.count).equal([1u, 2u, 3u]));
    static assert(is(typeof(cursor.front.count) == uint));
}

// Custom axes must provide an ownership-preserving const conversion when needed.
version(mir_stat_test_hist)
@safe pure nothrow
unittest
{
    import mir.stat.descriptive.histogram.axis: Bin;

    static struct ReferenceAxis
    {
        alias CountType = uint;
        alias BinType = double;
        double[] breaks;
        uint N_bin() const { return 1; }
        uint index(double x) const { return 0; }
        Bin!double bin(size_t i) const { return Bin!double(breaks[0], breaks[1]); }
    }
    static assert(!supportsBinView!(uint[], ReferenceAxis));

    static struct ReadOnlyAxis
    {
        alias CountType = uint;
        alias BinType = double;
        const(double)[] breaks;
        auto lightConst() const @property { return ReadOnlyAxis(breaks); }
        uint N_bin() const { return 1; }
        uint index(double x) const { return 0; }
        Bin!double bin(size_t i) const { return Bin!double(breaks[0], breaks[1]); }
    }
    static assert(supportsBinView!(uint[], ReadOnlyAxis));
    const view = HistogramBinView!(uint[], ReadOnlyAxis)(
        [2u], ReadOnlyAxis([0.0, 1.0]));
    assert(view.save.front.bin.high == 1.0 && view.front.count == 2);
}
