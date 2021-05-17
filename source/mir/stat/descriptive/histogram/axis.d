/++
This module contains algorithms for histogram axes.

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

module mir.stat.descriptive.histogram.axis;

import mir.stat.descriptive.histogram.traits: DefaultCountType, isBreakFunction;
import mir.ndslice.slice: isSlice;
import mir.ndslice.traits: isContiguousVector;
import std.meta: NoDuplicates;
import std.traits: EnumMembers;

///
struct IsRightClosed
{
    bool isRightClosed = false;
}

///
struct EnableOverflow
{
    bool enableOverflow = false;
}

///
struct EnableUnderflow
{
    bool enableUnderflow = false;
}

///
struct IsCircular
{
    bool isCircular = false;
}

/++
Options to provide to an `AxisType`.

+/
struct AxisOptions
{

private:

    /++
    Axis breaks are assumed to be non-overlappying. If `isRightClosed` equals
    `false` (default), then calculations assume the axis is left-closed and
    right-opened, as in `[a, b)` or `a <= x < b`; otherwise, if `isRightClosed`
    equals `true`, then the calculations assumeaxis is left-open and
    right-closed, as in `(a, b]` or `a < x <= b`.
    +/
    IsRightClosed value_isRightClosed = IsRightClosed();

    ///
    EnableOverflow value_enableOverflow = EnableOverflow();

    ///
    EnableUnderflow value_enableUnderflow = EnableUnderflow();
    
    /++
    Breaks are assumed to not wrap-around by default. If `isCircular` equals
    `true`, then the axis is circular and will wrap around. For instance, if the
    breaks are `[a, b)` and `[b, c)` then a value of `x = c` will be placed 
    into the first break instead of overflow (assuming it is enabled). One 
    use-case of circular breaks is data in polar coordinates. 
    +/
    IsCircular value_isCircular = IsCircular();

public:

    ///
    @safe pure nothrow @nogc
    bool isRightClosed() const
    {
        return value_isRightClosed.isRightClosed;
    }

    ///
    @safe pure nothrow @nogc
    bool enableOverflow() const
    {
        return value_enableOverflow.enableOverflow;
    }

    ///
    @safe pure nothrow @nogc
    bool enableUnderflow() const
    {
        return value_enableUnderflow.enableUnderflow;
    }

    ///
    @safe pure nothrow @nogc
    bool isCircular() const
    {
        return value_isCircular.isCircular;
    }

    ///
    @safe pure nothrow @nogc
    this(bool x) {
        value_isRightClosed = IsRightClosed(x);
    }

    ///
    @safe pure nothrow @nogc
    this(bool x, bool y) {
        value_isRightClosed = IsRightClosed(x);
        value_enableOverflow = EnableOverflow(y);
    }

    ///
    @safe pure nothrow @nogc
    this(bool x, bool y, bool z) {
        value_isRightClosed = IsRightClosed(x);
        value_enableOverflow = EnableOverflow(y);
        value_enableUnderflow = EnableUnderflow(z);
    }

    ///
    @safe pure nothrow @nogc
    this(bool w, bool x, bool y, bool z) {
        value_isRightClosed = IsRightClosed(w);
        value_enableOverflow = EnableOverflow(x);
        value_enableUnderflow = EnableUnderflow(y);
        value_isCircular = IsCircular(z);
    }

    @safe pure nothrow @nogc
    void set(Arg)(Arg arg) {
        static if (is(Arg == IsRightClosed)) {
            this.value_isRightClosed = arg;
        } else static if (is(Arg == EnableOverflow)) {
            this.value_enableOverflow = arg;
        } else static if (is(Arg == EnableUnderflow)) {
            this.value_enableUnderflow = arg;
        } else static if (is(Arg == IsCircular)) {
            this.value_isCircular = arg;
        } else {
            static assert(0, "AxisOptions.set: option not supported");
        }
    }

    ///
    @safe pure nothrow @nogc
    void set(Arg)(bool value) {
        set(Arg(value));
    }

    @safe pure nothrow @nogc
    Arg get(Arg)() const
    {
        static if (is(Arg == IsRightClosed)) {
            return value_isRightClosed;
        } else static if (is(Arg == EnableOverflow)) {
            return value_enableOverflow;
        } else static if (is(Arg == EnableUnderflow)) {
            return value_enableUnderflow;
        } else static if (is(Arg == IsCircular)) {
            return value_isCircular;
        } else {
            static assert(0, "AxisOptions.get: option not supported");
        }
    }

    ///
    @safe pure nothrow @nogc
    this(Arg)(Arg arg) {
        set(arg);
    }

    ///
    @safe pure nothrow @nogc
    this(Args...)(Args args) {
        foreach (arg; args)
        {
            set(arg);
        }
    }
}


// Example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    AxisOptions x1 = AxisOptions(IsRightClosed(true));
    assert(x1.isRightClosed == true);
    AxisOptions x2 = AxisOptions(EnableOverflow(true));
    assert(x2.enableOverflow == true);
    AxisOptions x3 = AxisOptions(EnableUnderflow(true));
    assert(x3.enableUnderflow == true);
    AxisOptions x4 = AxisOptions(IsCircular(true));
    assert(x4.isCircular == true);

    AxisOptions x5 = AxisOptions(IsRightClosed(true), IsCircular(true));
    assert(x5.isRightClosed == true);
    assert(x5.isCircular == true);

    x5.set!EnableOverflow(true);
    x5.set!EnableUnderflow(true);
    assert(x5.get!EnableOverflow == EnableOverflow(true));
    assert(x5.get!EnableUnderflow == EnableUnderflow(true));
    
    AxisOptions x6 = AxisOptions(true, true, true, true);
    assert(x6.isRightClosed == true);
    assert(x6.enableOverflow == true);
    assert(x6.enableUnderflow == true);
    assert(x6.isCircular == true);
}

// Complete checks
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    AxisOptions x1 = AxisOptions(IsRightClosed(true));
    assert(x1.isRightClosed == true);
    assert(x1.enableOverflow == false);
    assert(x1.enableUnderflow == false);
    assert(x1.isCircular == false);
    AxisOptions x2 = AxisOptions(EnableOverflow(true));
    assert(x2.isRightClosed == false);
    assert(x2.enableOverflow == true);
    assert(x2.enableUnderflow == false);
    assert(x2.isCircular == false);
    AxisOptions x3 = AxisOptions(EnableUnderflow(true));
    assert(x3.isRightClosed == false);
    assert(x3.enableOverflow == false);
    assert(x3.enableUnderflow == true);
    assert(x3.isCircular == false);
    AxisOptions x4 = AxisOptions(IsCircular(true));
    assert(x4.isRightClosed == false);
    assert(x4.enableOverflow == false);
    assert(x4.enableUnderflow == false);
    assert(x4.isCircular == true);

    AxisOptions x5 = AxisOptions(IsRightClosed(true), IsCircular(true));
    assert(x5.isRightClosed == true);
    assert(x5.enableOverflow == false);
    assert(x5.enableUnderflow == false);
    assert(x5.isCircular == true);

    x5.set!EnableOverflow(true);
    x5.set!EnableUnderflow(true);
    assert(x5.get!EnableOverflow == EnableOverflow(true));
    assert(x5.get!EnableUnderflow == EnableUnderflow(true));
    
    AxisOptions x6 = AxisOptions(true);
    assert(x6.isRightClosed == true);
    assert(x6.enableOverflow == false);
    assert(x6.enableUnderflow == false);
    assert(x6.isCircular == false);
    AxisOptions x7 = AxisOptions(true, true);
    assert(x7.isRightClosed == true);
    assert(x7.enableOverflow == true);
    assert(x7.enableUnderflow == false);
    assert(x7.isCircular == false);
    AxisOptions x8 = AxisOptions(true, true, true);
    assert(x8.isRightClosed == true);
    assert(x8.enableOverflow == true);
    assert(x8.enableUnderflow == true);
    assert(x8.isCircular == false);
    AxisOptions x9 = AxisOptions(true, true, true, true);
    assert(x9.isRightClosed == true);
    assert(x9.enableOverflow == true);
    assert(x9.enableUnderflow == true);
    assert(x9.isCircular == true);
}

///
struct Bin(T)
    if (!is(T == enum) && !isSlice!T)
{
    ///
    T low;
    ///
    T high;
}

///
struct Bin(T)
    if (is(T == enum))
{
    ///
    T slot;
}

///
struct Bin(T)
    if (isContiguousVector!T)
{
    import mir.primitives: DeepElementType;
    import mir.ndslice.slice: Slice;

    private T _payload;

    ///
    DeepElementType!T low()() const
    {
        return _payload[0];
    }
    
    ///
    DeepElementType!T high()() const
    {
        return _payload[1];
    }
    
    this(Iterator)(Slice!Iterator x)
    {
        assert(x.length == 2);
        _payload = x;
    }
}

/++
Axis for an interval of integral values with unit steps.

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF AxisOptions),
    $(LREF RegularAxis),
    $(LREF TransformAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
struct IntegralAxis(CountT, BinT, AxisOptions axisOptions)
{
private:
    CountType _N_bin;
    BinType _low;

public:
    ///
    alias CountType = CountT;

    ///
    alias BinType = BinT;
    
    ///
    alias options = axisOptions;

    ///
    this(CountType N_bin, BinType low)
    {
        _N_bin = N_bin;
        _low = low;
    }

    ///
    CountType N_bin()() const
    {
        return _N_bin;
    }

    ///
    BinType low()() const
    {
        return _low;
    }

    ///
    BinType high()() const
    {
        return _low + cast(BinType) _N_bin;
    }

    ///
    bool isUnderflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x <= _low;
        } else {
            return x < _low;
        }
    }

    ///
    bool isOverflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x > high();
        } else {
            return x >= high();
        }
    }

    ///
    CountType index()(BinType x) const
    {
        import mir.stat.descriptive.histogram.traits: checkOverUnderFlow;

        checkOverUnderFlow!(BinType, axisOptions)(x, _low, high());

        import mir.math.common: floor;
        import std.traits: isIntegral;

        static if (!axisOptions.isRightClosed) {
            static if (axisOptions.isCircular) {
                if (x == high()) {
                    return cast(CountType) 0;
                }
            }
            // Include a specialization for integral types because the behavior
            // is simpler here.
            static if (isIntegral!BinType) {
                return cast(CountType) (x - _low);
            } else {
                return cast(CountType) (floor(x) - _low);
            }
        } else {
            static if (axisOptions.isCircular) {
                if (x == _low) {
                    return cast(CountType) (_N_bin - 1);
                }
            }
            BinType binValue = x - _low;
            // Include a specialization for integral types because the behavior
            // is simpler here.
            static if (isIntegral!BinType) {
                return cast(CountType) (binValue - 1);
            } else {
                CountType output = cast(CountType) floor(binValue);
                // If binValue equals the floor of the binValue, then it is on integer, adjust for closed
                if (binValue != output) {
                    return output;
                } else {
                    return output - 1;
                }
            }
        }
    }

    ///
    Bin!BinType bin()(size_t x) const
    {
        assert(x < N_bin, "IntegralAxis.bin: input must be less than N_bin");
        return Bin!(BinType)(_low + x, _low + x + 1);
    }
}

/// Example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions())(10, 2.0);
    assert(integralAxis.high == 12);

    assert(!integralAxis.isOverflow(5.0));
    assert(!integralAxis.isUnderflow(5.0));
    assert(integralAxis.isOverflow(13.0));
    assert(integralAxis.isUnderflow(1.0));
    
    assert(integralAxis.index(2.0) == 0);
    assert(integralAxis.index(2.5) == 0);
    assert(integralAxis.index(3.0) == 1);
    assert(integralAxis.index(11.5) == 9);
    
    assert(integralAxis.bin(0) == Bin!double(2.0, 3.0));
    assert(integralAxis.bin(1) == Bin!double(3.0, 4.0));
    assert(integralAxis.bin(9) == Bin!double(11.0, 12.0));
}

// With isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions(true))(10, 2.0);

    assert(integralAxis.index(2.5) == 0);
    assert(integralAxis.index(3.0) == 0);
    assert(integralAxis.index(3.5) == 1);
    assert(integralAxis.index(4.0) == 1);
    assert(integralAxis.index(4.5) == 2);
    assert(integralAxis.index(5.0) == 2);
    assert(integralAxis.index(5.5) == 3);
    assert(integralAxis.index(6.0) == 3);
    assert(integralAxis.index(12.0) == 9);
    assert(integralAxis.index(11.5) == 9);
}

// Some more tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, double, AxisOptions())(10, 2.0);
    
    assert(integralAxis.index(3.5) == 1);
    assert(integralAxis.index(4.0) == 2);
    assert(integralAxis.index(4.5) == 2);
    assert(integralAxis.index(5.0) == 3);
    assert(integralAxis.index(5.5) == 3);
    assert(integralAxis.index(6.0) == 4);
}

// integral test
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, int, AxisOptions())(10, 2);

    assert(integralAxis.index(2) == 0);
    assert(integralAxis.index(4) == 2);
    assert(integralAxis.index(5) == 3);
    assert(integralAxis.index(6) == 4);
}

// integral test, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, int, AxisOptions(true))(10, 2);

    assert(integralAxis.index(4) == 1);
    assert(integralAxis.index(5) == 2);
    assert(integralAxis.index(6) == 3);
    assert(integralAxis.index(12) == 9);
}

// integral test, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, int, AxisOptions(IsCircular(true)))(10, 2);
    
    assert(integralAxis.index(2) == 0);
    assert(integralAxis.index(5) == 3);
    assert(integralAxis.index(12) == 0);
}

// integral test, isRightClosed = true, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto integralAxis = IntegralAxis!(size_t, int, AxisOptions(IsRightClosed(true), IsCircular(true)))(10, 2);
    
    assert(integralAxis.index(2) == 9);
    assert(integralAxis.index(5) == 2);
    assert(integralAxis.index(12) == 9);
}

/++
Factory function to produce $(LREF IntegralAxis) object

Params:
    N_bin = number of bins
    low = value of smallest bin

See_also:
    $(LREF IntegralAxis)
+/
IntegralAxis!(CountType, BinType, axisOptions)
    integralAxis(CountType, BinType, AxisOptions axisOptions = AxisOptions())(CountType N_bin, BinType low)
{
    return IntegralAxis!(CountType, BinType, axisOptions)(N_bin, low);
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
+/
IntegralAxis!(DefaultCountType, BinType, axisOptions)
    integralAxis(BinType, AxisOptions axisOptions = AxisOptions())(DefaultCountType N_bin, BinType low)
{
    return .integralAxis!(DefaultCountType, BinType, axisOptions)(N_bin, low);
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template integralAxis(CountType, BinType, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
    +/
    IntegralAxis!(CountType, BinType, axisOptions)
        integralAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low)
    {
        return .integralAxis!(CountType, BinType, axisOptions)(cast(CountType) breakFunction(slice.lightScope), low);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template integralAxis(BinType, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
    +/
    IntegralAxis!(DefaultCountType, BinType, axisOptions)
        integralAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low)
    {
        import core.lifetime: move;
        return .integralAxis!(DefaultCountType, BinType, breakFunction, axisOptions)(slice.move, low);
    }
}

/++
Params:
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template integralAxis(alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;
    import mir.primitives: DeepElementType;

    /++
    Params:
        slice = slice
        low = value of smallest bin
    +/
    IntegralAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), axisOptions)
        integralAxis(Iterator, size_t N, SliceKind kind, BinType)(Slice!(Iterator, N, kind) slice, BinType low)
            if (is(BinType : DeepElementType!(Slice!(Iterator, N, kind))))
    {
        import core.lifetime: move;
        return .integralAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), breakFunction, axisOptions)(slice.move, low);
    }
}

/++
Axis for an interval of values with equal width steps.

See $(LREF TransformAxis) for an alternative axis that allows for monotonic
transformations.

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF TransformAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
struct RegularAxis(CountT, BinT, AxisOptions axisOptions)
{
    import mir.math.common: fmamath;

private:
    CountType _N_bin;
    BinType _low;
    BinType _high;

public:
    ///
    alias CountType = CountT;

    ///
    alias BinType = BinT;
    
    ///
    alias options = axisOptions;

    ///
    this(CountType N_bin, BinType low, BinType high)
    {
        assert(high > low, "RegularAxis.this: high must be greater than low");
        _N_bin = N_bin;
        _low = low;
        _high = high;
    }

    ///
    CountType N_bin()() const
    {
        return _N_bin;
    }

    ///
    BinType low()() const
    {
        return _low;
    }

    ///
    BinType high()() const
    {
        return _high;
    }

    ///
    bool isUnderflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x <= _low;
        } else {
            return x < _low;
        }
    }

    ///
    bool isOverflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x > _high;
        } else {
            return x >= _high;
        }
    }

    ///
    @fmamath BinType stepSize()() const
    {
        return (_high - _low) / (cast(BinType) _N_bin);
    }

    ///
    @fmamath BinType value()(BinType x) const
    {
        return (x - _low) / (_high - _low);
    }

    ///
    @fmamath CountType index()(BinType x) const
    {
        import mir.stat.descriptive.histogram.traits: checkOverUnderFlow;

        checkOverUnderFlow!(BinType, axisOptions)(x, _low, _high);

        import mir.math.common: floor;
        static if (!axisOptions.isRightClosed) {
            static if (axisOptions.isCircular) {
                if (x == _high) {
                    return cast(CountType) 0;
                }
            }
            return cast(CountType) floor(_N_bin * this.value(x));
        } else {
            static if (axisOptions.isCircular) {
                if (x == _low) {
                    return cast(CountType) (_N_bin - 1);
                }
            }
            BinType binValue = _N_bin * this.value(x);
            CountType output = cast(CountType) floor(binValue);
            // If binValue equals the floor of the binValue, then it is on integer, adjust for closed
            if (binValue != output) {
                return output;
            } else {
                return output - 1;
            }
        }
    }

    ///
    @fmamath Bin!BinType bin()(size_t x) const
    {
        assert(x < N_bin, "RegularAxis.bin: input must be less than N_bin");
        return Bin!(BinType)(_low + x * stepSize(), _low + x * stepSize() + stepSize());
    }
}

/// Basic tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions())(10, 2.0, 12.0);
    assert(regularAxis.low == 2);
    assert(regularAxis.high == 12);

    assert(!regularAxis.isOverflow(5.0));
    assert(!regularAxis.isUnderflow(5.0));
    assert(regularAxis.isOverflow(13.0));
    assert(regularAxis.isUnderflow(1.0));
    
    assert(regularAxis.index(2.0) == 0);
    assert(regularAxis.index(2.5) == 0);
    assert(regularAxis.index(3.0) == 1);
    assert(regularAxis.index(11.5) == 9);

    assert(regularAxis.bin(0) == Bin!double(2.0, 3.0));
    assert(regularAxis.bin(1) == Bin!double(3.0, 4.0));
    assert(regularAxis.bin(9) == Bin!double(11.0, 12.0));
}

// isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(true))(10, 2.0, 12.0);

    assert(regularAxis.index(2.5) == 0);
    assert(regularAxis.index(3.0) == 0);
    assert(regularAxis.index(3.5) == 1);
    assert(regularAxis.index(12.0) == 9);
}

// isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(IsCircular(true)))(10, 2.0, 12.0);

    assert(regularAxis.index(2.5) == 0);
    assert(regularAxis.index(3.0) == 1);
    assert(regularAxis.index(3.5) == 1);
    assert(regularAxis.index(12.0) == 0);
}

// isRightClosed = true, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(IsRightClosed(true), IsCircular(true)))(10, 2.0, 12.0);

    assert(regularAxis.index(2.0) == 9);
    assert(regularAxis.index(2.5) == 0);
    assert(regularAxis.index(3.0) == 0);
    assert(regularAxis.index(3.5) == 1);
    assert(regularAxis.index(12.0) == 9);
}

// Some more tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions())(10, 2.0, 12.0);
    assert(regularAxis.stepSize == 1);
    assert(regularAxis.value(7.0) == 0.5);

    assert(regularAxis.index(3.5) == 1);
    assert(regularAxis.index(4.0) == 2);
    assert(regularAxis.index(4.5) == 2);
    assert(regularAxis.index(5.0) == 3);
    assert(regularAxis.index(5.5) == 3);
    assert(regularAxis.index(6.0) == 4);
}

// Double N_bin
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions())(20, 2.0, 12.0);
    
    assert(regularAxis.index(2.0) == 0);
    assert(regularAxis.index(2.25) == 0);
    assert(regularAxis.index(2.5) == 1);
    assert(regularAxis.index(2.75) == 1);
    assert(regularAxis.index(3.0) == 2);
    assert(regularAxis.index(3.5) == 3);
    assert(regularAxis.index(4.0) == 4);
    assert(regularAxis.index(4.5) == 5);
    assert(regularAxis.index(5.0) == 6);
    assert(regularAxis.index(5.5) == 7);
    assert(regularAxis.index(6.0) == 8);
    assert(regularAxis.index(11.5) == 19);
    assert(regularAxis.index(11.75) == 19);
}

// Double N_bin, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    auto regularAxis = RegularAxis!(size_t, double, AxisOptions(true))(20, 2.0, 12.0);

    assert(regularAxis.index(2.25) == 0);
    assert(regularAxis.index(2.5) == 0);
    assert(regularAxis.index(2.75) == 1);
    assert(regularAxis.index(3.0) == 1);
    assert(regularAxis.index(3.5) == 2);
    assert(regularAxis.index(4.0) == 3);
    assert(regularAxis.index(4.5) == 4);
    assert(regularAxis.index(5.0) == 5);
    assert(regularAxis.index(5.5) == 6);
    assert(regularAxis.index(6.0) == 7);
    assert(regularAxis.index(11.5) == 18);
    assert(regularAxis.index(11.75) == 19);
    assert(regularAxis.index(12.0) == 19);
}

/++
Factory function to produce $(LREF RegularAxis) object

Params:
    N_bin = number of bins
    low = value of smallest bin
    high = value of the largest bin

See_also:
    $(LREF RegularAxis)
+/
RegularAxis!(CountType, BinType, axisOptions)
    regularAxis(CountType, BinType, AxisOptions axisOptions = AxisOptions())(CountType N_bin, BinType low, BinType high)
{
    return RegularAxis!(CountType, BinType, axisOptions)(N_bin, low, high);
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
+/
RegularAxis!(DefaultCountType, BinType, axisOptions)
    regularAxis(BinType, AxisOptions axisOptions = AxisOptions())(DefaultCountType N_bin, BinType low, BinType high)
{
    return .regularAxis!(DefaultCountType, BinType, axisOptions)(N_bin, low, high);
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template regularAxis(CountType, BinType, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    RegularAxis!(CountType, BinType, axisOptions)
        regularAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        return .regularAxis!(CountType, BinType, axisOptions)(cast(CountType) breakFunction(slice.lightScope), low, high);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template regularAxis(BinType, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    RegularAxis!(DefaultCountType, BinType, axisOptions)
        regularAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        import core.lifetime: move;
        return .regularAxis!(DefaultCountType, BinType, breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Params:
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template regularAxis(alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;
    import mir.primitives: DeepElementType;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    RegularAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), axisOptions)
        regularAxis(Iterator, size_t N, SliceKind kind, BinType)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
            if (is(BinType : DeepElementType!(Slice!(Iterator, N, kind))))
    {
        import core.lifetime: move;
        return .regularAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Axis for an interval of values with equal underlying steps that may be
transformed to provide fast, unequal steps.

A $(LREF RegularAxis) is equivalent to a $(LREF TransformAxis) with an identity
`transform` function.

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins
    transform = function to transform axis
    axisOptions = options

See_also:
    $(MATHREF common, log10),
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
struct TransformAxis(CountT, BinT, alias transform, alias inverseTransform, AxisOptions axisOptions)
{
    import mir.math.common: fmamath;
    import mir.functional: naryFun;

private:
    RegularAxis!(CountType, BinType, axisOptions) regularAxis = void;
    BinType _low;
    BinType _high;

    alias transformFunction = naryFun!transform;
    alias inverseTransformFunction = naryFun!inverseTransform;
    alias transformType = typeof(transformFunction(BinType.init));
    alias inverseTransformType = typeof(inverseTransformFunction(BinType.init));
    static assert (is(BinType == inverseTransformType), "the return type of inverseTransform must match BinType");

public:
    ///
    alias CountType = CountT;

    ///
    alias BinType = BinT;
    
    ///
    alias options = regularAxis.options;

    ///
    this(CountType N_bin, BinType low, BinType high)
    {
        assert(high > low, "TransformAxis.this: high must be greater than low");
        regularAxis = RegularAxis!(CountType, BinType, axisOptions)(N_bin, transformFunction(low), transformFunction(high));
        _low = low;
        _high = high;
    }

    ///
    CountType N_bin()() const
    {
        return regularAxis._N_bin;
    }

    ///
    BinType low()() const
    {
        return _low;
    }

    ///
    BinType high()() const
    {
        return _high;
    }

    ///
    transformType lowTransform()() const
    {
        return regularAxis._low;
    }

    ///
    transformType highTransform()() const
    {
        return regularAxis._high;
    }

    ///
    bool isUnderflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x <= _low;
        } else {
            return x < _low;
        }
    }

    ///
    bool isOverflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x > _high;
        } else {
            return x >= _high;
        }
    }

    ///
    @fmamath BinType stepSize()() const
    {
        return regularAxis.stepSize();
    }

    ///
    @fmamath BinType value()(BinType x) const
    {
        return regularAxis.value(transformFunction(x));
    }

    ///
    @fmamath CountType index()(BinType x) const
    {
        return regularAxis.index(transformFunction(x));
    }
/*
    ///
    @fmamath Bin!BinType bin()(size_t x) const
    {
        import mir.math.common: log, log10, log2, sqrt;

        Bin!BinType regularBin = regularAxis.bin(x);

        static if (__traits(isSame, transform, log)) {
            import mir.math.common: exp;
            regularBin.low = exp(regularBin.low);
            regularBin.high = exp(regularBin.high);
        } else static if (__traits(isSame, transform, log10)) {
            regularBin.low = (cast(BinType) 10) ^^ regularBin.low;
            regularBin.high = (cast(BinType) 10) ^^ regularBin.high;
        } else static if (__traits(isSame, transform, log2)) {
            regularBin.low = (cast(BinType) 2) ^^ regularBin.low;
            regularBin.high = (cast(BinType) 2) ^^ regularBin.high;
        } else static if (__traits(isSame, transform, sqrt)) {
            regularBin.low = regularBin.low ^^ (cast(BinType) 2);
            regularBin.high = regularBin.high ^^ (cast(BinType) 2);
        } else {
            static assert(0, "TransformAxis.bin: Built-in inverse function not supplied (only defined for when transform is log, log10, log2, and sqrt), use overload to provide inverse transform function");
        }

        return regularBin;
    }
*/
    ///
    @fmamath Bin!BinType bin(size_t x) const
    {
        import mir.math.common: log, log10, log2, sqrt;

        Bin!BinType regularBin = regularAxis.bin(x);

        regularBin.low = inverseTransformFunction(regularBin.low);
        regularBin.high = inverseTransformFunction(regularBin.high);

        return regularBin;
    }
}

/// Basic tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.low == 10.0 ^^ 2.0);
    assert(transformAxis.high == 10.0 ^^ 12.0);

    assert(!transformAxis.isOverflow(10.0 ^^ 5.0));
    assert(!transformAxis.isUnderflow(10.0 ^^ 5.0));
    assert(transformAxis.isOverflow(10.0 ^^ 13.0));
    assert(transformAxis.isUnderflow(10.0 ^^ 1.0));
    
    assert(transformAxis.index(10.0 ^^ 2.0) == 0);
    assert(transformAxis.index(10.0 ^^ 2.5) == 0);
    assert(transformAxis.index(10.0 ^^ 3.0) == 1);
    assert(transformAxis.index(10.0 ^^ 11.5) == 9);

    assert(transformAxis.bin(0) == Bin!double(10.0 ^^ 2.0, 10.0 ^^ 3.0));
    assert(transformAxis.bin(1) == Bin!double(10.0 ^^ 3.0, 10.0 ^^ 4.0));
    assert(transformAxis.bin(9) == Bin!double(10.0 ^^ 11.0, 10.0 ^^ 12.0));

    // Can also supply lambda
    auto transformAxis2 = TransformAxis!(size_t, double, a => log10(a), a => 10.0 ^^ a, AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis2.index(10.0 ^^ 3.0) == 1);

    // Or string lambda
    auto transformAxis3 = TransformAxis!(size_t, double, "log10(a)", "10.0 ^^ a", AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis3.index(10.0 ^^ 3.0) == 1);
}

// isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions(true))(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.index(10.0 ^^ 2.5) == 0);
    assert(transformAxis.index(10.0 ^^ 3.0) == 0);
    assert(transformAxis.index(10.0 ^^ 3.5) == 1);
    assert(transformAxis.index(10.0 ^^ 12.0) == 9);
}

// isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions(IsCircular(true)))(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.index(10.0 ^^ 2.5) == 0);
    assert(transformAxis.index(10.0 ^^ 3.0) == 1);
    assert(transformAxis.index(10.0 ^^ 3.5) == 1);
    assert(transformAxis.index(10.0 ^^ 12.0) == 0);
}

// isRightClosed = true, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions(IsRightClosed(true), IsCircular(true)))(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.index(10.0 ^^ 2.0) == 9);
    assert(transformAxis.index(10.0 ^^ 2.5) == 0);
    assert(transformAxis.index(10.0 ^^ 3.0) == 0);
    assert(transformAxis.index(10.0 ^^ 3.5) == 1);
    assert(transformAxis.index(10.0 ^^ 12.0) == 9);
}

// Some more tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis.stepSize == 1);
    assert(transformAxis.value(10.0 ^^ 7.0) == 0.5);

    assert(transformAxis.index(10.0 ^^ 3.5) == 1);
    assert(transformAxis.index(10.0 ^^ 4.0) == 2);
    assert(transformAxis.index(10.0 ^^ 4.5) == 2);
    assert(transformAxis.index(10.0 ^^ 5.0) == 3);
    assert(transformAxis.index(10.0 ^^ 5.5) == 3);
    assert(transformAxis.index(10.0 ^^ 6.0) == 4);
}

// Double N_bin
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions())(20, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.index(10.0 ^^ 2.0) == 0);
    assert(transformAxis.index(10.0 ^^ 2.25) == 0);
    assert(transformAxis.index(10.0 ^^ 2.5) == 1);
    assert(transformAxis.index(10.0 ^^ 2.75) == 1);
    assert(transformAxis.index(10.0 ^^ 3.0) == 2);
    assert(transformAxis.index(10.0 ^^ 3.5) == 3);
    assert(transformAxis.index(10.0 ^^ 4.0) == 4);
    assert(transformAxis.index(10.0 ^^ 4.5) == 5);
    assert(transformAxis.index(10.0 ^^ 5.0) == 6);
    assert(transformAxis.index(10.0 ^^ 5.5) == 7);
    assert(transformAxis.index(10.0 ^^ 6.0) == 8);
    assert(transformAxis.index(10.0 ^^ 11.5) == 19);
    assert(transformAxis.index(10.0 ^^ 11.75) == 19);
}

// Double N_bin
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;

    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions(true))(20, 10.0 ^^ 2.0, 10.0 ^^ 12.0);

    assert(transformAxis.index(10.0 ^^ 2.25) == 0);
    assert(transformAxis.index(10.0 ^^ 2.5) == 0);
    assert(transformAxis.index(10.0 ^^ 2.75) == 1);
    assert(transformAxis.index(10.0 ^^ 3.0) == 1);
    assert(transformAxis.index(10.0 ^^ 3.5) == 2);
    assert(transformAxis.index(10.0 ^^ 4.0) == 3);
    assert(transformAxis.index(10.0 ^^ 4.5) == 4);
    assert(transformAxis.index(10.0 ^^ 5.0) == 5);
    assert(transformAxis.index(10.0 ^^ 5.5) == 6);
    assert(transformAxis.index(10.0 ^^ 6.0) == 7);
    assert(transformAxis.index(10.0 ^^ 11.5) == 18);
    assert(transformAxis.index(10.0 ^^ 11.75) == 19);
    assert(transformAxis.index(10.0 ^^ 12.0) == 19);
}

// Check bin
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: log10;
    
    double inverseLog10(double x) {
        return 10.0 ^^ x;
    }

    auto transformAxis1 = TransformAxis!(size_t, double, log10, inverseLog10, AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis1.bin(0) == Bin!double(10.0 ^^ 2.0, 10.0 ^^ 3.0));

    // lambda function
    auto transformAxis2 = TransformAxis!(size_t, double, a => log10(a), a => (10.0 ^^ a), AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis2.bin(0) == Bin!double(10.0 ^^ 2.0, 10.0 ^^ 3.0));

    // string lambda
    auto transformAxis3 = TransformAxis!(size_t, double, "log10(a)", "(10.0 ^^ a)", AxisOptions())(10, 10.0 ^^ 2.0, 10.0 ^^ 12.0);
    assert(transformAxis3.bin(0) == Bin!double(10.0 ^^ 2.0, 10.0 ^^ 3.0));
}

// Additional bin tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: exp, log, log2, sqrt;

    auto transformAxis1 = TransformAxis!(size_t, double, log, exp, AxisOptions())(10, exp(2.0), exp(12.0));
    assert(transformAxis1.bin(0) == Bin!double(exp(2.0), exp(3.0)));

    auto transformAxis2 = TransformAxis!(size_t, double, log2, "2.0 ^^ a", AxisOptions())(10, 2.0 ^^ 2.0, 2.0 ^^ 12.0);
    assert(transformAxis2.bin(0) == Bin!double(4.0, 8.0));

    auto transformAxis3 = TransformAxis!(size_t, double, sqrt, (a => a ^^ 2.0), AxisOptions())(10, 4.0, 144.0);
    assert(transformAxis3.bin(0) == Bin!double(4.0, 9.0));
}

private T tenPow(T)(T x) {
    return 10 ^^ x;
}

private T twoPow(T)(T x) {
    return 2 ^^ x;
}

private T square(T)(T x) {
    return x ^^ 2;
}

/// TODO
template inverseTransformMapping(alias transform) {
    import mir.math.common: exp, log, log2, log10, sqrt;

    static if (__traits(isSame, transform, exp)) {
        alias inverseTransformMapping = log;
    } else static if (__traits(isSame, transform, log)) {
        alias inverseTransformMapping = exp;
    } else static if (__traits(isSame, transform, log2)) {
        alias inverseTransformMapping = twoPow;
    } else static if (__traits(isSame, transform, log10)) {
        alias inverseTransformMapping = tenPow;
    } else static if (__traits(isSame, transform, sqrt)) {
        alias inverseTransformMapping = square;
    } else {
        static assert (0, "inverseTransformMapping: transform does not match built-in support");
    }
}

/// TODO
template hasInverseTransformMapping(alias transform)
{
    import mir.math.common: exp, log, log2, log10, sqrt;

    static if (__traits(isSame, transform, exp)) {
        enum bool hasInverseTransformMapping = true;
    } else static if (__traits(isSame, transform, log)) {
        enum bool hasInverseTransformMapping = true;
    } else static if (__traits(isSame, transform, log2)) {
        enum bool hasInverseTransformMapping = true;
    } else static if (__traits(isSame, transform, log10)) {
        enum bool hasInverseTransformMapping = true;
    } else static if (__traits(isSame, transform, sqrt)) {
        enum bool hasInverseTransformMapping = true;
    } else {
        enum bool hasInverseTransformMapping = false;
    }
}

/++
Factory function to produce $(LREF TransformAxis) object

Params:
    N_bin = number of bins
    low = value of smallest bin
    high = value of the largest bin

See_also:
    $(LREF TransformAxis)
+/
TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)
    transformAxis(CountType, BinType, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())(CountType N_bin, BinType low, BinType high)
{
    return TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
}

/// TODO
TransformAxis!(CountType, BinType, transform, inverseTransformMapping!transform, axisOptions)
    transformAxis(CountType, BinType, alias transform, AxisOptions axisOptions = AxisOptions())(CountType N_bin, BinType low, BinType high)
        if (hasInverseTransformMapping!transform)
{
    alias inverseTransform = inverseTransformMapping!transform;
    return TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    transform = function to transform axis
    inverseTransform = function to undo transform
+/
TransformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)
    transformAxis(BinType, alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())(DefaultCountType N_bin, BinType low, BinType high)
{
    return .transformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
}

/// TODO
TransformAxis!(DefaultCountType, BinType, transform, inverseTransformMapping!transform, axisOptions)
    transformAxis(BinType, alias transform, AxisOptions axisOptions = AxisOptions())(DefaultCountType N_bin, BinType low, BinType high)
        if (hasInverseTransformMapping!transform)
{
    alias inverseTransform = inverseTransformMapping!transform;
    return .transformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
}

/++
Params:
    transform = function to transform axis
    inverseTransform = function to undo transform
+/
template transformAxis(alias transform, alias inverseTransform, AxisOptions axisOptions = AxisOptions())
{
    /++
    Params:
        N_bin = number of bins
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)
        transformAxis(BinType)(DefaultCountType N_bin, BinType low, BinType high)
    {
        return .transformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
    }
}

/++
Params:
    transform = function to transform axis
+/
template transformAxis(alias transform, AxisOptions axisOptions = AxisOptions())
    if (hasInverseTransformMapping!transform)
{
    /++
    Params:
        N_bin = number of bins
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, BinType, transform, inverseTransformMapping!transform, axisOptions)
        transformAxis(BinType)(DefaultCountType N_bin, BinType low, BinType high)
    {
        alias inverseTransform = inverseTransformMapping!transform;
        return .transformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)(N_bin, low, high);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    transform = function to transform axis
    inverseTransform = function to undo transform
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(CountType, BinType, alias transform, alias inverseTransform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        return .transformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)(cast(CountType) breakFunction(slice.lightScope), low, high);
    }
}

/++
Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    transform = function to transform axis
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(CountType, BinType, alias transform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (hasInverseTransformMapping!transform && isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(CountType, BinType, transform, inverseTransformMapping!transform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        alias inverseTransform = inverseTransformMapping!transform;
        return .transformAxis!(CountType, BinType, transform, inverseTransform, axisOptions)(cast(CountType) breakFunction(slice.lightScope), low, high);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    transform = function to transform axis
    inverseTransform = function to undo transform
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(BinType, alias transform, alias inverseTransform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, BinType, transform, inverseTransform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        import core.lifetime: move;
        return .transformAxis!(DefaultCountType, BinType, transform, inverseTransfornm, breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Params:
    BinType = the type of the values that are compared in histogram bins
    transform = function to transform axis
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(BinType, alias transform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (hasInverseTransformMapping!transform && isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, BinType, transform, inverseTransformMapping!transform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
    {
        import core.lifetime: move;
        alias inverseTransform = inverseTransformMapping!transform;
        return .transformAxis!(DefaultCountType, BinType, transform, inverseTransform, breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Params:
    transform = function to transform axis
    inverseTransform = function to undo transform
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(alias transform, alias inverseTransform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;
    import mir.primitives: DeepElementType;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), transform, inverseTransform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind, BinType)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
            if (is(BinType : DeepElementType!(Slice!(Iterator, N, kind))))
    {
        import core.lifetime: move;
        return .transformAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), transform, inverseTransform, breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Params:
    transform = function to transform axis
    breakFunction = function used to determine breaks
    axisOptions = options
+/
template transformAxis(alias transform, alias breakFunction, AxisOptions axisOptions = AxisOptions())
    if (hasInverseTransformMapping!transform && isBreakFunction!breakFunction)
{
    import mir.ndslice.slice: Slice, SliceKind;
    import mir.primitives: DeepElementType;

    /++
    Params:
        slice = slice
        low = value of smallest bin
        high = value of the largest bin
    +/
    TransformAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), transform, inverseTransformMapping!transform, axisOptions)
        transformAxis(Iterator, size_t N, SliceKind kind, BinType)(Slice!(Iterator, N, kind) slice, BinType low, BinType high)
            if (is(BinType : DeepElementType!(Slice!(Iterator, N, kind))))
    {
        import core.lifetime: move;
        
        alias inverseTransform = inverseTransformMapping!transform;
        return .transformAxis!(DefaultCountType, DeepElementType!(Slice!(Iterator, N, kind)), transform, inverseTransform, breakFunction, axisOptions)(slice.move, low, high);
    }
}

/++
Axis where the bins are made up of values from an enum.

Does not allow for overflow or underflow for improved performance.

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF TransformAxis),
    $(LREF CategoryAxis),
    $(LREF VariableAxis)
+/
struct EnumAxis(CountT, BinT)
    if (is(BinT == enum) &&
        EnumMembers!BinT.length == NoDuplicates!(EnumMembers!BinT).length)
{
    ///
    alias CountType = CountT;

    ///
    alias BinType = BinT;

    ///
    CountType N_bin()() const
    {
        import std.traits: EnumMembers;

        return EnumMembers!(BinType).length;
    }
    
    ///
    CountType index()(BinType value) const
    {
        import std.traits: OriginalType, EnumMembers;
        import mir.stat.descriptive.histogram.traits: isSwitchable;

        static if (isSwitchable!(OriginalType!BinType) && EnumMembers!BinType.length <= 50)
        {
            final switch (value)
            {
                foreach (size_t i, member; EnumMembers!BinType)
                {
                    case member:
                        return cast(CountType) i;
                }
            }
        }
        else
        {
            foreach (size_t i, member; EnumMembers!BinType)
            {
                if (value == member) {
                    break;
                }
                return cast(CountType) i;
            }
        }
    }
    
    ///
    Bin!BinType bin(size_t x)() const
    {
        import std.traits: EnumMembers;

        assert(x <= N_bin(), "EnumAxis.bin: input must be less than or equal to N_bin()");
        return Bin!(BinType)(EnumMembers!BinType[x]);
    }

    ///
    Bin!BinType bin()(size_t x) const
    {
        assert(x <= N_bin(), "EnumAxis.bin: input must be less than or equal to N_bin()");

        import std.traits: OriginalType, EnumMembers;
        import mir.stat.descriptive.histogram.traits: isSwitchable;

        static if (isSwitchable!(OriginalType!BinType) && EnumMembers!BinType.length <= 50)
        {
            final switch (x)
            {
                foreach (size_t i, member; EnumMembers!BinType)
                {
                    case i:
                        return Bin!BinType(member);
                }
            }
        }
        else
        {
            foreach (size_t i, member; EnumMembers!BinType)
            {
                if (x == i) {
                    break;
                }
                return Bin!BinType(member);
            }
        }
    }
}

///
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    enum Foo
    {
        A,
        B,
        C
    }
    EnumAxis!(size_t, Foo) enumAxis;

    assert(enumAxis.index(Foo.A) == 0);
    assert(enumAxis.index(Foo.B) == 1);
    assert(enumAxis.index(Foo.C) == 2);

    assert(enumAxis.bin(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin(2) == Bin!Foo(Foo.C));

    assert(enumAxis.bin!(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin!(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin!(2) == Bin!Foo(Foo.C));
}

///
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    enum Foo : string
    {
        A = "Z",
        B = "Y",
        C = "X"
    }
    EnumAxis!(size_t, Foo) enumAxis;

    assert(enumAxis.index(Foo.A) == 0);
    assert(enumAxis.index(Foo.B) == 1);
    assert(enumAxis.index(Foo.C) == 2);

    assert(enumAxis.bin(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin(2) == Bin!Foo(Foo.C));

    assert(enumAxis.bin!(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin!(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin!(2) == Bin!Foo(Foo.C));
}

///
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    enum Foo
    {
        A = 0,
        B = 1,
        C = 3
    }
    EnumAxis!(size_t, Foo) enumAxis;

    assert(enumAxis.index(Foo.A) == 0);
    assert(enumAxis.index(Foo.B) == 1);
    assert(enumAxis.index(Foo.C) == 2);

    assert(enumAxis.bin(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin(2) == Bin!Foo(Foo.C));

    assert(enumAxis.bin!(0) == Bin!Foo(Foo.A));
    assert(enumAxis.bin!(1) == Bin!Foo(Foo.B));
    assert(enumAxis.bin!(2) == Bin!Foo(Foo.C));
}

/++
Factory function to produce $(LREF EnumAxis)

Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins

See_also:
    $(LREF EnumAxis)
+/
EnumAxis!(CountType, BinType) enumAxis(CountType, BinType)()
{
    return EnumAxis!(CountType, BinType)();
}

/// ditto
EnumAxis!(DefaultCountType, BinType) enumAxis(BinType)()
{
    return .enumAxis!(DefaultCountType, BinType)();
}

/++
Axis similar to EnumAxis, but allows for overflow for when a string is passed
that does not match with enum members of `BinT`. 

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF TransformAxis),
    $(LREF EnumAxis),
    $(LREF VariableAxis)
+/
struct CategoryAxis(CountT, BinT, AxisOptions axisOptions)
    if (is(BinT == enum) &&
        EnumMembers!BinT.length == NoDuplicates!(EnumMembers!BinT).length)
{
    import std.traits: isSomeString;

    ///
    EnumAxis!(CountT, BinT) enumAxis;

    ///
    alias CountType = CountT;

    ///
    alias BinType = BinT;

    ///
    alias options = axisOptions;

    ///
    CountType N_bin()() const
    {
        return enumAxis.N_bin;
    }

    ///
    CountType index()(BinType value) const
    {
        return enumAxis.index(value);
    }

    ///
    CountType index(A)(A value) const
        if (isSomeString!A)
    {
        import mir.conv: to;
        try
        {
            return this.index(value.to!BinType);
        }
        catch (Exception e)
        {
            assert(0, "CategoryAxis.index: string value does not convert to enum, index is invalid here, increment overflow instead");   
        }
    }

    ///
    Bin!BinType bin(size_t x)() const
    {
        return enumAxis.bin!x;
    }

    ///
    Bin!BinType bin()(size_t x) const
    {
        return enumAxis.bin(x);
    }

    ///
    bool isOverflow()(BinType value) const
    {
        return false;    
    }

    ///
    bool isOverflow(A : const(char)[])(A value) const
    {
        import mir.conv: to;
        try
        {
            BinType x = value.to!BinType;
            return false;
        }
        catch (Exception e)
        {
            return true;
        }
    }

    ///
    bool isOverflow(A : const(char))(A value) const
    {
        import mir.conv: to;
        try
        {
            BinType x = value.to!BinType;
            return false;
        }
        catch (Exception e)
        {
            return true;
        }
    }
}

///
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    enum Foo
    {
        A,
        B,
        C
    }
    CategoryAxis!(size_t, Foo, AxisOptions()) categoryAxis;

    assert(categoryAxis.index(Foo.A) == 0);
    assert(categoryAxis.index(Foo.B) == 1);
    assert(categoryAxis.index(Foo.C) == 2);
 
    assert(categoryAxis.index("A") == 0);
    assert(categoryAxis.index("B") == 1);
    assert(categoryAxis.index("C") == 2);
    
    assert(categoryAxis.N_bin == 3);

    assert(categoryAxis.bin(0) == Bin!Foo(Foo.A));
    assert(categoryAxis.bin(1) == Bin!Foo(Foo.B));
    assert(categoryAxis.bin(2) == Bin!Foo(Foo.C));

    assert(categoryAxis.bin!(0) == Bin!Foo(Foo.A));
    assert(categoryAxis.bin!(1) == Bin!Foo(Foo.B));
    assert(categoryAxis.bin!(2) == Bin!Foo(Foo.C));

    assert(!categoryAxis.isOverflow(Foo.B));
	assert(!categoryAxis.isOverflow("B"));
	assert(categoryAxis.isOverflow("D"));
}

// Check that assert thrown when string input does not match enum
version(mir_stat_test_hist)
unittest
{
    import core.exception: AssertError;
    import std.exception: assertThrown;

    enum Foo
    {
        A,
        B,
        C
    }
    CategoryAxis!(size_t, Foo, AxisOptions()) categoryAxis;
    
    assertThrown!AssertError(categoryAxis.index("D"));
}

/++
Factory function to produce $(LREF CategoryAxis)

Params:
    CountType = the type that is used to count in histogram bins
    BinType = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF CategoryAxis)
+/
CategoryAxis!(CountType, BinType, axisOptions)
    categoryAxis(CountType, BinType, AxisOptions axisOptions = AxisOptions())()
{
    return CategoryAxis!(CountType, BinType, axisOptions)();
}

/// ditto
CategoryAxis!(DefaultCountType, BinType, axisOptions)
    categoryAxis(BinType, AxisOptions axisOptions = AxisOptions())()
{
    return .categoryAxis!(DefaultCountType, BinType, axisOptions)();
}

/++
Axis for non-equidistant data.

Params:
    CountT = the type that is used to count in histogram bins
    BinT = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF TransformAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis)
+/
struct VariableAxis(CountT, Iterator, AxisOptions axisOptions)
{
    import mir.primitives: DeepElementType;
    import mir.ndslice.slice: Slice, SliceKind;

private:
    Slice!(Iterator) _payload;

public:

    ///
    alias CountType = CountT;

    ///
    alias BinType = DeepElementType!(Slice!(Iterator));
    
    ///
    alias options = axisOptions;

    ///
    this(It, SliceKind kind)(Slice!(It, 1LU, kind) slice)
    {
        _payload = slice;
    }

    ///
    CountType N_bin()() const
    {
        return cast(CountType) _payload.length - 1;
    }

    ///
    BinType low()() const
    {
        return _payload[0];
    }

    ///
    BinType high()() const
    {
        return _payload[$ - 1];
    }

    ///
    bool isUnderflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x <= low();
        } else {
            return x < low();
        }
    }

    ///
    bool isOverflow()(BinType x) const
    {
        static if (axisOptions.isRightClosed &&
                   !axisOptions.isCircular) {
            return x > high();
        } else {
            return x >= high();
        }
    }

    ///
    CountType index()(BinType x)
    {
        import mir.stat.descriptive.histogram.traits: checkOverUnderFlow;

        checkOverUnderFlow!(BinType, axisOptions)(x, low(), high());

        static if (!axisOptions.isRightClosed) {
            static if (axisOptions.isCircular) {
                if (x == high()) {
                    return cast(CountType) 0;
                }
            }
            import std.range: assumeSorted;
            return _payload.assumeSorted!("a <= b").lowerBound(x).length - 1;
        } else {
            static if (axisOptions.isCircular) {
                if (x == low()) {
                    return cast(CountType) (N_bin() - 1);
                }
            }
            import std.range: assumeSorted;
            return _payload.assumeSorted!("a < b").lowerBound(x).length - 1;
        }
    }

    ///
    @trusted Bin!(Slice!(Iterator)) bin()(size_t x)
    {
        assert(x < _payload.length - 1, "VariableAxis.bin: input must be less than the length of _payload minus one");
        return Bin!(Slice!(Iterator))(_payload.select!0(x, (x + 2)));
    }
}

/// Example
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(double)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = i + 2.0;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(double), AxisOptions())(counts.asSlice);
    assert(variableAxis.N_bin == (counts.length - 1));
    assert(variableAxis.low == 2.0);
    assert(variableAxis.high == 12.0);

    assert(!variableAxis.isOverflow(5.0));
    assert(!variableAxis.isUnderflow(5.0));
    assert(variableAxis.isOverflow(13.0));
    assert(variableAxis.isUnderflow(1.0));

    assert(variableAxis.index(2.0) == 0);
    assert(variableAxis.index(2.5) == 0);
    assert(variableAxis.index(3.0) == 1);
    assert(variableAxis.index(11.5) == 9);

    assert(variableAxis.bin(0).low == 2.0);
    assert(variableAxis.bin(0).high == 3.0);
    assert(variableAxis.bin(1).low == 3.0);
    assert(variableAxis.bin(1).high == 4.0);
    assert(variableAxis.bin(9).low == 11.0);
    assert(variableAxis.bin(9).high == 12.0);
}

// With isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(double)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = i + 2.0;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(double), AxisOptions(true))(counts.asSlice);

    assert(variableAxis.index(2.5) == 0);
    assert(variableAxis.index(3.0) == 0);
    assert(variableAxis.index(3.5) == 1);
    assert(variableAxis.index(4.0) == 1);
    assert(variableAxis.index(4.5) == 2);
    assert(variableAxis.index(5.0) == 2);
    assert(variableAxis.index(5.5) == 3);
    assert(variableAxis.index(6.0) == 3);
    assert(variableAxis.index(12.0) == 9);
    assert(variableAxis.index(11.5) == 9);
}

// Some more tests
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(double)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = i + 2.0;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(double), AxisOptions())(counts.asSlice);

    assert(variableAxis.index(3.5) == 1);
    assert(variableAxis.index(4.0) == 2);
    assert(variableAxis.index(4.5) == 2);
    assert(variableAxis.index(5.0) == 3);
    assert(variableAxis.index(5.5) == 3);
    assert(variableAxis.index(6.0) == 4);
}

// integral test
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(int)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = cast(int) i + 2;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(int), AxisOptions())(counts.asSlice);

    assert(variableAxis.index(2) == 0);
    assert(variableAxis.index(4) == 2);
    assert(variableAxis.index(5) == 3);
    assert(variableAxis.index(6) == 4);
}

// integral test, isRightClosed = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(int)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = cast(int) i + 2u;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(int), AxisOptions(true))(counts.asSlice);

    assert(variableAxis.index(4) == 1);
    assert(variableAxis.index(5) == 2);
    assert(variableAxis.index(6) == 3);
    assert(variableAxis.index(12) == 9);
}

// integral test, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(int)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = cast(int) i + 2u;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(int), AxisOptions(IsCircular(true)))(counts.asSlice);

    assert(variableAxis.index(2) == 0);
    assert(variableAxis.index(5) == 3);
    assert(variableAxis.index(12) == 0);
}

// integral test, isRightClosed = true, isCircular = true
version(mir_stat_test_hist)
@safe pure nothrow @nogc
unittest
{
    import mir.rc.array;

    size_t len = 11;
    auto counts = mininitRcarray!(int)(len);
    size_t i = 0;
    while (i < len)
    {
        counts[i] = cast(int) i + 2u;
        i++;
    }
    auto variableAxis = VariableAxis!(size_t, RCI!(int), AxisOptions(IsRightClosed(true), IsCircular(true)))(counts.asSlice);

    assert(variableAxis.index(2) == 9);
    assert(variableAxis.index(5) == 2);
    assert(variableAxis.index(12) == 9);
}

/++
Factory function to produce $(LREF VariableAxis) object

Params:
    CountT = the type that is used to count in histogram bins
    Iterator = the type of the values that are compared in histogram bins
    axisOptions = options

See_also:
    $(LREF VariableAxis)
+/
template variableAxis(CountType, Iterator, AxisOptions axisOptions = AxisOptions())
{
    import core.lifetime: move;
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
    +/
    VariableAxis!(CountType, Iterator, axisOptions)
        variableAxis(size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice)
    {
        return VariableAxis!(CountType, Iterator, axisOptions)(slice.move);
    }
}

/++
Params:
    Iterator = the type of the values that are compared in histogram bins
    axisOptions = options
+/
template variableAxis(Iterator, AxisOptions axisOptions = AxisOptions())
{
    import core.lifetime: move;
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
    +/
    VariableAxis!(DefaultCountType, Iterator, axisOptions)
        variableAxis(size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice)
    {
        return .variableAxis!(DefaultCountType, Iterator, axisOptions)(slice.move);
    }
}

/++
Params:
    axisOptions = options
+/
template variableAxis(AxisOptions axisOptions = AxisOptions())
{
    import core.lifetime: move;
    import mir.ndslice.slice: Slice, SliceKind;

    /++
    Params:
        slice = slice
    +/
    VariableAxis!(DefaultCountType, Iterator, axisOptions)
        variableAxis(Iterator, size_t N, SliceKind kind)(Slice!(Iterator, N, kind) slice)
    {
        return .variableAxis!(DefaultCountType, Iterator, axisOptions)(slice.move);
    }
}
