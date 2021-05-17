/++
This module contains algorithms for frequency statistics.

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


module mir.stat.descriptive.histogram.frequency;

import mir.stat.descriptive.histogram.traits: isAxis;

/++
Accumulator used to generate frequency.

If the `AxisType` has an `options` member, the histogram may optionally allow
for overflow and underflow members.

Params:
    CountType = the type that is used to count in bins
    FrequencyType = the type that is used to output frequency bins
    AxisType = the type of the axis used to create the frequency bins

See_also:
    $(LREF AxisOptions),
    $(LREF IntegralAxis),
    $(LREF RegularAxis),
    $(LREF EnumAxis),
    $(LREF CategoryAxis),
    $(LREF HistogramAccumulator)
+/
struct FrequencyAccumulator(CountType, FrequencyType, AxisType)
    if (isAxis!AxisType)
{
    import std.traits: hasMember, isIterable, isSomeString;
    import mir.primitives: hasShape;
    import mir.stat.descriptive.histogram.traits: includeOverflow, includeUnderflow,
        BinTypeOf, isValidIndexOfAxis, isCategoryAxis;

    ///
    HistogramAccumulator!(CountType, AxisType) histogramAccumulator;

    ///
    CountType count;

    ///
    this(CountType[] x, AxisType y)
    {
        histogramAccumulator = HistogramAccumulator(x, y);
    }

    ///
    void put(Range)(Range r)
        if (isIterable!Range && 
            !(isCategoryAxis!AxisType && isSomeString!Range))
    {
        
        static if (hasShape!Range)
        {
            import mir.primitives: elementCount;

            count += r.elementCount;
            histogramAccumulator.put(x);
        }
        else
        {
            foreach(x; r)
            {
                put(x);
            }
        }
    }

    ///
    void put(T)(T x)
        if (is(T == BinTypeOf!AxisType) || 
            (isCategoryAxis!AxisType && isSomeString!T))
    {
        count++;
        histogramAccumulator.put(x);
    }

    ///
    void put(FrequencyAccumulator!(CountType, FrequencyType, AxisType) f)
    {
        import mir.stat.descriptive.histogram.traits: hasAxisOptions;

        assert(axis == f.axis);
        counts[] += f.counts[];
        count += f.count;
        static if (hasAxisOptions!AxisType) {
            static if (AxisType.options.enableOverflow && hasMember!(typeof(f), "overflow")) {
                overflow += f.overflow;
            }
            static if (AxisType.options.enableUnderflow && hasMember!(typeof(f), "underflow")) {
                overflow += f.underflow;
            }
        }
    }
    
    //needs frequency function
}

// how to handle overflow/underflow in frequency?
