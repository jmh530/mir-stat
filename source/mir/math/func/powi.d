/++
License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

Authors: John Michael Hall

Copyright: 2020 Mir Stat Authors.
+/

module mir.math.func.powi;

import mir.internal.utility: isFloatingPoint;
import std.traits: isSigned, isUnsigned, Unqual;

package(mir)
Unqual!T powi(T, U)(T x, U i)
    if (isFloatingPoint!(Unqual!(T)) && isSigned!(Unqual!(U)))
{
    import mir.math.common: powi;
    return powi(cast(Unqual!(T)) x, cast(Unqual!(U)) i);
}

package(mir)
double powi(T, U)(T x, U i)
    if (!isFloatingPoint!(Unqual!(T)) && is(Unqual!(T) : double) && isSigned!(Unqual!(U)))
{
    import mir.math.common: powi;
    return powi(cast(double) x, cast(Unqual!(U)) i);
}

package(mir)
Unqual!T powi(T, U)(T x, U i)
    if (isFloatingPoint!(Unqual!(T)) && isUnsigned!(Unqual!(U)))
{
    assert(i < int.max, "powi: converting unsigned i to signed will result in overflow");

    return powi(cast(Unqual!(T)) x, cast(int) i);
}

package(mir)
double powi(T, U)(T x, U i)
    if (!isFloatingPoint!(Unqual!(T)) && is(Unqual!(T) : double) && isUnsigned!(Unqual!(U)))
{
    assert(i < int.max, "powi: converting unsigned i to signed will result in overflow");

    return powi(cast(double) x, cast(int) i);
}

package(mir)
Unqual!T powi(T)(T x, size_t i)
    if (!isFloatingPoint!(Unqual!(T)) && !is(Unqual!(T) : double))
{
    if (i == 0) {
        return cast(Unqual!T) 0;
    } else if (i == 1) {
        return cast(Unqual!T) x;
    } else {
        auto output = cast(Unqual!T) x;
        for (size_t j = 1; j < i; j++) {
            output *= x;
        }
        return output;
    }
}

version(mir_stat_test)
@safe pure nothrow
unittest
{
    cdouble x = 1.0 + 2.0i;
    assert(x.powi(0) == 0);
    assert(x.powi(1) == x);
}