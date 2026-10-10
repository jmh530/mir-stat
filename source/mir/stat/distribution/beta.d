/++
This module contains algorithms for the $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution).

An alternate parameterization of this distribution is provided in $(MREF mir,stat,distribution,beta_proportion).

License: $(HTTP www.apache.org/licenses/LICENSE-2.0, Apache-2.0)

Authors: John Michael Hall

Copyright: 2022-3 Mir Stat Authors.

+/

module mir.stat.distribution.beta;

import mir.internal.utility: isFloatingPoint;

/++
Computes the beta probability density function (PDF).

Params:
    x = value to evaluate PDF
    alpha = shape parameter #1
    beta = shape parameter #2

See_also:
    $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution)
+/
@safe pure nothrow @nogc
T betaPDF(T)(const T x, const T alpha, const T beta)
    if (isFloatingPoint!T)
    in (x >= 0, "x must be greater than or equal to 0")
    in (x <= 1, "x must be less than or equal to 1")
    in (alpha > 0, "alpha must be greater than zero")
    in (beta > 0, "beta must be greater than zero")
{
    import mir.math.common: pow;
    import std.mathspecial: betaFunc = beta;

    return pow(x, (alpha - 1)) * pow((1 - x), (beta - 1)) / betaFunc(alpha, beta);
}

///
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: approxEqual;

    assert(0.5.betaPDF(1, 1) == 1);
    assert(0.75.betaPDF(1, 2).approxEqual(0.5));
    assert(0.25.betaPDF(0.5, 4).approxEqual(0.9228516));
}

/++
Computes the beta cumulatve distribution function (CDF).

Params:
    x = value to evaluate CDF
    alpha = shape parameter #1
    beta = shape parameter #2

See_also:
    $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution)
+/
@safe pure nothrow @nogc
T betaCDF(T)(const T x, const T alpha, const T beta)
    if (isFloatingPoint!T)
    in (x >= 0, "x must be greater than or equal to 0")
    in (x <= 1, "x must be less than or equal to 1")
    in (alpha > 0, "alpha must be greater than zero")
    in (beta > 0, "beta must be greater than zero")
{
    import std.mathspecial: betaIncomplete;

    return betaIncomplete(alpha, beta, x);
}

///
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: approxEqual;

    assert(0.5.betaCDF(1, 1).approxEqual(0.5));
    assert(0.75.betaCDF(1, 2).approxEqual(0.9375));
    assert(0.25.betaCDF(0.5, 4).approxEqual(0.8588867));
}

/++
Computes the beta complementary cumulative distribution function (CCDF).

For small $(D x), forming $(D 1 - x) can discard significant digits. Simply
subtracting the CDF from one is also inaccurate when the remaining upper tail
is tiny. For $(D x < 1.0 / 16), $(D alpha <= 1), and $(D beta <= 128), a
fallback can instead add the probability between $(D x) and an interior
point to the upper tail starting at that point. Other parameter regimes use
Phobos's reflected incomplete-beta calculation and retain its accuracy limits.

Params:
    x = value to evaluate CCDF
    alpha = shape parameter #1
    beta = shape parameter #2

See_also:
    $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution)
+/
@safe pure nothrow @nogc
T betaCCDF(T)(const T x, const T alpha, const T beta)
    if (isFloatingPoint!T)
    in (x >= 0, "x must be greater than or equal to 0")
    in (x <= 1, "x must be less than or equal to 1")
    in (alpha > 0, "alpha must be greater than zero")
    in (beta > 0, "beta must be greater than zero")
{
    import std.mathspecial: betaIncomplete;

    if (x < T(0.0625) && alpha <= 1 && beta <= 128)
        return cast(T) betaSmallInputCCDF(x, alpha, beta);
    return betaIncomplete(beta, alpha, 1 - x);
}

///
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: approxEqual;

    assert(0.5.betaCCDF(1, 1).approxEqual(0.5));
    assert(0.75.betaCCDF(1, 2).approxEqual(0.0625));
    assert(0.25.betaCCDF(0.5, 4).approxEqual(0.1411133));
}

// Compute the area to the right of x as two pieces:
//
//     x -------- c ---------------- 1
//       interval      upper tail
//
// Q(x) = integral(x..c, density) + Q(c). Neither piece is obtained by
// subtracting a probability close to one. This helper is deliberately bounded
// to alpha <= 1 and beta <= 128; it is not a replacement incomplete-beta API.
private @safe pure nothrow @nogc
real betaSmallInputCCDF(real x, real a, real b)
{
    import std.math: exp, expm1, log, fabs;
    import std.mathspecial: betaIncomplete, gamma;

    if (x == 0)
        return 1;

    // Choose a power-of-two anchor so 1-c is exact. Requiring b*c <= 1
    // also keeps the interval's binomial-series coefficients manageable.
    // In this parameter range c is at least 1/128.
    real c = 0.5L;
    while (b * c > 1)
        c *= 0.5L;
    if (x >= c)
        return betaIncomplete(b, a, 1 - x);

    // With subnormal b the distribution's endpoint mass dominates. This
    // ratio also handles both shapes being subnormal without overflowing B.
    // For these small x, the omitted correction is below real precision.
    if (b < real.min_normal)
        return a / (a + b);

    // 1/B(a,b) = a * b/(a+b) * Gamma(1+a+b)/Gamma(1+a)/Gamma(1+b).
    // Shifted gamma arguments avoid Gamma(a) overflow for tiny a. The shape
    // bounds keep these gamma values representable even in double precision.
    const real gammaRatio = gamma(1 + a + b) / gamma(1 + b) / gamma(1 + a);
    const real logRatio = log(x) - log(c);

    // Expand (1-t)^(b-1) and integrate each term from x to c. Factor a out
    // of the entire result: intermediate series terms then stay representable
    // even when the final probability is subnormal.
    // The first term is (1-(x/c)^a)/a. expm1 preserves it when a is tiny;
    // its limit -log(x/c) avoids division of an already rounded subnormal.
    real sum = fabs(a * logRatio) < real.epsilon
        ? -logRatio : -expm1(a * logRatio) / a;
    real ratioPower = exp(a * logRatio);
    const real ratio = x / c;
    real coefficient = 1;
    foreach (n; 1 .. 10_000)
    {
        coefficient *= (n - b) * c / n;
        ratioPower *= ratio;
        const real term = coefficient / (a + n) * (1 - ratioPower);
        sum += term;
        if (fabs(term) > real.epsilon * fabs(sum))
            continue;

        real anchorPerA;
        if (a <= real.epsilon && b <= real.epsilon)
        {
            // Here c=1/2. The limiting tail is a/(a+b); the first-order
            // correction at one half vanishes. Avoid a fragile Phobos call.
            anchorPerA = 1 / (a + b);
        }
        else if (a < real.min_normal)
        {
            // Evaluate the anchor with a normal surrogate, then restore the
            // endpoint-mass ratio. Changing this tiny shape has negligible
            // higher-order effect, but avoids overflow inside Phobos.
            const real normalA = 16 * real.min_normal;
            anchorPerA = betaIncomplete(b, normalA, 1 - c) / normalA
                * ((normalA + b) / (a + b));
        }
        else
            anchorPerA = betaIncomplete(b, a, 1 - c) / a;

        const real result = a * (anchorPerA
            + exp(a * log(c)) * gammaRatio * (b / (a + b)) * sum);
        // The two pieces may round slightly outside the probability bounds.
        return result > 1 ? 1 : result < 0 ? 0 : result;
    }
    // Do not silently return an unconverged probability.
    return real.nan;
}

// Recover small-input upper tails without rounding away tiny probabilities.
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import std.meta: AliasSeq;
    import std.math: exp, expm1, log, log1p, nextUp, nextDown;
    import mir.math.common: approxEqual;

    static foreach (T; AliasSeq!(float, double, real))
    {{
        // For beta=1, Q(x) = 1-x^alpha, evaluated independently with expm1.
        foreach (a; [T(0.01), T.epsilon * T.epsilon])
        foreach (x; [T(1e-20L), nextDown(T(0.0625)), T(0.0625),
            nextUp(T(0.0625)), T(0.25), nextDown(T(1))])
        {
            const T expected = -expm1(a * log(x));
            const T value = betaCCDF(x, a, T(1));
            assert(value >= 0 && value <= 1);
            assert(approxEqual(value, expected, 512 * T.epsilon, T(0)));
        }
        // The result itself may be subnormal; allow its final rounding step.
        const T tiny = nextUp(T(0));
        const T expected = -expm1(tiny * log(T(1e-20L)));
        assert(approxEqual(betaCCDF(T(1e-20L), tiny, T(1)), expected,
            512 * T.epsilon, 2 * tiny));
        assert(approxEqual(betaCCDF(T(1e-20L), tiny, tiny), T(0.5),
            512 * T.epsilon, T(0)));
        assert(betaCCDF(T(0), T(0.01), T(3)) == 1);
        assert(betaCCDF(T(1), T(0.01), T(3)) == 0);

        // Independent 380-digit incomplete-beta references; decimal input
        // rounding in T is included in the error allowance.
        assert(approxEqual(betaCCDF(T(1e-20L), T(0.03), T(7.25)),
            T(0.7295273104045321893412536502921225024L),
            1024 * T.epsilon, T(0)));
        assert(approxEqual(betaCCDF(T(0.001), T(0.1), T(128)),
            T(0.1540456317238734585412937362254863490L),
            1024 * T.epsilon, T(0)));
        // For alpha=1, Q(x)=(1-x)^beta. Check the anchor and shape-bound
        // transitions, including the unchanged path outside the fallback.
        foreach (b; [nextDown(T(128)), T(128), nextUp(T(128))])
        foreach (x; [nextDown(T(0.0078125)), T(0.0078125), nextUp(T(0.0078125))])
            // Evaluate the identity in real: older Phobos float/double
            // log1p implementations are inaccurate near this boundary.
            assert(approxEqual(betaCCDF(x, T(1), b),
                cast(T) exp(cast(real) b * log1p(-cast(real) x)),
                1024 * T.epsilon, T(0)));
    }}
}

/++
Computes the beta inverse cumulative distribution function (InvCDF).

Params:
    p = value to evaluate InvCDF
    alpha = shape parameter #1
    beta = shape parameter #2

See_also:
    $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution)
+/
@safe pure nothrow @nogc
T betaInvCDF(T)(const T p, const T alpha, const T beta)
    if (isFloatingPoint!T)
    in (p >= 0, "p must be greater than or equal to 0")
    in (p <= 1, "p must be less than or equal to 1")
    in (alpha > 0, "alpha must be greater than zero")
    in (beta > 0, "beta must be greater than zero")
{
    import std.mathspecial: betaIncompleteInverse;

    return betaIncompleteInverse(alpha, beta, p);
}

///
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: approxEqual;

    assert(0.5.betaInvCDF(1, 1).approxEqual(0.5));
    assert(0.9375.betaInvCDF(1, 2).approxEqual(0.75));
    assert(0.8588867.betaInvCDF(0.5, 4).approxEqual(0.25));
}

/++
Computes the beta log probability density function (LPDF).

Params:
    x = value to evaluate LPDF
    alpha = shape parameter #1
    beta = shape parameter #2

See_also:
    $(LINK2 https://en.wikipedia.org/wiki/Beta_distribution, Beta Distribution)
+/
@safe pure nothrow @nogc
T betaLPDF(T)(const T x, const T alpha, const T beta)
    if (isFloatingPoint!T)
    in (x >= 0, "x must be greater than or equal to 0")
    in (x <= 1, "x must be less than or equal to 1")
    in (alpha > 0, "alpha must be greater than zero")
    in (beta > 0, "beta must be greater than zero")
{
    import mir.math.internal.log_beta: logBeta;
    import mir.math.internal.xlogy: xlogy, xlog1py;

    return xlogy(alpha - 1, x) + xlog1py(beta - 1, -x) - logBeta(alpha, beta);
}

///
version(mir_stat_test)
@safe pure nothrow @nogc
unittest
{
    import mir.math.common: approxEqual, log;

    assert(0.5.betaLPDF(1, 1).approxEqual(log(betaPDF(0.5, 1, 1))));
    assert(0.75.betaLPDF(1, 2).approxEqual(log(betaPDF(0.75, 1, 2))));
    assert(0.25.betaLPDF(0.5, 4).approxEqual(log(betaPDF(0.25, 0.5, 4))));
}
