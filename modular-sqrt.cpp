/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Jun 3, 2021
*******************************************************************/

#include "modular-sqrt.h"

//#define ENABLE_ASSERTS

#ifdef ENABLE_ASSERTS
#define ASSERT(a) do { if (!(a)) { printf("ASSERTION FAILED: %s\n", #a); __debugbreak(); } } while (0)
#else
#define ASSERT(a) do {} while (0)
#endif

using U32 = uint32_t;
using U64 = uint64_t;

class MpzJanitor
{
    mpz_t& m_x;

public:
    MpzJanitor(mpz_t& x) : m_x(x)
    {
        mpz_init(x);
    }
    ~MpzJanitor()
    {
        mpz_clear(m_x);
    }
};

#define MPZ_CREATE(x) mpz_t x; MpzJanitor mpzJanitor##x((x))

// returns (not a deep copy) whichever mpz (x or p-x) is smaller
// t is a temporary scratch variable that can be the same as p
// but must not be the same as x
static inline mpz_srcptr MinWithNeg(mpz_ptr t, mpz_srcptr x, mpz_srcptr p)
{
    ASSERT(t != x);
    ASSERT(mpz_cmp_ui(x, 0) > 0);
    ASSERT(mpz_cmp(x, p) < 0);
    mpz_sub(t, p, x);
    return mpz_cmp(x, t) < 0 ? x : t;
}

bool SqrtModPrime(mpz_ptr ret, mpz_srcptr a, mpz_srcptr p)
{
    // ret = a (mod p)
    mpz_mod(ret, a, p);

    // sqrt(0) (mod p) == 0
    // sqrt(1) (mod p) == 1
    if (mpz_cmp_ui(ret, 1) <= 0)
        return true;

    // the above two cases are the only possible cases for p == 2,
    // so at this point, 1 < b < p and p is an odd prime, e.g. p == 1 or 3 (mod 4)

    // solution exists if and only if a is a quadratic residue (mod p), i.e. a is a perfect square (mod p)
    // because p is now known to be an odd prime, legendre symbol can be used
    // to quickly check this

    if (mpz_legendre(ret, p) != 1)
        return false;

    // solution exists

    // for sufficiently small p, perform direct search
    if (mpz_cmp_ui(p, 1300) < 0)
    {
        const U32 sP = U32(mpz_get_ui(p));
        const U32 sB = U32(mpz_get_ui(ret));
        for (U32 r = 2;; ++r)
        {
            if ((r * r) % sP == sB)
            {
                mpz_set_ui(ret, r);
                return true;
            }
        }
    }

    // b = a (mod p)
    mpz_ptr b = ret;

    MPZ_CREATE(t);

    // t = p (mod 4)
    mpz_mod_2exp(t, p, 2);
    if (mpz_cmp_ui(t, 3) == 0)
    {
        // case p == 3 (mod 4)
        // x = +/- b^((p+1)/4) (mod p)

        // t = p + 1
        mpz_add_ui(t, p, 1);
        // t = (p+1)/4
        mpz_div_2exp(t, t, 2);
        // t = b^((p+1)/4) (mod p)
        mpz_powm(t, b, t, p);

        mpz_set(ret, MinWithNeg(b, t, p));

        return true;
    }

    MPZ_CREATE(u);

    // p == 1 or 5 (mod 8)

    // t = p (mod 8)
    mpz_mod_2exp(t, p, 3);
    if (mpz_cmp_ui(t, 5) == 0)
    {
        // case p == 5 (mod 8)
        // v = (2b)^((p-5)/8)
        // x = +/- bv(2bv^2-1) (mod p)

        // u = 2b
        mpz_mul_2exp(u, b, 1);
        // t = (p-5)/8
        mpz_div_2exp(t, p, 3);
        // t = v = (2b)^((p-5)/8) (mod p)
        mpz_powm(t, u, t, p);
        // u = 2bv
        mpz_mul(u, u, t);
        // u = 2bv^2
        mpz_mul(u, u, t);
        // u = 2bv^2-1
        mpz_sub_ui(u, u, 1);
        // u = 2bv^2-1 (mod p)
        mpz_mod(u, u, p);
        // u = v(2bv^2-1)
        mpz_mul(u, u, t);
        // u = bv(2bv^2-1)
        mpz_mul(u, u, b);
        // u = bv(2bv^2-1) (mod p)
        mpz_mod(u, u, p);

        mpz_set(ret, MinWithNeg(b, u, p));

        return true;
    }

    // p == 1 (mod 8)

    mpz_sub_ui(t, p, 1);
    const U64 S = mpz_scan1(t, 0);
    const U64 m = mpz_sizeinbase(p, 2);

    // select between Cipolla or Tonelli-Shanks methods
    if (S * (S - 1) > (m << 3) + 20)
    {
        // Cipolla

        // Step 1: find t such that t^2 - b is not square (mod p)

        // t = 0
        mpz_set_ui(t, 0);

        do
        {
            // ++t
            mpz_add_ui(t, t, 1);
            // u = t^2
            mpz_mul(u, t, t);
            // u = t^2 - b
            mpz_sub(u, u, b);
        } while (mpz_legendre(u, p) == 1);

        if (mpz_cmp_ui(u, 0) == 0)
        {
            // t^2 - b == 0 (mod p)
            // t^2 == b (mod p)
            mpz_set(ret, t);
            return true;
        }

        // Step 2: Compute x = (t + Sqrt(t^2 - b))^((p+1)/2)
        // trickier than it seems because t^2 - b is not square,
        // so requires exponentiation by squaring, with effectively complex arithmetic
        // to separately track the constant and root terms

        // x = (t+Sqrt(u))^((p+1)/2)

        // because p == 1 (mod 8), p ends in 0b001. p+1 ends in 0b010. (p+1)/2 ends in 0b01.
        // it always has first bit set, so we can start with a power of 1 instead of 0,
        // and omit the first bit from the loop itself. therefore we can just use p
        // instead of (p+1)/2, and start from the second bit.

        // x = (t + Sqrt(u))^((p+1)/2)

        MPZ_CREATE(v);
        MPZ_CREATE(s);
        MPZ_CREATE(f);
        MPZ_CREATE(g);

        // current result x = t + s*Sqrt(u)
        mpz_set_ui(s, 1);

        // current squaring level = f + g*Sqrt(u)
        mpz_set(f, t);
        mpz_set_ui(g, 1);

        for (U64 i = 1; i < m; ++i)
        {
            if (mpz_tstbit(p, i))
            {
                // multiply the current squaring level into result:
                //
                // (t + s*Sqrt(u))*(f + g*Sqrt(u))
                // t*f + t*g*Sqrt(u) + s*f*Sqrt(u) + s*g*u
                // t*f + s*g*u + t*g*Sqrt(u) + s*f*Sqrt(u)
                // (t*f + s*g*u) + (t*g + s*f)*Sqrt(u)
                // t = t*f + s*g*u, s = t*g + s*f

                mpz_mul(b, s, g);
                mpz_mul(v, t, f);
                mpz_addmul(v, b, u);
                mpz_mul(b, t, g);
                mpz_mod(t, v, p);
                mpz_addmul(b, s, f);
                mpz_mod(s, b, p);
            }

            // up the squaring level
            //
            // (f + g*Sqrt(u))(f + g*Sqrt(u))
            // f^2 + 2*f*g*Sqrt(u) + g^2*u
            // (f^2 + g^2*u) + (2*f*g)*Sqrt(u)
            // f = f^2 + g^2*u, g = 2*f*g

            mpz_mul(b, g, g);
            mpz_mul(v, f, f);
            mpz_addmul(v, b, u);
            mpz_mul(g, f, g);
            mpz_mod(f, v, p);
            mpz_mul_2exp(g, g, 1);
            mpz_mod(g, g, p);
        }

        ASSERT(mpz_cmp_ui(s, 0) == 0);

        mpz_set(ret, MinWithNeg(b, t, p));

        return true;
    }
    else
    {
        // Tonelli-Shanks

        // Step 1: set e and an odd t such that p = (2^e)t + 1
        U64 e = mpz_scan1(p, 1);

        mpz_div_2exp(t, p, e);

        // Step 2: find u such that u is not square (mod p)

        mpz_set_ui(u, 0);

        do
        {
            // ++u
            mpz_add_ui(u, u, 1);
        } while (mpz_legendre(u, p) == 1);

        // Step 3: init

        mpz_powm(u, u, t, p);

        mpz_div_2exp(t, t, 1);

        mpz_powm(t, b, t, p);

        mpz_mul(b, b, t);
        mpz_mod(b, b, p);

        mpz_mul(t, t, b);
        mpz_mod(t, t, p);

        MPZ_CREATE(v);

        // Step 4: loop until t == 1

        while (mpz_cmp_ui(t, 1) != 0)
        {
            U64 i = 0;
            mpz_set(v, t);
            do
            {
                ++i;
                mpz_mul(v, v, v);
                mpz_mod(v, v, p);
            } while (mpz_cmp_ui(v, 1) != 0);

            mpz_set_ui(v, 0);
            mpz_setbit(v, e - i - 1);
            mpz_powm(v, u, v, p);

            e = i;

            mpz_mul(u, v, v);
            mpz_mod(u, u, p);

            mpz_mul(t, t, u);
            mpz_mod(t, t, p);

            mpz_mul(b, b, v);
            mpz_mod(b, b, p);
        }

        mpz_set(ret, MinWithNeg(t, b, p));
        return true;
    }
}

bool SqrtModPrimePower(mpz_ptr ret, mpz_ptr retMod, mpz_srcptr a, mpz_srcptr p, U64 k)
{
    ASSERT(k >= 1);

    if (k == 1)
    {
        const bool exists = SqrtModPrime(ret, a, p);
        if (exists)
            mpz_set(retMod, p);
        return exists;
    }

    if (mpz_cmp_ui(p, 2) == 0)
    {
        // q = p^k

        mpz_ptr b = ret;

        // b = a (mod q)
        mpz_mod_2exp(b, a, k);

        if (mpz_cmp_ui(b, 0) == 0)
        {
            mpz_set_ui(retMod, 0);
            mpz_setbit(retMod, (k + 1) >> 1);
            return true;
        }

        const U64 e = mpz_scan1(b, 0);

        // b == f*(2^e) where f is odd
        // looking for x such that x*x == b (mod 2^k)
        // looking for x such that x*x == f*(2^e) (mod 2^k)
        // looking for x such that x*x == f*(2^e) + n(2^k)
        // looking for x such that x*x == 2^e * (f + n(2^(k-e)))
        // if e is odd, the 2^e term cannot split evenly between the two x's, so there is no solution
        if (e & 1)
            return false;

        // divide out factor of 2^e
        mpz_div_2exp(b, b, e);
        k -= e;

        // b is now odd and q is a power of 2, so a solution exists if and only if b == 1 (mod 8)

        mpz_ptr t = retMod;

        mpz_mod_2exp(t, b, 3);
        if (mpz_cmp_ui(t, 1) != 0)
            return false;

        if (k < 4)
        {
            mpz_set_ui(ret, 0);
            mpz_setbit(ret, e >> 1);
            mpz_set_ui(retMod, 0);
            mpz_setbit(retMod, (e >> 1) + 1);
            return true;
        }

        if (mpz_cmp_ui(b, 1) == 0)
        {
            mpz_set_ui(ret, 0);
            mpz_setbit(ret, e >> 1);
            mpz_set_ui(retMod, 0);
            mpz_setbit(retMod, k + (e >> 1) - 1);
            return true;
        }

        // solution is x = 1 (mod 4), must be lifted up to (mod 2^k)

        // invert: x^-1 = 1 (mod 4)

        // now lift modulus using inverse square root:
        // y' = y(3-ay^2)/2

        MPZ_CREATE(y);
        mpz_set_ui(y, 1);

        for (U64 i = 2; i < k - 1;)
        {
            i = (i << 1) - 1;

            mpz_mul(t, y, y);
            mpz_mul(t, t, b);
            mpz_ui_sub(t, 3, t);
            mpz_mul(y, y, t);
            ASSERT(!mpz_tstbit(y, 0));
            mpz_div_2exp(y, y, 1);
            mpz_mod_2exp(y, y, i);
        }

        // recover square root from inverse square root by sqrt(x) = x/sqrt(x)
        mpz_mul(y, y, b);
        mpz_mod_2exp(y, y, k - 1);
        mpz_mul_2exp(y, y, e >> 1);

        mpz_set_ui(retMod, 0);
        mpz_setbit(retMod, k + (e >> 1) - 1);

        mpz_set(ret, MinWithNeg(b, y, retMod));

        return true;
    }

    // p > 2

    mpz_ptr b = ret;

    // q = p^k
    mpz_pow_ui(b, p, k);

    // b = a (mod q)
    mpz_mod(b, a, b);

    if (mpz_cmp_ui(b, 0) == 0)
    {
        mpz_set_ui(ret, 0);
        mpz_pow_ui(retMod, p, (k + 1) >> 1);
        return true;
    }

    // b == f*(p^e) where f is odd
    // looking for x such that x*x == b (mod p^k)
    // looking for x such that x*x == f*(p^e) (mod p^k)
    // looking for x such that x*x == f*(p^e) + n(p^k)
    // looking for x such that x*x == p^e * (f + n(p^(k-e)))
    // if e is odd, the p^e term cannot split evenly between the two x's, so there is no solution
    const U64 e = mpz_remove(b, b, p);
    if (e & 1)
        return false;

    // divide out factor of p^e

    k -= e;

    mpz_ptr y = retMod;

    if (!SqrtModPrime(y, b, p))
        return false;

    // solution is now known (mod p), must be lifted up to (mod p^k)

    // invert
    {
        const bool inverseExists = mpz_invert(y, y, p);
        ASSERT(inverseExists);
        static_cast<void>(inverseExists);
    }

    // now lift modulus using inverse square root:
    // y' = y(3-ay^2)/2

    MPZ_CREATE(m);
    MPZ_CREATE(f);

    mpz_set(m, p);

    for (U64 i = 1; i < k;)
    {
        i <<= 1;
        mpz_mul(m, m, m);

        mpz_mul(f, y, y);
        mpz_mul(f, f, b);
        mpz_ui_sub(f, 3, f);
        mpz_mul(y, y, f);

        if (mpz_tstbit(y, 0))
            mpz_add(y, y, m);

        mpz_div_2exp(y, y, 1);
        mpz_mod(y, y, m);
    }

    // recover square root from inverse square root by sqrt(x) = x/sqrt(x)
    mpz_mul(y, y, b);
    mpz_pow_ui(m, p, k);
    mpz_mod(y, y, m);

    mpz_pow_ui(f, p, e >> 1);

    mpz_mul(ret, f, MinWithNeg(b, y, m));
    mpz_mul(retMod, f, m);

    return true;
}

#ifdef ALLOW_MODULAR_SQRT_OF_COMPOSITES

void SqrtModComposite::Init()
{
    ASSERT(mpz_cmp_ui(m_n, 0) > 0);

    mpz_init(m_sol);
    mpz_init(m_b);
    mpz_init(m_t);

    m_done = false;
    m_partialSolsNeedCleanup = true;

    // compute partial sols m_s[0] and m_s[1] of sqrt(a) (both modulo m_n) for each nontrivial factor of n
    // if at any point these do not exist, the overall solution does not exist either, and we return (empty) sols.
    m_partialSols = std::vector<PartialSol>(m_facN.size(), PartialSol());

    for (U64 i = 0; i < m_facN.size(); ++i)
    {
        const FactorInfo& info = m_facN[i];
        if (info.m_exp == 0)
            continue;

        PartialSol& partialSol = m_partialSols[i];

        mpz_init(partialSol.m_s[0]);
        mpz_init(partialSol.m_s[1]);
        mpz_init(partialSol.m_n);

        if (!SqrtModPrimePower(partialSol.m_s[0], partialSol.m_n, m_a, info.m_factor, info.m_exp))
        {
            // no solutions exist

            for (U64 j = 0; j <= i; ++j)
            {
                if (m_facN[j].m_exp)
                {
                    mpz_clear(m_partialSols[j].m_s[0]);
                    mpz_clear(m_partialSols[j].m_s[1]);
                    mpz_clear(m_partialSols[j].m_n);
                }
            }

            m_done = true;
            m_partialSolsNeedCleanup = false;
            return;
        }

        if (mpz_cmp_ui(partialSol.m_s[0], 0) != 0)
            mpz_sub(partialSol.m_s[1], partialSol.m_n, partialSol.m_s[0]);
    }

    // prepare counters for generating all possible composite solutions using Chinese Remainder Theorem
    m_counters = std::vector<U64>(m_facN.size(), 0ULL);

    ComputeCurrentSol();
}

void SqrtModComposite::ComputeCurrentSol()
{
    // m_sol = sum of contributions, one from each factor, for all permutations of all valid partial solutions:
    // let a be the partial sol
    // b = n / q (i.e. the product of all moduli except this one)
    // then the contribution from this factor is a * b * b^-1 (mod q)

    mpz_set_ui(m_sol, 0);

    for (U64 i = 0; i < m_facN.size(); ++i)
    {
        if (m_facN[i].m_exp == 0)
            continue;

        mpz_pow_ui(m_t, m_facN[i].m_factor, m_facN[i].m_exp);
        mpz_divexact(m_b, m_n, m_t);

        const bool hasInverse = mpz_invert(m_t, m_b, m_t);
        ASSERT(hasInverse);
        static_cast<void>(hasInverse);
        mpz_mul(m_b, m_b, m_t);

        mpz_mul_ui(m_t, m_partialSols[i].m_n, m_counters[i] >> 1);
        mpz_add(m_t, m_t, m_partialSols[i].m_s[m_counters[i] & 1]);

        mpz_addmul(m_sol, m_t, m_b);
    }

    mpz_mod(m_sol, m_sol, m_n);
}

void SqrtModComposite::Advance()
{
    // if we have no counters (e.g. because n == 1), we're done.
    if (m_facN.empty())
    {
        m_done = true;
        return;
    }

    // move to next permutation
    for (U64 i = 0;;)
    {
        if (m_facN[i].m_exp)
        {
            // increment counter
            m_counters[i] += mpz_cmp(m_partialSols[i].m_s[0], m_partialSols[i].m_s[1]) == 0 ? 2 : 1;

            mpz_mul_ui(m_b, m_partialSols[i].m_n, m_counters[i] >> 1);
            mpz_pow_ui(m_t, m_facN[i].m_factor, m_facN[i].m_exp);

            // if not too large, we're done incrementing. break out and process permutation
            if (mpz_cmp(m_b, m_t) < 0)
                break;

            // otherwise, cycle counter back to 0
            m_counters[i] = 0;
        }

        // advance to next counter
        ++i;

        // if we're out of counters, we're done with all solutions.
        if (i >= m_facN.size())
        {
            m_done = true;
            return;
        }
    }

    ComputeCurrentSol();
}

SqrtModComposite::~SqrtModComposite()
{
    if (m_partialSolsNeedCleanup)
    {
        for (U64 i = 0; i < m_facN.size(); ++i)
        {
            if (m_facN[i].m_exp)
            {
                mpz_clear(m_partialSols[i].m_s[0]);
                mpz_clear(m_partialSols[i].m_s[1]);
                mpz_clear(m_partialSols[i].m_n);
            }
        }
    }

    for (FactorInfo& info : m_computedFacN)
        mpz_clear(info.m_factor);

    mpz_clear(m_sol);
    mpz_clear(m_b);
    mpz_clear(m_t);
}

#endif //ALLOW_MODULAR_SQRT_OF_COMPOSITES
