/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Jun 3, 2021
*******************************************************************/

#pragma once


// enabling this requires https://github.com/komrad36/EllipticCurveFactorization or equivalent,
// because modular square root modulo a composite integer 'n' requires factoring 'n'.
#define ALLOW_MODULAR_SQRT_OF_COMPOSITES


#include <cstdint>
#include <gmp.h>
#include <vector>

#ifdef ALLOW_MODULAR_SQRT_OF_COMPOSITES
#include "factorize.h"
#endif //ALLOW_MODULAR_SQRT_OF_COMPOSITES

// computes square root of 'a' modulo any prime number 'p'
// i.e. Sqrt(a) (mod p)
//
// 'p' MUST be a prime number. not a composite, and not a prime power. 2 is allowed.
// this is not checked; it is the caller's responsibility to ensure.
//
// 'ret' MUST be a valid, initialized mpz_t.
//
// if Sqrt(a) (mod p) does not exist, false is returned and 'ret' is undefined.
//
// if Sqrt(a) (mod p) does exist, true is returned and 'ret' is set to the SMALLEST valid solution.
//
// -ret (mod p) is the second valid solution, which the caller can compute if desired
// note that this second solution is the same as the first when a == 0 (ret == -ret == 0)
// or when p == 2 (ret == -ret). in all other cases, it is unique, and can be computed cheaply
// by the caller if desired, using -ret (mod p) == p - ret.
bool SqrtModPrime(mpz_ptr ret, mpz_srcptr a, mpz_srcptr p);


// computes square root 'a' modulo any prime power 'p^k'
// i.e. Sqrt(a) (mod p^k)
//
// 'p' MUST be a prime number. not a composite, and not a prime power. 2 is allowed.
// this is not checked; it is the caller's responsibility to ensure.
//
// 'ret' and 'retMod' MUST be valid, initialized mpz_ts.
//
// 'k' must be >= 1, as otherwise, the modulus p^k is not a prime power.
//
// if Sqrt(a) (mod p^k) does not exist, false is returned, and 'ret' and 'retMod' are undefined.
//
// if Sqrt(a) (mod p^k) does exist, true is returned.
// 'ret' is set to the SMALLEST valid solution
// 'retMod' is set to the modulus of the solution space, which may be smaller than the original modulus p^k.
//
// -ret (mod retMod) is the second valid solution, which the caller can compute if desired
// note that this second solution is the same as the first when a == 0 (ret == -ret == 0)
// or when retMod is even and ret == -ret == retMod / 2. in all other cases, it is unique, and can be computed cheaply
// by the caller if desired, using -ret (mod retMod) == retMod - ret.
//
bool SqrtModPrimePower(mpz_ptr ret, mpz_ptr retMod, mpz_srcptr a, mpz_srcptr p, uint64_t k);

#ifdef ALLOW_MODULAR_SQRT_OF_COMPOSITES

// class that computes square root 'a' modulo any positive integer 'n', even composites
// i.e. Sqrt(a) (mod n)
//
// because the number of solutions to such a problem can get very large very quickly,
// this class is provided, allowing iteration through the solution space, as opposed to
// returning a vector of ALL solutions.
//
// solutions are not traversed in any particular order.
//
// the caller may provide an existing factorization:
// SqrtModComposite(mpz_srcptr a, mpz_srcptr n, const std::vector<FactorInfo>& facN)
//
// or, if one is not provided, one will be computed automatically:
// SqrtModComposite(mpz_srcptr a, mpz_srcptr n)
//
// Usage:
//
// for (mpz_srcptr sol : SqrtModComposite(a, n))
// {
//     // do stuff with sol
// }
//
class SqrtModComposite
{
private:
    struct PartialSol
    {
        mpz_t m_s[2];
        mpz_t m_n;
    };

    // forward decls
    class FwdIterator;
    class FwdIteratorEndSentinel;

public:
    SqrtModComposite(mpz_srcptr a, mpz_srcptr n) : m_a(a), m_n(n), m_computedFacN(Factorize(n)), m_facN(m_computedFacN)
    {
        Init();
    }

    SqrtModComposite(mpz_srcptr a, mpz_srcptr n, const std::vector<FactorInfo>& facN) : m_a(a), m_n(n), m_facN(facN)
    {
        Init();
    }

    ~SqrtModComposite();

    FwdIterator begin()
    {
        return *this;
    }

    FwdIteratorEndSentinel end()
    {
        return FwdIteratorEndSentinel();
    }

private:
    void Init();
    void ComputeCurrentSol();
    void Advance();

    class FwdIteratorEndSentinel {};

    class FwdIterator
    {
        friend class SqrtModComposite;

        FwdIterator(SqrtModComposite& parent) : m_parent(parent) {}

    public:
        bool operator==(const SqrtModComposite::FwdIteratorEndSentinel&)
        {
            return m_parent.m_done;
        }

        bool operator!=(const SqrtModComposite::FwdIteratorEndSentinel&)
        {
            return !m_parent.m_done;
        }

        void operator++()
        {
            m_parent.Advance();
        }

        mpz_srcptr operator*() const
        {
            return m_parent.m_sol;
        }

    private:
        SqrtModComposite& m_parent;
    };

private:
    mpz_srcptr m_a;
    mpz_srcptr m_n;
    std::vector<FactorInfo> m_computedFacN;
    const std::vector<FactorInfo>& m_facN;
    std::vector<PartialSol> m_partialSols;
    std::vector<uint64_t> m_counters;
    mpz_t m_sol;
    mpz_t m_b;
    mpz_t m_t;
    bool m_done;
    bool m_partialSolsNeedCleanup;
};

#endif //ALLOW_MODULAR_SQRT_OF_COMPOSITES
