/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Mar 25, 2021
*******************************************************************/

#pragma once

#include <cstdint>
#include <gmp.h>

// computes square root of 'a' modulo any prime number 'p'
// i.e. Sqrt(a) (mod p)
//
// 'p' MUST be a prime number. not a composite, and not a prime power. 2 is allowed.
// this is not checked; it is the caller's responsibility to ensure.
//
// 'ret' MUST be a valid, initialized mpz_t.
//
// if Sqrt(a) (mod p) does not exist, false is returned and 'ret' is undefined
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
// if Sqrt(a) (mod p^k) does not exist, false is returned, and 'ret' and 'retMod' are undefined
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
