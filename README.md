# ModularSqrt
Fast modular square root of primes and prime powers, including 2. Interface uses GMP bigints.

Most implementations use slower approaches for at least some cases, or only support odd primes.

Composites are not supported because that problem is a superset of large integer factorization, and a variable number of square roots must be returned. However, it would be simple to use, e.g., https://github.com/komrad36/EllipticCurveFactorization to first break the modulus into prime power factors, then use this library to separately compute the roots for each prime power, then use the Chinese Remainder Theorem (CRT) to combine them and reconstitute the overall roots.

Usage:

```cpp
bool SqrtModPrime(mpz_ptr ret, mpz_srcptr a, mpz_srcptr p);
```

```cpp
bool SqrtModPrimePower(mpz_ptr ret, mpz_ptr retMod, mpz_srcptr a, mpz_srcptr p, uint64_t k);
```
