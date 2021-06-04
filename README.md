# ModularSqrt
Fast modular square root, modulo ANY positive integer. Primes, including 2; prime powers; and composites are all supported with custom handling. Interface uses GMP bigints.

Most implementations use slower approaches for at least some cases, or only support odd primes, or do not support composites, or do not support prime powers.

Computing a square root modulo a composite is a superset of integer factorization, so to enable that functionality, https://github.com/komrad36/EllipticCurveFactorization is required as a dependency.

Usage:

```cpp
bool SqrtModPrime(mpz_ptr ret, mpz_srcptr a, mpz_srcptr p);
```

```cpp
bool SqrtModPrimePower(mpz_ptr ret, mpz_ptr retMod, mpz_srcptr a, mpz_srcptr p, uint64_t k);
```

```cpp
 for (mpz_srcptr sol : SqrtModComposite(a, n))
{
    // do stuff with sol
}
```
