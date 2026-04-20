use num_bigint::BigUint;

/// Arithmetic operations in a prime finite field **F_p**.
///
/// All operations take values and a prime modulus `p` and return a result
/// in the range `[0, p)`.  The modulus is passed explicitly so the struct
/// itself is stateless — there is no field element type; raw `BigUint` values
/// are used throughout.
pub struct FiniteField {}

impl FiniteField {
    /// Returns `(c + d) mod p`.
    pub fn add(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        (c + d) % p
    }

    /// Returns `(c * d) mod p`.
    pub fn mul(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        (c * d) % p
    }

    /// Returns the additive inverse of `c` in **F_p**, i.e. `p - c`.
    ///
    /// # Panics
    /// Panics if `c >= p` (the value must already be reduced).
    pub fn inverse_add(c: &BigUint, p: &BigUint) -> BigUint {
        assert!(c < p, "number: {} is bigger or equal than p: {}", c, p);
        p - c
    }

    /// Returns `(c - d) mod p`.
    ///
    /// Computed as `c + (-d)` where `-d` is the additive inverse of `d`.
    /// Requires `d < p`.
    pub fn sub(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        let d_inv = FiniteField::inverse_add(d, p);
        FiniteField::add(c, &d_inv, p)
    }

    /// Returns the multiplicative inverse of `c` in **F_p** using Fermat's
    /// Little Theorem: `c^(-1) ≡ c^(p-2) (mod p)`.
    ///
    /// **Requires `p` to be prime.**  If `p` is composite the result is
    /// incorrect; use the extended Euclidean algorithm instead.
    pub fn inverse_mul_prime(c: &BigUint, p: &BigUint) -> BigUint {
        c.modpow(&(p - BigUint::from(2u32)), p)
    }

    /// Returns `(c / d) mod p`, i.e. `c * d^(-1) mod p`.
    ///
    /// Requires `p` to be prime (delegates to [`inverse_mul_prime`]).
    ///
    /// [`inverse_mul_prime`]: FiniteField::inverse_mul_prime
    pub fn div(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint {
        let d_inv = FiniteField::inverse_mul_prime(d, p);
        FiniteField::mul(c, &d_inv, p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_above_p() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(3u32));
    }

    #[test]
    fn test_add_below_p() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(31u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(14u32));
    }

    #[test]
    fn test_mul_above_p() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::mul(&c, &d, &p);

        assert_eq!(r, BigUint::from(7u32));
    }

    #[test]
    fn test_mul_below_p() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(53u32);

        let r = FiniteField::mul(&c, &d, &p);

        assert_eq!(r, BigUint::from(40u32));
    }

    #[test]
    fn test_inverse_add() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(53u32);

        let c_inv = FiniteField::inverse_add(&c, &p);

        assert_eq!(c_inv, BigUint::from(49u32));
        assert_eq!(FiniteField::add(&c, &c_inv, &p), BigUint::from(0u32));
    }

    #[test]
    #[should_panic(expected = "bigger or equal")]
    fn test_fail_inverse_add() {
        let c = BigUint::from(54u32);
        let p = BigUint::from(53u32);

        FiniteField::inverse_add(&c, &p);
    }

    #[test]
    fn test_inverse_mul() {
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let c_inv = FiniteField::inverse_mul_prime(&c, &p);

        assert_eq!(c_inv, BigUint::from(3u32));
        assert_eq!(FiniteField::mul(&c, &c_inv, &p), BigUint::from(1u32));
    }

    #[test]
    fn test_sub_above_p() {
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::sub(&c, &d, &p);

        assert_eq!(r, BigUint::from(5u32));
    }

    #[test]
    fn test_sub_below_p() {
        let c = BigUint::from(10u32);
        let d = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::sub(&c, &d, &p);

        assert_eq!(r, BigUint::from(6u32));
    }

    #[test]
    fn test_div_above_p() {
        let c = BigUint::from(10u32);
        let d = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::div(&c, &d, &p);

        assert_eq!(r, BigUint::from(8u32));
    }

    #[test]
    fn test_div_below_p() {
        let c = BigUint::from(2u32);
        let d = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::div(&c, &d, &p);

        assert_eq!(r, BigUint::from(6u32));
    }
}
