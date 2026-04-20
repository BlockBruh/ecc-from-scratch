use num_bigint::BigUint;

use crate::finite_fields::FiniteField;

/// A point on an elliptic curve — either a pair of affine coordinates or
/// the point at infinity (the group identity element).
#[derive(PartialEq, Clone, Debug)]
pub enum Point {
    Coordonate { x: BigUint, y: BigUint },
    Identity,
}

/// An elliptic curve over a prime field **F_p** in short Weierstrass form:
///
/// ```text
/// y² ≡ x³ + a·x + b  (mod p)
/// ```
#[derive(Debug)]
pub struct EllipticCurve {
    pub a: BigUint,
    pub b: BigUint,
    pub p: BigUint,
}

impl EllipticCurve {
    /// Adds two **distinct** points on the curve using the chord-and-tangent rule.
    ///
    /// Use [`double`] when the two points are equal.
    ///
    /// # Panics
    /// Panics if either point is not on the curve, or if `c == d`
    /// (use [`double`] for that case).
    ///
    /// [`double`]: EllipticCurve::double
    pub fn add(&self, c: &Point, d: &Point) -> Point {
        assert!(self.is_on_curve(c), "First point is not in curve");
        assert!(self.is_on_curve(d), "Second point is not in curve");
        assert!(c != d, "Points should be different from each other");

        match (c, d) {
            (Point::Identity, _) => d.clone(),
            (_, Point::Identity) => c.clone(),
            (Point::Coordonate { x: x1, y: y1 }, Point::Coordonate { x: x2, y: y2 }) => {
                let y1plusy2 = FiniteField::add(y1, y2, &self.p);
                if x1 == x2 && y1plusy2 == BigUint::from(0u32) {
                    return Point::Identity;
                }

                let numerator = FiniteField::sub(y2, y1, &self.p);
                let denominator = FiniteField::sub(x2, x1, &self.p);
                let s = FiniteField::div(&numerator, &denominator, &self.p);

                self.point_from_slope(&s, x1, x2, y1)
            }
        }
    }

    /// Doubles a point on the curve, i.e. computes `P + P` using the tangent line.
    ///
    /// Returns [`Point::Identity`] when `P` is the identity or when `y = 0`
    /// (the tangent is vertical).
    ///
    /// # Panics
    /// Panics if the point is not on the curve.
    pub fn double(&self, c: &Point) -> Point {
        assert!(self.is_on_curve(c), "Point is not in curve");
        match c {
            Point::Identity => Point::Identity,
            Point::Coordonate { x, y } => {
                if *y == BigUint::from(0u32) {
                    return Point::Identity;
                }

                let xpow2 = x.modpow(&BigUint::from(2u32), &self.p);
                let numerator = FiniteField::add(
                    &FiniteField::mul(&BigUint::from(3u32), &xpow2, &self.p),
                    &self.a,
                    &self.p,
                );
                let denominator = FiniteField::mul(&BigUint::from(2u32), y, &self.p);
                let s = FiniteField::div(&numerator, &denominator, &self.p);

                self.point_from_slope(&s, x, x, y)
            }
        }
    }

    /// Computes a new point from a slope `s`, two x-coordinates, and `y1`:
    ///
    /// ```text
    /// x3 = s² - x1 - x2  (mod p)
    /// y3 = s·(x1 - x3) - y1  (mod p)
    /// ```
    fn point_from_slope(&self, s: &BigUint, x1: &BigUint, x2: &BigUint, y1: &BigUint) -> Point {
        let s2 = s.modpow(&BigUint::from(2u32), &self.p);
        let x3 = FiniteField::sub(&FiniteField::sub(&s2, x1, &self.p), x2, &self.p);
        let y3 = FiniteField::sub(
            &FiniteField::mul(s, &FiniteField::sub(x1, &x3, &self.p), &self.p),
            y1,
            &self.p,
        );
        Point::Coordonate { x: x3, y: y3 }
    }

    /// Computes the scalar multiple `d·P` using the double-and-add algorithm.
    ///
    /// The algorithm processes the bits of `d` from most-significant to
    /// least-significant.  Each step doubles the accumulator; when a bit is 1
    /// it also adds the base point `c`.
    pub fn scalar_mul(&self, c: &Point, d: &BigUint) -> Point {
        let mut t = c.clone();
        for i in (0..(d.bits() - 1)).rev() {
            t = self.double(&t);
            if d.bit(i) {
                t = self.add(&t, c);
            }
        }
        t
    }

    /// Returns `true` if the point satisfies the curve equation `y² ≡ x³ + ax + b (mod p)`.
    /// The identity point always returns `true`.
    pub fn is_on_curve(&self, c: &Point) -> bool {
        match c {
            Point::Coordonate { x, y } => {
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let ax = FiniteField::mul(&self.a, x, &self.p);

                y2 == FiniteField::add(&FiniteField::add(&x3, &ax, &self.p), &self.b, &self.p)
            }
            Point::Identity => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ec_point_addition() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        let p1 = Point::Coordonate {
            x: BigUint::from(6u32),
            y: BigUint::from(3u32),
        };
        let p2 = Point::Coordonate {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };
        let pr = Point::Coordonate {
            x: BigUint::from(10u32),
            y: BigUint::from(6u32),
        };

        let r = ec.add(&p1, &p2);
        assert_eq!(r, pr);

        let r = ec.add(&p2, &p1);
        assert_eq!(r, pr);
    }

    #[test]
    fn test_ec_point_double() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        let p1 = Point::Coordonate {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };
        let expected = Point::Coordonate {
            x: BigUint::from(6u32),
            y: BigUint::from(3u32),
        };

        assert_eq!(ec.double(&p1), expected);
    }

    #[test]
    fn test_ec_point_double_identity() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        assert_eq!(ec.double(&Point::Identity), Point::Identity);
    }

    #[test]
    fn test_ec_point_addition_opposite() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        let p1 = Point::Coordonate {
            x: BigUint::from(5u32),
            y: BigUint::from(16u32),
        };
        let p2 = Point::Coordonate {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };

        let r = ec.add(&p1, &p2);
        assert_eq!(r, Point::Identity);
    }

    #[test]
    fn test_ec_point_scalar_mul() {
        let ec = EllipticCurve {
            a: BigUint::from(2u32),
            b: BigUint::from(2u32),
            p: BigUint::from(17u32),
        };

        let c = Point::Coordonate {
            x: BigUint::from(5u32),
            y: BigUint::from(1u32),
        };
        let expected = Point::Coordonate {
            x: BigUint::from(6u32),
            y: BigUint::from(3u32),
        };
        let r = ec.scalar_mul(&c, &BigUint::from(2u32));
        assert_eq!(r, expected);

        let expected = Point::Coordonate {
            x: BigUint::from(7u32),
            y: BigUint::from(11u32),
        };
        let r = ec.scalar_mul(&c, &BigUint::from(10u32));
        assert_eq!(r, expected);

        let r = ec.scalar_mul(&c, &BigUint::from(19u32));
        assert_eq!(r, Point::Identity);
    }

    #[test]
    fn test_ec_secp256r1() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
            16,
        )
        .unwrap();

        let n = BigUint::parse_bytes(
            b"FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551",
            16,
        )
        .unwrap();

        let gx = BigUint::parse_bytes(
            b"6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296",
            16,
        )
        .unwrap();

        let gy = BigUint::parse_bytes(
            b"4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5",
            16,
        )
        .unwrap();

        let a = BigUint::parse_bytes(
            b"FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC",
            16,
        )
        .unwrap();

        let b = BigUint::parse_bytes(
            b"5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B",
            16,
        )
        .unwrap();

        let ec = EllipticCurve { a, b, p };
        let g = Point::Coordonate { x: gx, y: gy };

        assert!(ec.is_on_curve(&g));
        assert_eq!(ec.scalar_mul(&g, &n), Point::Identity);
    }

    #[test]
    fn test_ec_secp256k1() {
        let p = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap();

        let n = BigUint::parse_bytes(
            b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap();

        let gx = BigUint::parse_bytes(
            b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .unwrap();

        let gy = BigUint::parse_bytes(
            b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .unwrap();

        let ec = EllipticCurve {
            a: BigUint::from(0u32),
            b: BigUint::from(7u32),
            p,
        };

        let g = Point::Coordonate { x: gx, y: gy };

        assert_eq!(ec.scalar_mul(&g, &n), Point::Identity);
    }
}
