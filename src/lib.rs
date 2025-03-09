//! Conversion between geodetic and geocentric (cartesian) coordinates.

use core::f64::consts::FRAC_PI_6;

/// Convert from geodetic coordinates to geocentric (cartesian) coordinates.
///
/// # Arguments
/// * `a` - Semi-major axis
/// * `e_sq` - Eccentricity squared
/// * `lng` - Longitude in degrees
/// * `lat` - Latitude in degrees
/// * `height` - Height in meters above the ellipsoid
#[inline]
pub fn geodetic_to_geocentric(
    a: f64,
    e_sq: f64,
    lng: f64,
    lat: f64,
    height: f64,
) -> (f64, f64, f64) {
    let (lam, phi) = (lng.to_radians(), lat.to_radians());
    let n = if e_sq == 0.0 {
        a // optimization for sphere (e=0)
    } else {
        a / (1. - e_sq * phi.sin() * phi.sin()).sqrt()
    };
    let x = (n + height) * phi.cos() * lam.cos();
    let y = (n + height) * phi.cos() * lam.sin();
    let z = (n * (1. - e_sq) + height) * phi.sin(); // updated to use e_sq instead of ellips.e_sq()
    (x, y, z)
}

/// Convert from geocentric (cartesian) coordinates to geodetic coordinates.
///
/// This function implements Hugues Vermeille's *[An analytical method to transform geocentric into geodetic coordinates](https://DOI.org/10.1007/s00190-010-0419-x)*, J. Geodesy (2011) 85, page 105-117.
///
/// # Arguments
/// * `a` - Semi-major axis
/// * `e_sq` - Eccentricity squared
/// * `x` - X coordinate
/// * `y` - Y coordinate
/// * `z` - Z coordinate
#[inline]
pub fn geocentric_to_geodetic(a: f64, e_sq: f64, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let e_quad = e_sq * e_sq;

    let p = (x * x + y * y) / (a * a);
    let q = (1. - e_sq) * z * z / (a * a);
    let r = (p + q - e_quad) / 6.;

    let evol = 8. * r.powi(3) + e_quad * p * q;

    let (phi, h) = if evol > 0. || q != 0. {
        let u = if evol > 0. {
            // Outside the evolute
            let l = (evol.sqrt() + (e_quad * p * q).sqrt()).cbrt();
            (3. * r * r) / (2. * l * l) + 0.5 * (l + r / l) * (l + r / l)
        } else {
            // On or inside the evolute and not in the singular disc
            let t = 2. / 3.
                * ((e_quad * p * q).sqrt()).atan2((-evol).sqrt() + (-8. * r * r * r).sqrt());
            -4. * r * (t).sin() * (FRAC_PI_6 + t).cos()
        };
        let v = (u * u + e_quad * q).sqrt();
        let w = e_sq * (u + v - q) / (2. * v);
        let k = (u + v) / ((w * w + u + v).sqrt() + w);
        let d = k * (x * x + y * y).sqrt() / (k + e_sq);
        let h = (k + e_sq - 1.) / k * (d * d + z * z).sqrt();
        let phi = 2. * z.atan2(d + (d * d + z * z).sqrt());
        (phi, h)
    } else {
        // In the singular disc
        let h = -a * ((1. - e_sq) * (e_sq - p) / e_sq).sqrt();
        let phi =
            2. * ((e_quad - p).sqrt()).atan2((e_sq * (e_sq - p)).sqrt() + ((1. - e_sq) * p).sqrt());
        (phi, h)
    };

    let lam = f64::atan2(y, x);

    (lam.to_degrees(), phi.to_degrees(), h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip() {
        let a = 6378137.;
        let inv_f = 298.257223563;
        let f = 1. / inv_f;
        let b = a * (1. - f);
        let e_sq = f * (2. - f);

        // Outside the evolute
        {
            let (lng, lat, height) = (140., 37., 50.);
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, lng, lat, height);
            let (lng2, lat2, height2) = geocentric_to_geodetic(a, e_sq, x, y, z);
            assert!((lng - lng2).abs() < 1e-10);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((height - height2).abs() < 1e-7);
        }

        // On or inside the evolute and not in the singular disc
        {
            let (lng, lat, height) = (45., 74.58501644931525, -6344866.234164982);
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, lng, lat, height);
            let (lng2, lat2, height2) = geocentric_to_geodetic(a, e_sq, x, y, z);
            assert!((lng - lng2).abs() < 1e-10);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((height - height2).abs() < 1e-7);
        }

        // In the singular disc
        {
            let (lng, lat, height) = (120., 88.10828645, -6356728.972246517);
            let (x, y, _) = geodetic_to_geocentric(a, e_sq, lng, lat, height);
            let z = 0.;
            let (lng2, lat2, height2) = geocentric_to_geodetic(a, e_sq, x, y, z);
            assert!((lng - lng2).abs() < 1e-10);
            assert!((lat - lat2).abs() < 1e-10);
            assert!((height - height2).abs() < 30.);
        }

        // Earth's center
        {
            let (lng, lat, height) = (0., 90., b);
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, lng, lat, height);
            let (lng2, lat2, height2) = geocentric_to_geodetic(a, e_sq, x, y, z);
            assert!((lng - lng2).abs() < 1e-9);
            assert!((lat - lat2).abs() < 1e-9);
            assert!((height - height2).abs() < 1e-7);
        }
    }

    #[test]
    fn to_geocentric() {
        let a = 6378137.;
        let inv_f = 298.257223563;
        let f = 1. / inv_f;
        let b = a * (1. - f);
        let e_sq = f * (2. - f);

        {
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, 140., 37., 50.);
            assert!((x - -3906851.9770472576).abs() < 1e-10);
            assert!((y - 3278238.0530045824).abs() < 1e-10);
            assert!((z - 3817423.251099322).abs() < 1e-10);
        }

        // north pole
        {
            let height = 150.;
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, 123., 90., height);
            assert!((x - 0.).abs() < 1e-9);
            assert!((y - 0.).abs() < 1e-9);
            assert!((z - (b + height)).abs() < 1e-9);
        }

        // null island
        {
            let height = 100.;
            let (x, y, z) = geodetic_to_geocentric(a, e_sq, 0., 0., height);
            assert!((x - (a + height)).abs() < 1e-9);
            assert!((y - 0.).abs() < 1e-9);
            assert!((z - 0.).abs() < 1e-9);
        }
    }
}
