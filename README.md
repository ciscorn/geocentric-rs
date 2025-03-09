# geocentric-rs

[![Test](https://github.com/ciscorn/geocentric-rs/actions/workflows/Test.yml/badge.svg)](https://github.com/ciscorn/geocentric-rs/actions/workflows/Test.yml)
[![codecov](https://codecov.io/gh/ciscorn/geocentric-rs/graph/badge.svg?token=zWu8h2egPG)](https://codecov.io/gh/ciscorn/geocentric-rs)

Conversion between geodetic (geographic) coordinates and geocentric (cartesian) coordinates in Rust.

Keywords: EPSG:4326, EPSG:4979, EPSG:4978

## Usage

```rust
// WGS 84 Ellipsoid
let a = 6378137.; // Semi-major axis
let inv_f = 298.257223563; // Inverse flattening
let f = 1. / inv_f; // Flattening
let e_sq = f * (2. - f); // Eccentricity squared

// Convert from geodetic to geocentric
let (x, y, z) = geodetic_to_geocentric(a, e_sq, 140., 37., 50.);

// Convert from geocentric to geodetic
let (lng, lat, height) = geocentric_to_geodetic(a, e_sq, x, y, z);
```

## References

The `geocentric_to_geodetic` function implements the algorithm described in:

- Hugues Vermeille, *"An analytical method to transform geocentric into geodetic coordinates"*, Journal of Geodesy (2011) 85, pages 105-117. [DOI:10.1007/s00190-010-0419-x](https://doi.org/10.1007/s00190-010-0419-x)

## License

MIT

## Author

Taku Fukada ([@ciscorn](https://github.com/ciscorn))
