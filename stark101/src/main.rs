use lambdaworks_math::field::fields::fft_friendly::stark_252_prime_field::Stark252PrimeField;
use lambdaworks_math::field::element::FieldElement;
use lambdaworks_math::polynomial::Polynomial;

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type FE = FieldElement<Stark252PrimeField>;

// interpolation domain of size 1024 = 2^10
const INT_DOM_SIZE: usize = 0b10000000000;

// evaluation domain of size 8192 = 2^13 (blow-up factor is 2^3)
const EVAL_DOM_SIZE: usize = 0b10000000000000;

fn main() {
    // create vec to hold fibonacci square sequence
    let mut fib_sq = Vec::<FE>::with_capacity(INT_DOM_SIZE);
    fib_sq.push(FE::from(1 as u64));
    fib_sq.push(FE::from(3141592 as u64));

    for i in 2..INT_DOM_SIZE {
        let x = fib_sq[i-2];
        let y = fib_sq[i-1];
        fib_sq.push(x.pow(2_u64) + y.pow(2_u64));
    } 

    // fft-interpolate the fibonacci square sequence
    let trace_poly = match Polynomial::interpolate_fft::<Stark252PrimeField>(&fib_sq) {
        Ok(poly) => poly,
        Err(e) => panic!("{:?}", e),
    };
    assert_eq!(trace_poly.coefficients.len(), INT_DOM_SIZE);
}
