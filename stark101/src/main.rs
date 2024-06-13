use lambdaworks_math::field::fields::u64_prime_field::U64FieldElement;

// modulus is 3221225473 = 3 * 2^30 + 1
const MODULUS: u64 = 0b11000000000000000000000000000001;
type FE = U64FieldElement<MODULUS>;

// interpolation domain of size 1024 = 2^10
const INT_DOM_SIZE: usize = 0b10000000000;

// evaluation domain of size 8192 = 2^13 (blow-up factor is 2^3)
const EVAL_DOM_SIZE: usize = 0b10000000000000;

fn main() {
    // create vec to hold fibonacci square sequence
    let mut fib_sq = Vec::<FE>::with_capacity(INT_DOM_SIZE);
    fib_sq.push(FE::new(1));
    fib_sq.push(FE::new(3141592));

    for i in 2..INT_DOM_SIZE {
        let x = fib_sq[i-2];
        let y = fib_sq[i-1];
        fib_sq.push(x.pow(2_u64) + y.pow(2_u64));
    }

    // define generator of evaluation domain (order 2^10)
    let generator_eval_dom = FE::new(1855261384);
}
