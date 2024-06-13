use lambdaworks_math::field::fields::u64_prime_field::U64FieldElement;

const MODULUS: u64 = 0xc0000001;
type FE = U64FieldElement<MODULUS>;    

fn main() {
    // create vec to hold fibonacci square sequence
    let mut fib_sq = Vec::<FE>::with_capacity(1024);
    fib_sq.push(FE::new(1));
    fib_sq.push(FE::new(3141592));
}
