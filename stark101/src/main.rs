use lambdaworks_math::field::{
    fields::montgomery_backed_prime_fields::IsModulus,
    fields::fft_friendly::stark_252_prime_field::{
        Stark252PrimeField,
        MontgomeryConfigStark252PrimeField
    },
    element::FieldElement
};

mod poly;
mod common;
mod fri;
mod prover;
mod verifier;

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FConfig = MontgomeryConfigStark252PrimeField;
type FE = FieldElement<F>;

// interpolation domain of size 1024 = 2^10
const INT_DOM_SIZE: usize = 0b10000000000;
// evaluation domain of size 8192 = 2^13 (blow-up factor is 2^3)
const EVAL_DOM_SIZE: usize = 0b10000000000000;
// number of queries in FRI
const NUM_QUERIES: usize = 10;


fn main() {

    // public input //
    // field properties
    let modulus = FConfig::MODULUS;

    // trace properties
    let fib_squared_0 = FE::one();
    let fib_squared_1022 = FE::from_hex_unchecked("6A317721EF632FF24FB815C9BBD4D4582BC7E21A43CFBDD89A8B8F0BDA68252");

    let public_input = common::PublicInput(
        modulus,
        INT_DOM_SIZE,
        EVAL_DOM_SIZE,
        NUM_QUERIES,
        fib_squared_0,
        fib_squared_1022,
    );

    // generate valid proof
    let proof = prover::generate_proof(public_input.clone());

    // simulate invalid proof
    let mut invalid_proof = proof.clone();
    invalid_proof.0[0] += 1;

    if verifier::verify_proof(public_input.clone(), proof) {
        println!("Valid Proof: successfully verified.");
    } else {
        println!("Valid Proof: could not be verified.");
    }

    if verifier::verify_proof(public_input.clone(), invalid_proof) {
        println!("Invalid Proof: successfully verified.");
    } else {
        println!("Invalid Proof: could not be verified.");
    }
}