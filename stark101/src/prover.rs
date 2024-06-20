use lambdaworks_math::traits::ByteConversion;
use lambdaworks_math::field::{
    traits::IsFFTField,
    fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    element::FieldElement
};
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_crypto::merkle_tree::{
    merkle::MerkleTree,
    backends::types::Keccak256Backend
};
use lambdaworks_crypto::fiat_shamir::{
    is_transcript::IsTranscript,
    default_transcript::DefaultTranscript
};

use crate::poly;
use crate::common::PublicInput;
use crate::fri;

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FE = FieldElement<F>;

pub fn generate_proof(public_input: PublicInput<F>) {
    println!(
        "
        ===================================
        ==========|   PROVER   |===========
        ===================================
        "
    );
    // ===================================
    // ==========|    Part 1:   |=========
    // === Statement, LDE & Commitment ===
    // ===================================
    // extract public input
    let PublicInput(
        modulus,
        int_dom_size,
        eval_dom_size,
        num_queries,
        fib_squared_0,
        fib_squared_1022
    ) = public_input;

    // initialize transcript and append all public inputs
    let mut transcript = DefaultTranscript::<F>::new(&[]);
    transcript.append_bytes(&modulus.to_bytes_be());
    transcript.append_bytes(&int_dom_size.to_be_bytes());
    transcript.append_bytes(&eval_dom_size.to_be_bytes());
    transcript.append_bytes(&num_queries.to_be_bytes());
    transcript.append_bytes(&fib_squared_0.to_bytes_be());
    transcript.append_bytes(&fib_squared_1022.to_bytes_be());
    println!("Appending public inputs to transcript.");

    // define example parameters
    let one = FE::one();
    let witness = FE::from(3141592_u64);

    // define primitive root
    let power_of_two = usize::BITS - int_dom_size.leading_zeros() - 1;
    let g = F::get_primitive_root_of_unity(power_of_two as u64).unwrap();
    let g_to_the_1021 = g.pow(1021_u64);
    let g_to_the_1022 = g * g_to_the_1021;
    let g_to_the_1023 = g * g_to_the_1022;


    // create vec to hold fibonacci square sequence
    let mut fib_squared = Vec::<FE>::with_capacity(int_dom_size);
    fib_squared.push(fib_squared_0);
    fib_squared.push(witness);

    for i in 2..int_dom_size {
        let x = fib_squared[i-2];
        let y = fib_squared[i-1];
        fib_squared.push(x.square() + y.square());
    }

    // fft-interpolate the fibonacci square sequence
    let trace_poly = match Polynomial::interpolate_fft::<F>(&fib_squared) {
        Ok(p) => p,
        Err(e) => panic!("{:?}", e),
    };
    assert_eq!(trace_poly.degree(), int_dom_size - 1);
    assert_eq!(trace_poly.evaluate(&one), fib_squared_0);
    assert_eq!(trace_poly.evaluate(&g), witness);
    assert_eq!(trace_poly.evaluate(&g_to_the_1022), fib_squared_1022);

    // fft-evaluate the fibonacci square sequence over a larger domain
    // of size (blow-up factor) * (interpolation domain size)
    // the offset is obtained as an outside not in the interpolation domain
    let offset = FE::from(2_u64);
    assert!(offset.pow(int_dom_size as u64) != one);
    let trace_eval = Polynomial::evaluate_offset_fft::<F>(
        &trace_poly, 1, Some(eval_dom_size), &offset
    ).unwrap();

    // commit to the trace evaluations over the larger domain using a merkle tree
    let trace_poly_merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&trace_eval);
    transcript.append_bytes(&trace_poly_merkle_tree.root);
    println!("Appending root of trace polynomial to transcript.");

    // ===================================
    // =========|    Part 2:   |==========
    // ===== Polynomial Constraints ======
    // ===================================
    let x = Polynomial::new_monomial(one, 1);
    let x_to_the_1024 = Polynomial::new_monomial(one, int_dom_size);

    // initial element constraint
    let contraint_0_poly = poly::polynomial_division(
        &(&trace_poly - fib_squared_0),
        &(&x - one),
        eval_dom_size,
        &offset
    );
    assert_eq!(contraint_0_poly.degree(), int_dom_size - 2);

    // result element constraint
    let contraint_1022_poly = poly::polynomial_division(
        &(&trace_poly - fib_squared_1022),
        &(&x - g_to_the_1022),
        eval_dom_size,
        &offset
    );
    assert_eq!(contraint_1022_poly.degree(), int_dom_size - 2);

    // trace transition constraint
    // numerator
    let trace_poly_scaled_once = trace_poly.scale(&g);
    let trace_poly_scaled_twice = trace_poly_scaled_once.scale(&g);
    let trace_poly_squared = poly::polynomial_power(
        &trace_poly,
        2_u64,
        eval_dom_size,
        &offset
    );
    let trace_poly_scaled_once_squared = poly::polynomial_power(
        &trace_poly_scaled_once,
        2_u64,
        eval_dom_size,
        &offset
    );
    assert_eq!(trace_poly_squared.degree(), 2 * int_dom_size - 2);
    assert_eq!(trace_poly_scaled_once_squared.degree(), 2 * int_dom_size - 2);
    assert_eq!(
        trace_poly_scaled_twice.evaluate(&g_to_the_1021),
        trace_poly_scaled_once_squared.evaluate(&g_to_the_1021) + trace_poly_squared.evaluate(&g_to_the_1021)
    );

    let numerator = poly::polynomial_multiplication(
        &[
            &(trace_poly_scaled_twice - trace_poly_scaled_once_squared - trace_poly_squared),
            &(&x - g_to_the_1021), 
            &(&x - g_to_the_1022),
            &(&x - g_to_the_1023)
        ],
        eval_dom_size,
        &offset
    );
    // denominator
    let denominator = &x_to_the_1024 - one;
    // polynomial
    let transition_constraint_poly = poly::polynomial_division(
        &numerator,
        &denominator,
        eval_dom_size,
        &offset
    );
    assert_eq!(transition_constraint_poly.degree(), int_dom_size + 1);

    // composition polynomial
    let a = transcript.sample_field_element();
    let b = transcript.sample_field_element();
    let c = transcript.sample_field_element();
    let comp_poly = a * contraint_0_poly + b * contraint_1022_poly + c * transition_constraint_poly;
    assert_eq!(comp_poly.degree(), int_dom_size + 1);

    // ===================================
    // =========|    Part 4:   |==========
    // ========= FRI Commitment ==========
    // ===================================
    let fri_layers = fri::commit_and_fold(
        &comp_poly,
        eval_dom_size,
        &offset,
        num_queries,
        &mut transcript
    );

}
