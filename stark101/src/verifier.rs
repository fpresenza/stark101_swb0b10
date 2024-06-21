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
use crate::common::{self, PublicInput, InclusionProof, StarkProof};
use crate::fri::{self, FriLayer};

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FE = FieldElement<F>;

pub fn verify_proof(public_input: PublicInput<F>, proof: StarkProof<F>) -> bool {
    println!(
        "
        ===================================
        =========|   VERIFIER   |==========
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

    let StarkProof(
        trace_poly_tree_root,
        trace_poly_incl_proofs,
        fri_layers
    ) = proof;

    // initialize result
    let mut result = false;

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
    let offset = FE::from(2_u64); 

    /*
        TODO: OFFSET IS PUBLIC INPUT
    */

    // define primitive root
    let power_of_two = usize::BITS - int_dom_size.leading_zeros() - 1;
    let g = F::get_primitive_root_of_unity(power_of_two as u64).unwrap();
    let g_to_the_1021 = g.pow(1021_u64);
    let g_to_the_1022 = g * g_to_the_1021;
    let g_to_the_1023 = g * g_to_the_1022;


    transcript.append_bytes(&trace_poly_tree_root);
    println!("Appending root of trace polynomial to transcript.");

    // ===================================
    // =========|    Part 2:   |==========
    // ===== Polynomial Constraints ======
    // ===================================
    let a = transcript.sample_field_element();
    let b = transcript.sample_field_element();
    let c = transcript.sample_field_element();
    // println!("a = {:?}", a);
    // println!("b = {:?}", b);
    // println!("c = {:?}", c);


    // ===================================
    // =========|    Part 3:   |==========
    // ========= FRI Commitment ==========
    // ===================================
    // get queries evaluations and add to transcript
    let query_indices = common::sample_queries(num_queries, eval_dom_size, &mut transcript);
    println!("Sampling Query indices and appending to transcript: {:?}", query_indices);

    for i in 0..query_indices.len() {
        let indices = [
            query_indices[i],
            (query_indices[i] + 1) % eval_dom_size,
            (query_indices[i] + 2) % eval_dom_size
        ];
        let incl_proofs = &trace_poly_incl_proofs[i];

        let incl_proofs_result = incl_proofs
            .iter()
            .zip(indices)
            .map(|(InclusionProof(eval, proof), j)| {
                proof.verify::<Keccak256Backend<F>>(
                    &trace_poly_tree_root,
                    j,
                    &eval
                )
            })
            .fold(true, |agg, res| {agg && res});
        if !incl_proofs_result {
            println!("Verification of composition polynomial inclusion proofs did not pass");
            return false
        }
    }

    true
}
