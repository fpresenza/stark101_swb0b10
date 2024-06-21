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

pub fn verify_proof(public_input: PublicInput<F>, stark_proof: StarkProof<F>) -> bool {
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
        trace_poly_proofs,
        fri_layers
    ) = stark_proof;

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

    let power_of_two = usize::BITS - eval_dom_size.leading_zeros() - 1;
    let w = F::get_primitive_root_of_unity(power_of_two as u64).unwrap();

    transcript.append_bytes(&trace_poly_tree_root);
    println!("Appending root of trace polynomial to transcript.");

    // ===================================
    // =========|    Part 2:   |==========
    // ===== Polynomial Constraints ======
    // ===================================
    let a = transcript.sample_field_element();
    let b = transcript.sample_field_element();
    let c = transcript.sample_field_element();

    // get queries evaluations and add to transcript
    let query_indices = common::sample_queries(num_queries, eval_dom_size, &mut transcript);
    println!("Sampling Query indices and appending to transcript: {:?}", query_indices);

    // verify trace inclusion proofs
    for (index, t_proof) in query_indices.iter().zip(trace_poly_proofs) {
        let result = t_proof
            .iter()
            .enumerate()
            .map(|(k, InclusionProof(eval, proof))| {
                proof.verify::<Keccak256Backend<F>>(
                    &trace_poly_tree_root,
                    (index + 8*k) % eval_dom_size,
                    &eval
                )
            })
            .fold(true, |res, valid| {res && valid});
        if !result {
            println!("Verification of trace polynomial inclusion proofs did not pass");
            return false
        }

        let t = [t_proof[0].0, t_proof[1].0, t_proof[2].0];
        let x0 = offset * w.pow(index.to_owned());
        let cp_eval = 
            a * (t[0] - fib_squared_0) / (x0 - one) +
            b * (t[0] - fib_squared_1022) / (x0 - g_to_the_1022) +
            c * (
                    (t[2] - t[1].square() - t[0].square()) * 
                    (x0 - g_to_the_1021) * 
                    (x0 - g_to_the_1022) * 
                    (x0 - g_to_the_1023) / 
                    (x0.pow(1024_u64) - one)
            );
    }

    // verify composition polynomial inclusion proofs

    // ===================================
    // =========|    Part 3:   |==========
    // ========= FRI Decommitment ==========
    // ===================================

    true
}
