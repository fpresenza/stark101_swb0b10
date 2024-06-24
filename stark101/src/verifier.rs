use lambdaworks_math::traits::ByteConversion;
use lambdaworks_math::field::{
    traits::IsFFTField,
    fields::fft_friendly::stark_252_prime_field::Stark252PrimeField,
    element::FieldElement
};
use lambdaworks_crypto::fiat_shamir::{
    is_transcript::IsTranscript,
    default_transcript::DefaultTranscript
};

use crate::common::{self, PublicInput, VectorCommitment, StarkProof};
use crate::fri;

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FE = FieldElement<F>;

pub fn verify_proof(public_input: PublicInput<F>, stark_proof: StarkProof<F>) -> bool {

    // ===================================
    // ==========|    Part 1:   |=========
    // === Statement, LDE & Commitment ===
    // ===================================
    // extract public input
    let PublicInput(
        modulus,
        interp_two_power,
        eval_two_power,
        num_queries,
        fib_squared_0,
        fib_squared_1022
    ) = public_input;

    let StarkProof(
        VectorCommitment(
            trace_poly_root,
            trace_poly_proofs
        ),
        fri_layers
    ) = stark_proof;

    // initialize transcript and append all public inputs
    let mut transcript = DefaultTranscript::<F>::new(&[]);
    transcript.append_bytes(&modulus.to_bytes_be());
    transcript.append_bytes(&interp_two_power.to_be_bytes());
    transcript.append_bytes(&eval_two_power.to_be_bytes());
    transcript.append_bytes(&num_queries.to_be_bytes());
    transcript.append_bytes(&fib_squared_0.to_bytes_be());
    transcript.append_bytes(&fib_squared_1022.to_bytes_be());

    // define example parameters
    let one = FE::one();
    let offset = FE::from(2_u64); 
    let eval_order: usize = 1 << eval_two_power;

    /*
        TODO: OFFSET IS PUBLIC INPUT
    */

    // define primitive root
    let g = F::get_primitive_root_of_unity(interp_two_power as u64).unwrap();
    let g_to_the_1021 = g.pow(1021_u64);
    let g_to_the_1022 = g * g_to_the_1021;
    let g_to_the_1023 = g * g_to_the_1022;

    let w = F::get_primitive_root_of_unity(eval_two_power as u64).unwrap();

    transcript.append_bytes(&trace_poly_root);

    // ===================================
    // =========|    Part 2:   |==========
    // ===== Polynomial Constraints ======
    // ===================================
    let a = transcript.sample_field_element();
    let b = transcript.sample_field_element();
    let c = transcript.sample_field_element();

    // get queries evaluations and add to transcript
    let query_indices = common::sample_queries(num_queries, eval_order, &mut transcript);
    let aux_indices = [0_usize, 8, 16];
    let aux_indices_len = aux_indices.len();
    let all_indices = query_indices
        .iter()
        .map(|i| {
            aux_indices
                .iter()
                .map(|j| (i + j) % eval_order)
                .collect::<Vec<usize>>()
    }).collect::<Vec<Vec<usize>>>()
    .concat();

    // verify trace inclusion proofs
    if !common::verify_inclusion_proofs(&all_indices, &trace_poly_proofs, trace_poly_root) {
        return false
    }

    // compute queries
    let queries = query_indices
        .iter()
        .map(|idx| offset * w.pow(idx.to_owned()))
        .collect::<Vec<FE>>();

    // compute composition polynomial evaluations
    let comp_poly_query_evals = queries
        .iter()
        .enumerate()
        .map(|(i, x0)| {
            let t = (0..aux_indices_len).map(|k| {
                trace_poly_proofs[aux_indices_len * i + k].0
            }).collect::<Vec<FE>>();
            a * (t[0] - fib_squared_0) / (x0 - one) +
            b * (t[0] - fib_squared_1022) / (x0 - g_to_the_1022) +
            c * (
                    (t[2] - t[1].square() - t[0].square()) * 
                    (x0 - g_to_the_1021) * 
                    (x0 - g_to_the_1022) * 
                    (x0 - g_to_the_1023) / 
                    (x0.pow(1024_u64) - one)
            )
        }).collect::<Vec<FE>>();

    // ===================================
    // =========|    Part 3:   |==========
    // ======== FRI Decommitment =========
    // ===================================
    // build fri layers
    fri::decommit_and_fold(
        &fri_layers,
        &eval_order,
        &query_indices,
        &queries,
        &comp_poly_query_evals,
        &mut transcript
    )
}
