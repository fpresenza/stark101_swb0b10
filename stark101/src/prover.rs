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
use crate::common::{self, PublicInput, VectorCommitment, StarkProof};
use crate::fri;

// the stark252 field has 2-adicity of 192, i.e., the largest
// multiplicative subgroup whose order is a power of two has order 2^192
type F = Stark252PrimeField;
type FE = FieldElement<F>;

pub fn generate_proof(public_input: PublicInput<F>) -> StarkProof<F> {

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
    let witness = FE::from(3141592_u64);
    let interp_order: usize = 1 << interp_two_power;
    let eval_order: usize = 1 << eval_two_power;


    // define primitive root
    let g = F::get_primitive_root_of_unity(interp_two_power as u64).unwrap();
    let g_to_the_1021 = g.pow(1021_u64);
    let g_to_the_1022 = g * g_to_the_1021;
    let g_to_the_1023 = g * g_to_the_1022;


    // create vec to hold fibonacci square sequence
    let mut fib_squared = Vec::<FE>::with_capacity(interp_order);
    fib_squared.push(fib_squared_0);
    fib_squared.push(witness);

    for i in 2..interp_order {
        let x = fib_squared[i-2];
        let y = fib_squared[i-1];
        fib_squared.push(x.square() + y.square());
    }

    // fft-interpolate the fibonacci square sequence
    let trace_poly = match Polynomial::interpolate_fft::<F>(&fib_squared) {
        Ok(p) => p,
        Err(e) => panic!("{:?}", e),
    };

    // fft-evaluate the fibonacci square sequence over a larger domain
    // of size (blow-up factor) * (interpolation domain size)
    // the offset is obtained as an outside not in the interpolation domain
    let offset = FE::from(2_u64);
    let trace_poly_eval = Polynomial::evaluate_offset_fft::<F>(
        &trace_poly, 1, Some(eval_order), &offset
    ).unwrap();

    // commit to the trace evaluations over the larger domain using a merkle tree
    let trace_poly_tree = MerkleTree::<Keccak256Backend<F>>::build(&trace_poly_eval);
    transcript.append_bytes(&trace_poly_tree.root);
    let mut trace_commitment = VectorCommitment::<F> {
        root: trace_poly_tree.root,
        inclusion_proofs: vec![]
    };

    // ===================================
    // =========|    Part 2:   |==========
    // ===== Polynomial Constraints ======
    // ===================================
    let x = Polynomial::new_monomial(one, 1);
    let x_to_the_1024 = Polynomial::new_monomial(one, interp_order);

    // initial element constraint
    let constraint_0_poly = poly::polynomial_division(
        &(&trace_poly - fib_squared_0),
        &(&x - one),
        eval_order,
        &offset
    );

    // result element constraint
    let constraint_1022_poly = poly::polynomial_division(
        &(&trace_poly - fib_squared_1022),
        &(&x - g_to_the_1022),
        eval_order,
        &offset
    );

    // trace transition constraint
    // numerator
    let trace_poly_scaled_once = trace_poly.scale(&g);
    let trace_poly_scaled_twice = trace_poly_scaled_once.scale(&g);
    let trace_poly_squared = poly::polynomial_power(
        &trace_poly,
        2_u64,
        eval_order,
        &offset
    );
    let trace_poly_scaled_once_squared = poly::polynomial_power(
        &trace_poly_scaled_once,
        2_u64,
        eval_order,
        &offset
    );

    let numerator = poly::polynomial_multiplication(
        &[
            &(trace_poly_scaled_twice - trace_poly_scaled_once_squared - trace_poly_squared),
            &(&x - g_to_the_1021), 
            &(&x - g_to_the_1022),
            &(&x - g_to_the_1023)
        ],
        eval_order,
        &offset
    );
    // denominator
    let denominator = &x_to_the_1024 - one;
    // polynomial
    let transition_constraint_poly = poly::polynomial_division(
        &numerator,
        &denominator,
        eval_order,
        &offset
    );

    // composition polynomial
    let a = transcript.sample_field_element();
    let b = transcript.sample_field_element();
    let c = transcript.sample_field_element();
    let comp_poly = a * constraint_0_poly + b * constraint_1022_poly + c * transition_constraint_poly;

    // ===================================
    // =========|    Part 3:   |==========
    // ========= FRI Commitment ==========
    // ===================================
    // get queries evaluations and add to transcript
    let query_indices = common::sample_queries(num_queries, eval_order, &mut transcript);
    let aux_indices = [0_usize, 8, 16];
    let all_indices = query_indices
        .iter()
        .map(|i| {
            aux_indices
                .iter()
                .map(|j| (i + j) % eval_order)
                .collect::<Vec<usize>>()
    }).collect::<Vec<Vec<usize>>>()
    .concat();

    trace_commitment.inclusion_proofs.extend(common::generate_inclusion_proofs(
        &all_indices,
        &trace_poly_eval,
        &trace_poly_tree,
    ));
    // let trace_commitment = VectorCommitment(trace_poly_tree.root, trace_poly_incl_proofs);
    
    // build fri layers
    let composition_commitment = fri::commit_and_fold(
        &comp_poly,
        eval_order,
        &offset,
        query_indices,
        &mut transcript
    );


    StarkProof {
        trace_commitment,
        composition_commitment
    }

}
