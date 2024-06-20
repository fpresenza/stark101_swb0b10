use lambdaworks_math::traits::{AsBytes, ByteConversion};
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
};
use lambdaworks_math::polynomial::Polynomial;
use lambdaworks_crypto::merkle_tree::{
    merkle::MerkleTree,
    backends::types::Keccak256Backend, 
    proof::Proof
};
use lambdaworks_crypto::fiat_shamir::{
    is_transcript::IsTranscript,
    default_transcript::DefaultTranscript
};
use lambdaworks_math::unsigned_integer::element::U256;

use crate::poly;

pub struct FriLayer<F: IsField> {
    merkle_root: [u8; 32],
    queries_eval: Vec<FieldElement<F>>,
    // queries_merkle_proof: Vec<Proof<[u8; 32]>>,
    queries_sym_eval: Vec<FieldElement<F>>,
    // queries_sym_merkle_proof: Vec<Proof<[u8; 32]>>,
}

impl<F: IsField> FriLayer<F> {
    fn new(
        merkle_root: [u8; 32],
        queries_eval: Vec<FieldElement<F>>,
        // queries_merkle_proof: Vec<Proof<[u8; 32]>>,
        queries_sym_eval: Vec<FieldElement<F>>,
        // queries_sym_merkle_proof: Vec<Proof<[u8; 32]>>
    ) -> Self {
        Self {
            merkle_root,
            queries_eval,
            // queries_merkle_proof,
            queries_sym_eval,
            // queries_sym_merkle_proof
        }
    }
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        num_queries: usize,
        transcript: &mut DefaultTranscript<F>
    ) -> (Vec<FriLayer<F>>, FieldElement<F>)
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let mut polynomial = polynomial.clone();
    let mut offset = offset.clone();
    let number_of_layers = (usize::BITS - polynomial.degree().leading_zeros() + 1) as usize;
    let mut fri_layers = Vec::<FriLayer<F>>::with_capacity(number_of_layers);

    // commit to evaluations
    let (eval, root) = commit(&polynomial, domain_size, &offset, transcript);

     // sample queries
    let query_indices = sample_queries(num_queries, domain_size, transcript);

    /*
        TODO: 
        get queries evaluations and generate proofs
    */

    // append layer
    fri_layers.push(
        FriLayer::<F>::new(
            root,
            vec![FieldElement::<F>::from(0_u64)],
            // vec![Proof { merkle_path [0u8; 32]}],
            vec![FieldElement::<F>::from(0_u64)],
            // vec![Proof { merkle_path [0u8; 32]}],
        )
    );

    // recursive foldings
    for _ in 1..number_of_layers {
        let beta = transcript.sample_field_element();
    
        polynomial = poly::fold_polynomial(&polynomial, &beta);
        domain_size /= 2;
        offset = offset.square();

        let (eval, root) = commit(&polynomial, domain_size, &offset, transcript);
        fri_layers.push(
            FriLayer::<F>::new(
                root,
                vec![FieldElement::<F>::from(0_u64)],
                // vec![Proof { merkle_path [0u8; 32]}],
                vec![FieldElement::<F>::from(0_u64)],
                // vec![Proof { merkle_path [0u8; 32]}],
            )
        );
    }

    let constant_poly = polynomial.coefficients.first().unwrap();
    transcript.append_bytes(&constant_poly.to_bytes_be());

    (fri_layers, constant_poly.clone())
}


fn commit<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        domain_size: usize,
        offset: &FieldElement<F>,
        transcript: &mut DefaultTranscript<F>
    ) -> (Vec<FieldElement<F>>, [u8; 32])
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    let eval = Polynomial::evaluate_offset_fft::<F>(
        polynomial,
        1, 
        Some(domain_size),
        offset
    ).unwrap();
    
    let merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&eval);
    transcript.append_bytes(&merkle_tree.root);
    (eval, merkle_tree.root)
}

fn sample_queries<F>(
        num_queries: usize,
        domain_size: usize,
        transcript: &mut DefaultTranscript<F>
    ) -> Vec<u64> 
    where 
        F: IsField,
        FieldElement<F>: AsBytes + ByteConversion {

    let mut query_indices = Vec::<u64>::with_capacity(num_queries);

    for _ in 0..num_queries {
        let (_, query_index) = U256::from_bytes_be(
            &transcript.sample()
        ).unwrap()
        .div_rem(
            &U256::from(domain_size as u64)
        );
        query_indices.push(query_index.limbs[3]);
    }
    query_indices
}