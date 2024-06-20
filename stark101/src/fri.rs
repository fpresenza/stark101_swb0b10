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

pub struct Query<F: IsField> {
    index: usize,
    eval: FieldElement<F>,
    sym_eval: FieldElement<F>,
    // proof: Proof<[u8; 32]>,
    // sym_proof: Proof<[u8; 32]>,
}

impl<F: IsField> Query<F> {
    fn new(
        index: usize,
        eval: FieldElement<F>,
        sym_eval: FieldElement<F>,
    ) -> Self {
        Self {
            index,
            eval,
            sym_eval
        }
    }
}

pub struct FriLayer<F: IsField> {
    root: [u8; 32],
    queries: Vec<Query<F>>,
}

impl<F: IsField> FriLayer<F> {
    fn new(
        root: [u8; 32],
        queries: Vec<Query<F>>,
    ) -> Self {
        Self {
            root,
            queries
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
    println!(
        "Layer 0: appending root of composition polynomial (degree {:?}) to transcript.",
         polynomial.degree()
    );

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
            vec![Query::<F>::new(1, FieldElement::one(), FieldElement::one())]
        )
    );

    // recursive foldings
    for layer in 1..number_of_layers {
        let beta = transcript.sample_field_element();
    
        polynomial = poly::fold_polynomial(&polynomial, &beta);
        domain_size /= 2;
        offset = offset.square();

        let (eval, root) = commit(&polynomial, domain_size, &offset, transcript);
        println!(
            "Layer {:?}: appending root of folded polynomial (degree {:?}) to transcript.",
            layer,
            polynomial.degree()
        );

        fri_layers.push(
            FriLayer::<F>::new(
                root,
                vec![Query::<F>::new(1, FieldElement::one(), FieldElement::one())]
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