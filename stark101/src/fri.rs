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

use crate::poly;

pub struct FriLayer<F: IsField> {
    merkle_root: [u8; 32],
    query_value: FieldElement<F>,
    // query_merkle_proof,
    query_sym_value: FieldElement<F>,
    // query_sym_merkle_proof,
}

impl<F: IsField> FriLayer<F> {
    fn new(
        merkle_root: [u8; 32],
        query_value: FieldElement<F>,
        query_sym_value: FieldElement<F>,
    ) -> Self {
        Self {
            merkle_root,
            query_value,
            query_sym_value
        }
    }
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        transcript: &mut DefaultTranscript<F>
    ) -> (Vec<FriLayer<F>>, FieldElement<F>)
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion,
        <F as IsField>::BaseType: Sync + Send,
    {
        let mut polynomial = polynomial.clone();
        let mut offset = offset.clone();
        let number_of_layers = (usize::BITS - polynomial.degree().leading_zeros() + 1) as usize;
        let mut fri_layers = Vec::<FriLayer<F>>::with_capacity(number_of_layers);

        // before first folding
        let root = commit(&polynomial, domain_size, &offset, transcript);
        fri_layers.push(
            FriLayer::<F>::new(
                root,
                FieldElement::<F>::from(0_u64),
                FieldElement::<F>::from(0_u64),
            )
        );
    
        // recursive foldings
        for _ in 1..number_of_layers {
            let beta = transcript.sample_field_element();
    
            polynomial = poly::fold_polynomial(&polynomial, &beta);
            domain_size /= 2;
            offset = offset.square();

            let root = commit(&polynomial, domain_size, &offset, transcript);
            fri_layers.push(
                FriLayer::<F>::new(
                    root,
                    FieldElement::<F>::from(0_u64),
                    FieldElement::<F>::from(0_u64),
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
    ) -> [u8; 32]
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion,
        <F as IsField>::BaseType: Sync + Send,
    {
        let eval = Polynomial::evaluate_offset_fft::<F>(
            &polynomial,
            1, 
            Some(domain_size),
            &offset
        ).unwrap();
        
        let merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&eval);
        transcript.append_bytes(&merkle_tree.root);
        merkle_tree.root
}