use lambdaworks_math::traits::{AsBytes, ByteConversion};
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
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

pub struct FriLayer {
    merkle_root: [u8; 32],
}

impl FriLayer {
    fn new(
        merkle_root: [u8; 32],
    ) -> Self {
        Self {
            merkle_root,
        }
    }
}

pub fn commit_and_fold<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        transcript: &mut DefaultTranscript<F>
    ) -> (Vec<FriLayer>, FieldElement<F>)
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion,
        <F as IsField>::BaseType: Sync + Send,
    {
        let mut curr_poly = polynomial.clone();
        let mut curr_offset = offset.clone();
        let number_of_layers = (usize::BITS - curr_poly.degree().leading_zeros() + 1) as usize;
        let mut fri_layers = Vec::<FriLayer>::with_capacity(number_of_layers);

        // before first folding
        let curr_eval = Polynomial::evaluate_offset_fft::<F>(
            &curr_poly,
            1,
            Some(domain_size),
            &curr_offset
        ).unwrap();
        let merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&curr_eval);
        transcript.append_bytes(&merkle_tree.root);
        fri_layers.push(
            FriLayer::new(
                merkle_tree.root,
            )
        );
    
        // recursive foldings
        for _ in 1..number_of_layers {
            let beta = transcript.sample_field_element();
    
            curr_poly = poly::fold_polynomial(&curr_poly, &beta);
            domain_size /= 2;
            curr_offset = curr_offset.square();
    
            let curr_eval = Polynomial::evaluate_offset_fft::<F>(
                &curr_poly,
                1, 
                Some(domain_size),
                &curr_offset
            ).unwrap();
    
            let curr_merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&curr_eval);
            transcript.append_bytes(&curr_merkle_tree.root);
        }

        let constant_poly = curr_poly.coefficients.first().unwrap();
        transcript.append_bytes(&constant_poly.to_bytes_be());

        (fri_layers, constant_poly.clone())
}