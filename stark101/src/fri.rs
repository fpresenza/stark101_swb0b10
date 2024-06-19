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

pub fn commit<F>(
        polynomial: &Polynomial<FieldElement<F>>,
        mut domain_size: usize,
        offset: &FieldElement<F>,
        transcript: &mut DefaultTranscript<F>
    ) 
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion,
        <F as IsField>::BaseType: Sync + Send,
    {
        let mut curr_poly = polynomial.clone();
        // let mut domain_size = domain_size;
        let mut curr_offset = offset.clone();
    
        // before folding
        let eval = Polynomial::evaluate_offset_fft::<F>(
            &curr_poly,
            1,
            Some(domain_size),
            &curr_offset
        ).unwrap();
        let merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&eval);
        transcript.append_bytes(&merkle_tree.root);
    
        // recursive foldings
        let n_of_folds = usize::BITS - curr_poly.degree().leading_zeros();
        for _ in 1..=n_of_folds {
            let beta = transcript.sample_field_element();
    
            curr_poly = poly::fold_polynomial(&curr_poly, &beta);
            domain_size /= 2;
            curr_offset = curr_offset.square();
    
            let eval = Polynomial::evaluate_offset_fft::<F>(
                &curr_poly,
                1, 
                Some(domain_size),
                &curr_offset
            ).unwrap();
    
            let curr_merkle_tree = MerkleTree::<Keccak256Backend<F>>::build(&eval);
            transcript.append_bytes(&curr_merkle_tree.root);
        }
        println!("{:?}", curr_poly);
}