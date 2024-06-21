use lambdaworks_math::unsigned_integer::element::U256;
use lambdaworks_math::field::{
    element::FieldElement,
    traits::{IsField, IsFFTField}
};
use lambdaworks_math::traits::{AsBytes, ByteConversion};
use lambdaworks_crypto::merkle_tree::{
    merkle::MerkleTree,
    backends::types::Keccak256Backend, 
    proof::Proof
};
use lambdaworks_crypto::fiat_shamir::{
    default_transcript::DefaultTranscript
};

use crate::fri::FriLayer;

#[derive(Clone)]
pub struct PublicInput<F: IsField> (
	pub U256,
	pub usize,
	pub usize,
	pub usize,
	pub FieldElement<F>,
	pub FieldElement<F>
);

#[derive(Clone)]
pub struct InclusionProof<F: IsField> (
    pub FieldElement<F>,
    pub Proof<[u8; 32]>
);

#[derive(Clone)]
pub struct StarkProof<F: IsField> (
	pub [u8; 32],
	pub Vec<InclusionProof<F>>,
	pub Vec<FriLayer<F>>
);

pub fn sample_queries<F>(
        num_queries: usize,
        domain_size: usize,
        transcript: &mut DefaultTranscript<F>
    ) -> Vec<usize> 
    where 
        F: IsField,
        FieldElement<F>: AsBytes + ByteConversion {

        (0..num_queries)
        .map(|_| {
            let query_index = U256::from_bytes_be(&transcript.sample()).unwrap();
            let(_, query_index) = query_index.div_rem(&U256::from(domain_size as u64));
            query_index.limbs[3] as usize
        })
        .collect::<Vec<usize>>()
}

pub fn generate_inclusion_proofs<F>(
        indices: &[usize],
        poly_eval: &[FieldElement<F>],
        poly_tree: &MerkleTree<Keccak256Backend<F>>,
    ) -> Vec<InclusionProof<F>> 
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    indices
        .iter()
        .map(|i| {
        	InclusionProof(poly_eval[*i].to_owned(), poly_tree.get_proof_by_pos(*i).unwrap())
        })
        .collect()
}

pub fn verify_inclusion_proofs<F>(
        indices: &[usize],
        proofs: &Vec<InclusionProof<F>>,
        root: [u8; 32],
    ) -> bool
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    indices
        .iter()
        .zip(proofs)
        .map(|(index, InclusionProof(eval, proof))| {
            proof.verify::<Keccak256Backend<F>>(
                &root,
                *index,
                eval
            )
        }).all(|valid| valid)
}