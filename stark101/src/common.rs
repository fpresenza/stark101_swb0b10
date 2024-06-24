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

use crate::fri::FriCommitment;

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
pub struct VectorCommitment<F: IsField> {
	pub root: [u8; 32],
	pub inclusion_proofs: Vec<InclusionProof<F>>
}

#[derive(Clone)]
pub struct StarkProof<F: IsField> {
	pub trace_commitment: VectorCommitment<F>,
	pub composition_commitment: FriCommitment<F>
}

impl<F> VectorCommitment<F>
    where
        F: IsField + IsFFTField,
        FieldElement<F>: AsBytes + ByteConversion + Sync + Send {

    pub fn new_from_tree(tree: &MerkleTree<Keccak256Backend<F>>) -> Self {
        Self {
            root: tree.root,
            inclusion_proofs: vec![],
        }
    }

    pub fn generate_inclusion_proofs(
        &mut self,
        indices: &[usize],
        poly_eval: &[FieldElement<F>],
        poly_tree: &MerkleTree<Keccak256Backend<F>>,
    ) {

    self.inclusion_proofs.extend(
        indices
            .iter()
            .map(|i| {
                InclusionProof(poly_eval[*i].to_owned(), poly_tree.get_proof_by_pos(*i).unwrap())
            })
            .collect::<Vec<InclusionProof<F>>>()
        );
    }

    pub fn verify_inclusion_proofs(
            &self,
            indices: &[usize],
        ) -> bool {
    
        indices
            .iter()
            .zip(&self.inclusion_proofs)
            .map(|(index, InclusionProof(eval, proof))| {
                proof.verify::<Keccak256Backend<F>>(
                    &self.root,
                    *index,
                    eval
                )
            }).all(|valid| valid)
    }
}

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