use lambdaworks_math::unsigned_integer::element::U256;
use lambdaworks_math::field::{
	traits::IsField,
    element::FieldElement
};

use crate::fri::{InclusionProof, FriLayer};

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
pub struct StarkProof<F: IsField> (
	pub [u8; 32],
	pub Vec<[InclusionProof<F>; 3]>,
	pub Vec<FriLayer<F>>
);

